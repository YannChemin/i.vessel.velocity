"""Tests of i.vessel.velocity on synthetic Sentinel-2-like scenes."""

import importlib.machinery
import importlib.util
import math
import os
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pytest

import grass.script as gs
from grass.exceptions import CalledModuleError
from grass.experimental import TemporaryMapsetSession
from grass.script import array as garray
from grass.tools import Tools

_path = Path(__file__).resolve().parents[1] / "i.vessel.velocity.py"
_loader = importlib.machinery.SourceFileLoader("i_vessel_velocity", str(_path))
_spec = importlib.util.spec_from_loader("i_vessel_velocity", _loader)
ivv = importlib.util.module_from_spec(_spec)
_loader.exec_module(ivv)

# Sentinel-2 MSI acquisition times after B02, seconds.
LAGS = {"blue": 0.0, "green": 0.527, "red": 1.005, "nir": 0.264}
N, RES = 120, 10.0
KN = 1852.0 / 3600.0


def scene(ds, vessels, noise=0.003, seed=1):
    """Write blue, green, red and nir bands with moving hulls and wakes.

    :param vessels: list of (x, y, speed_kn, heading_deg, length_m)
    """
    gs.run_command("g.region", n=N * RES, s=0, w=0, e=N * RES, res=RES, env=ds.env)
    rng = np.random.default_rng(seed)
    yy, xx = np.mgrid[0:N, 0:N]
    for band, t in LAGS.items():
        a = np.full((N, N), 0.03) + rng.normal(scale=noise, size=(N, N))
        for x, y, kn, hdg, length in vessels:
            ux, uy = math.sin(math.radians(hdg)), math.cos(math.radians(hdg))
            v = kn * KN
            cx = (x + v * ux * t) / RES - 0.5
            cy = (N * RES - (y + v * uy * t)) / RES - 0.5
            ex, ey = xx - cx, -(yy - cy)
            along = ex * ux + ey * uy
            across = -ex * uy + ey * ux
            half = length / RES / 2
            a += 0.25 * np.exp(-0.5 * ((along / (half * 0.9)) ** 2 + (across / 0.5) ** 2))
            if kn > 0:
                a += np.where((along < -half - 1) & (along > -4 * half - 4)
                              & (np.abs(across) < 1), 0.03, 0.0)
        arr = garray.array(dtype=np.float32, env=ds.env)
        arr[...] = a
        arr.write(band, overwrite=True)


def points(ds, rows):
    """Vessel points with length and hull-axis columns."""
    txt = "\n".join(f"{x}|{y}|{i + 1}|{length}|{hdg % 180}"
                    for i, (x, y, _kn, hdg, length) in enumerate(rows))
    gs.write_command("v.in.ascii", input="-", output="pts", separator="pipe", x=1, y=2,
                     cat=3, columns="x double precision, y double precision, cat integer, "
                     "length double precision, orient double precision",
                     stdin=txt, env=ds.env)


def run(ds, **kw):
    ds.tools.i_vessel_velocity(input="pts", blue="blue", green="green", red="red",
                               nir="nir", output="out", length_column="length",
                               orient_column="orient", **kw)
    sel = gs.vector_db_select("out", env=ds.env)
    cols = sel["columns"]
    return [dict(zip(cols, v)) for _, v in sorted(sel["values"].items())]


@pytest.fixture
def ds(tmp_path):
    """Tools in a temporary mapset of an XY project."""
    gs.create_project(tmp_path / "xy")
    with (
        gs.setup.init(tmp_path / "xy", env=os.environ.copy()) as session,
        TemporaryMapsetSession(env=session.env) as mapset,
        Tools(session=mapset) as tools,
    ):
        yield SimpleNamespace(tools=tools, env=mapset.env)


def test_fit_velocity_exact():
    """A target moving at 5 m/s east and 2 m/s south is fitted exactly."""
    t = [0.0, 0.264, 0.527, 1.005]
    vx, vy, se = ivv.fit_velocity(t, [100 + 5 * s for s in t], [50 - 2 * s for s in t])
    assert vx == pytest.approx(5.0)
    assert vy == pytest.approx(-2.0)
    assert se == pytest.approx(0.0, abs=1e-9)


def test_moving_and_stationary(ds):
    """Speed, heading, the moving flag and the wake side are recovered."""
    rows = [(400.0, 700.0, 12.0, 60.0, 30.0), (800.0, 300.0, 0.0, 120.0, 30.0)]
    scene(ds, rows)
    points(ds, rows)
    mov, sta = run(ds)
    assert mov["kin_status"] == "ok"
    assert float(mov["speed_kn"]) == pytest.approx(12.0, abs=1.5)
    assert float(mov["heading_deg"]) == pytest.approx(60.0, abs=5.0)
    assert int(mov["moving"]) == 1
    assert float(mov["axis_diff_deg"]) < 10
    assert float(mov["wake_behind"]) > float(mov["wake_ahead"])
    assert float(sta["speed_kn"]) < 1.5
    assert int(sta["moving"]) == 0
    assert sta["wake_behind"] == ""


def test_bias_is_subtracted(ds):
    """A calibration bias equal to the true velocity leaves a still vessel."""
    rows = [(500.0, 500.0, 10.0, 0.0, 30.0)]
    scene(ds, rows)
    points(ds, rows)
    (r,) = run(ds, bias=[0.0, 10.0 * KN])
    assert float(r["speed_kn"]) < 1.5
    assert int(r["moving"]) == 0


def test_point_near_edge(ds):
    """A vessel whose window leaves the region is flagged, not measured."""
    rows = [(30.0, 600.0, 10.0, 90.0, 30.0)]
    scene(ds, rows)
    points(ds, rows)
    (r,) = run(ds)
    assert r["kin_status"] == "edge"
    assert r["speed_kn"] == ""


def test_wrong_number_of_lags(ds):
    """Three time lags for four bands fail loudly."""
    rows = [(500.0, 500.0, 10.0, 0.0, 30.0)]
    scene(ds, rows)
    points(ds, rows)
    with pytest.raises(CalledModuleError):
        run(ds, time_lags=[0, 0.5, 1.0])


def test_input_left_unchanged(ds):
    """The input map gets no new columns and the region is kept."""
    rows = [(500.0, 500.0, 10.0, 0.0, 30.0)]
    scene(ds, rows)
    points(ds, rows)
    run(ds)
    assert "speed_kn" not in gs.vector_columns("pts", env=ds.env)
    assert gs.region(env=ds.env)["rows"] == N


@pytest.mark.parametrize("kn,hdg", [(6.0, 0.0), (12.0, 135.0), (20.0, 250.0)])
def test_speed_and_heading_sweep(ds, kn, hdg):
    """Speed within 15 % or 1.5 knots, the documented precision, and
    heading within 5 degrees, wake included."""
    rows = [(600.0, 600.0, kn, hdg, 35.0)]
    scene(ds, rows, seed=int(kn))
    points(ds, rows)
    (r,) = run(ds)
    assert float(r["speed_kn"]) == pytest.approx(kn, abs=max(0.15 * kn, 1.5))
    d = abs(float(r["heading_deg"]) - hdg) % 360
    assert min(d, 360 - d) < 5
    assert int(r["moving"]) == 1
