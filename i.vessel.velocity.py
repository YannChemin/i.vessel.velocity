#!/usr/bin/env python3
############################################################################
#
# MODULE:       i.vessel.velocity
# AUTHOR(S):    Yann Chemin
# PURPOSE:      Estimates the speed and heading of vessels from the
#               displacement between the bands of a push-broom image
#               (Sentinel-2 MSI by default): each band is acquired at a
#               different time, so a moving vessel shifts between bands.
#
# SPDX-License-Identifier: Unlicense
#
############################################################################

# %Module
# % description: Estimates vessel speed and heading from the inter-band displacement of a multispectral push-broom image (Sentinel-2 MSI by default).
# % keyword: imagery
# % keyword: vector
# % keyword: vessel
# % keyword: velocity
# % keyword: Sentinel
# % keyword: maritime
# %end

# %option G_OPT_V_INPUT
# % description: Vessel points (barycentres of the detections)
# %end

# %option G_OPT_V_FIELD
# %end

# %option G_OPT_R_INPUT
# % key: blue
# % description: Blue band (Sentinel-2 B02)
# % guisection: Bands
# %end

# %option G_OPT_R_INPUT
# % key: green
# % description: Green band (Sentinel-2 B03)
# % guisection: Bands
# %end

# %option G_OPT_R_INPUT
# % key: red
# % description: Red band (Sentinel-2 B04)
# % guisection: Bands
# %end

# %option G_OPT_R_INPUT
# % key: nir
# % description: Near-infrared band (Sentinel-2 B08)
# % guisection: Bands
# %end

# %option G_OPT_V_OUTPUT
# % description: Vessel points with speed, heading and quality columns
# %end

# %option
# % key: time_lags
# % type: double
# % required: no
# % multiple: yes
# % answer: 0,0.527,1.005,0.264
# % description: Acquisition time of the blue, green, red and NIR bands, in seconds after the earliest (default: Sentinel-2 MSI, ESA)
# % guisection: Sensor
# %end

# %option G_OPT_DB_COLUMN
# % key: length_column
# % required: no
# % description: Column with the vessel length in metres (sets the measurement window; default 20 m)
# % guisection: Attributes
# %end

# %option G_OPT_DB_COLUMN
# % key: orient_column
# % required: no
# % description: Column with the hull axis as a bearing from grid north, 0-180 degrees (for the axis check)
# % guisection: Attributes
# %end

# %option
# % key: bias
# % type: double
# % required: no
# % multiple: yes
# % key_desc: vx,vy
# % description: Systematic apparent velocity to subtract, east and north components in m/s (calibration of the inter-band registration)
# % guisection: Sensor
# %end

# %option
# % key: method
# % type: string
# % required: no
# % options: barycentre,outline,both
# % answer: both
# % description: Position of the vessel in each band
# % descriptions: barycentre;intensity-weighted barycentre of the hull;outline;shift of the saturated hull silhouette between bands (insensitive to brightness within the hull);both;both, combined by inverse variance, with an agreement flag
# % guisection: Sensor
# %end

# %option
# % key: sigma
# % type: double
# % required: no
# % answer: 3
# % description: A vessel is moving when its speed exceeds this many standard errors
# % guisection: Sensor
# %end

# %option
# % key: min_speed
# % type: double
# % required: no
# % answer: 2
# % description: A vessel is moving only above this speed, in knots
# % guisection: Sensor
# %end

import math
import sys

import grass.script as gs

MS_TO_KN = 3600.0 / 1852.0
# Power of the brightness above background used as barycentre weight.
HULL_POWER = 2
BAND_KEYS = ("blue", "green", "red", "nir")
NEW_COLUMNS = [
    ("n_bands", "INTEGER"),
    ("snr_min", "DOUBLE PRECISION"),
    ("vx_ms", "DOUBLE PRECISION"),
    ("vy_ms", "DOUBLE PRECISION"),
    ("speed_kn", "DOUBLE PRECISION"),
    ("speed_se_kn", "DOUBLE PRECISION"),
    ("heading_deg", "DOUBLE PRECISION"),
    ("axis_diff_deg", "DOUBLE PRECISION"),
    ("moving", "INTEGER"),
    ("bary_vx", "DOUBLE PRECISION"),
    ("bary_vy", "DOUBLE PRECISION"),
    ("bary_se_kn", "DOUBLE PRECISION"),
    ("outl_vx", "DOUBLE PRECISION"),
    ("outl_vy", "DOUBLE PRECISION"),
    ("outl_se_kn", "DOUBLE PRECISION"),
    ("agree", "INTEGER"),
    ("wake_behind", "DOUBLE PRECISION"),
    ("wake_ahead", "DOUBLE PRECISION"),
    ("kin_status", "VARCHAR(12)"),
]


def barycentres(chips, yy, xx, r0, c0, core, ring_in, ring_out):
    """Intensity-weighted barycentre (row, col) of the target in each band.

    :param chips: dict band key -> 2-D array around the target
    :param yy: row index grid of the chips
    :param xx: column index grid of the chips
    :param r0: row of the target barycentre
    :param c0: column of the target barycentre
    :param core: radius in pixels of the disc that holds the hull
    :param ring_in: inner radius of the background ring
    :param ring_out: outer radius of the background ring
    :return: (positions, signal-to-noise ratios, background-subtracted chips)
    """
    import numpy as np
    from scipy import ndimage

    dist = np.hypot(yy - r0, xx - c0)
    ring = (dist > ring_in) & (dist < ring_out)
    disc = dist <= core
    pos, snr, sub = {}, {}, {}
    for key, chip in chips.items():
        v = chip[ring & np.isfinite(chip)]
        if v.size < 20 or not np.isfinite(chip[disc]).all():
            continue
        bg = float(np.median(v))
        sig = 1.4826 * float(np.median(np.abs(v - bg))) or 1e-6
        # Weights grow as a power of the brightness above background, so the
        # bright hull outweighs its wake: the foam of a wake stays where it
        # was made while the vessel moves on, and would pull the barycentre
        # towards zero motion. A smooth weight, unlike a threshold, does not
        # jump as the hull moves within a pixel.
        wts = np.where(disc, np.clip(chip - bg - 2 * sig, 0, None), 0.0) ** HULL_POWER
        if wts.sum() <= 0:
            continue
        # Keep the bright blob connected to the brightest pixel of the disc,
        # not a neighbouring vessel inside the window.
        lab, _ = ndimage.label(wts > 0)
        pk = np.unravel_index(np.argmax(wts), wts.shape)
        wts = np.where(lab == lab[pk], wts, 0.0)
        pos[key] = (float((yy * wts).sum() / wts.sum()), float((xx * wts).sum() / wts.sum()))
        snr[key] = float((chip[disc].max() - bg) / sig)
        sub[key] = chip - bg
    return pos, snr, sub


def fit_velocity(times, xs, ys):
    """Least-squares velocity of positions against acquisition times.

    :param times: acquisition times, seconds
    :param xs: easting of the target in each band, metres
    :param ys: northing of the target in each band, metres
    :return: (vx, vy, standard error of the speed), m/s
    """
    import numpy as np

    t = np.asarray(times, float)
    A = np.vstack([t, np.ones_like(t)]).T
    cx, *_ = np.linalg.lstsq(A, np.asarray(xs, float), rcond=None)
    cy, *_ = np.linalg.lstsq(A, np.asarray(ys, float), rcond=None)
    resid = np.concatenate([np.asarray(xs) - A @ cx, np.asarray(ys) - A @ cy])
    dof = max(2 * len(t) - 4, 1)
    s_pos = math.sqrt(float((resid**2).sum()) / dof)
    se = s_pos / math.sqrt(float(((t - t.mean()) ** 2).sum()))
    return float(cx[0]), float(cy[0]), se


def outline_image(chip, bg, sig, dist, core):
    """Saturated silhouette of the hull, tapered to its footprint.

    The absolute contrast makes the silhouette the same whether the hull is
    brighter or darker than the sea. Contrast below a floor counts as 0 and
    above a ceiling as 1, so the inside of the hull is flat whatever its
    brightness pattern (superstructure, deck, shadow), the faint wake falls
    below the floor, and only the edge pixels, through their partial
    coverage, carry the position.
    """
    from scipy import ndimage

    import numpy as np

    # Thresholds chosen on synthetic hulls of 30-90 m: smoothing 1 pixel,
    # floor 15 % and ceiling 50 % of the peak contrast.
    a = ndimage.gaussian_filter(np.abs(np.nan_to_num(chip - bg)), 1.0)
    peak = float(a[dist <= core].max()) if (dist <= core).any() else 0.0
    lo = max(3 * sig, 0.15 * peak)
    hi = max(1.5 * lo, 0.5 * peak)
    s = np.clip((a - lo) / (hi - lo), 0, 1) if hi > lo else np.zeros_like(a)
    # Flat over the hull and its motion, so the window does not pull the
    # shift towards zero; it only cuts neighbours further out.
    return s * np.exp(-((dist / (core + 4)) ** 8))


def register(ref, mov, mask, iters=6):
    """Sub-pixel shift (rows, cols) with mov(p) = a ref(p - d) + c.

    Gauss-Newton on the linearised model over the pixels of mask; the gain
    a and offset c absorb the different contrast of each band.
    """
    from scipy import ndimage

    import numpy as np

    gr, gc = np.gradient(ref)
    X = np.column_stack([ref[mask], np.ones(int(mask.sum())), -gc[mask], -gr[mask]])
    d = np.zeros(2)
    for _ in range(iters):
        m = ndimage.shift(mov, (-d[0], -d[1]), order=3, mode="nearest")
        coef, *_ = np.linalg.lstsq(X, m[mask], rcond=None)
        a = coef[0] if abs(coef[0]) > 1e-12 else 1.0
        d += np.array([coef[3], coef[2]]) / a
        if math.hypot(*d) > 3:
            break
    return d


def outline_shifts(chips, dist, core, ref="nir"):
    """Shift of each band's hull outline relative to the reference band, pixels."""
    import numpy as np

    ring = (dist > core + 4) & (dist < core + 8)
    imgs = {}
    for k, c in chips.items():
        v = c[ring & np.isfinite(c)]
        if v.size < 20:
            continue
        bg = float(np.median(v))
        sig = 1.4826 * float(np.median(np.abs(v - bg))) or 1e-6
        imgs[k] = outline_image(c, bg, sig, dist, core)
    if ref not in imgs:
        return {}
    mask = dist <= core + 5
    return {k: (np.zeros(2) if k == ref else register(imgs[ref], im, mask))
            for k, im in imgs.items()}


def wake_strips(chip, r0, c0, ux, uy, dist_px):
    """Mean brightness on 3-pixel strips behind and ahead of the target.

    :param chip: background-subtracted NIR chip
    :param r0: row of the target in the chip
    :param c0: column of the target in the chip
    :param ux: east component of the unit velocity
    :param uy: north component of the unit velocity
    :param dist_px: distances along the track to sample, pixels
    :return: (behind, ahead)
    """
    import numpy as np
    from scipy import ndimage

    off = np.array([-1.0, 0.0, 1.0])
    out = []
    for sign in (-1, 1):
        cols = c0 + sign * dist_px[:, None] * ux + off[None, :] * uy
        rows = r0 - sign * dist_px[:, None] * uy + off[None, :] * ux
        vals = ndimage.map_coordinates(chip, [rows.ravel(), cols.ravel()], order=1,
                                       mode="nearest")
        out.append(float(np.mean(vals)))
    return out[0], out[1]


def measure(bands, times, x, y, length, orient, bias, sigma, min_kn, env, method="both"):
    """Velocity of one vessel; returns a dict of the new columns."""
    import numpy as np
    from grass.script import array as garray

    reg = gs.region(env=env)
    res = float(reg["ewres"])
    length = length if length and length == length else 20.0
    wake_l = max(length, 20.0)
    core = int(math.ceil(length / res / 2)) + 3
    half = int(max(core + 8, math.ceil(3.5 * wake_l / res)))
    # Chip of (2 half + 1) cells centred on the cell that holds the point.
    col = math.floor((x - float(reg["w"])) / res)
    row = math.floor((float(reg["n"]) - y) / res)
    if (row - half < 0 or col - half < 0 or row + half >= int(reg["rows"])
            or col + half >= int(reg["cols"])):
        return {"kin_status": "edge"}
    n = float(reg["n"]) - (row - half) * res
    w = float(reg["w"]) + (col - half) * res
    size = (2 * half + 1) * res
    with gs.RegionManager(n=n, s=n - size, w=w, e=w + size, res=res, env=env.copy()) as rm:
        chips = {k: np.asarray(garray.array(m, env=rm.env), dtype=float) for k, m in bands.items()}
    for c in chips.values():
        c[~np.isfinite(c)] = np.nan
    yy, xx = np.mgrid[0:2 * half + 1, 0:2 * half + 1]
    r0 = (n - y) / res - 0.5
    c0 = (x - w) / res - 0.5
    pos, snr, sub = barycentres(chips, yy, xx, r0, c0, core, core + 4, core + 8)
    out = {"n_bands": len(pos), "snr_min": round(min(snr.values()), 1) if snr else None}
    est = {}
    if method in ("barycentre", "both") and len(pos) >= 3:
        keys = list(pos)
        xs = [w + (pos[k][1] + 0.5) * res for k in keys]
        ys = [n - (pos[k][0] + 0.5) * res for k in keys]
        est["bary"] = fit_velocity([times[k] for k in keys], xs, ys)
    if method in ("outline", "both"):
        dist = np.hypot(yy - r0, xx - c0)
        sh = outline_shifts(chips, dist, core)
        if len(sh) >= 3:
            keys = list(sh)
            est["outl"] = fit_velocity([times[k] for k in keys], [sh[k][1] * res for k in keys],
                                       [-sh[k][0] * res for k in keys])
    for tag, (evx, evy, ese) in est.items():
        out.update({f"{tag}_vx": round(evx - bias[0], 3), f"{tag}_vy": round(evy - bias[1], 3),
                    f"{tag}_se_kn": round(ese * MS_TO_KN, 2)})
    if not est:
        out["kin_status"] = "few_bands"
        return out
    if len(est) == 2:
        (bx, by, bs), (ox, oy, os_) = est["bary"], est["outl"]
        # The two positions share the band registration but not the hull
        # appearance; they agree when their difference is within two
        # combined standard errors (with a floor of 1 knot).
        tol = 2 * max(math.hypot(bs, os_), 1.0 / MS_TO_KN)
        out["agree"] = int(math.hypot(bx - ox, by - oy) <= tol)
        # Inverse-variance mean; a floor keeps a perfect four-point fit from
        # taking all the weight.
        wb, wo = 1 / max(bs, 0.15) ** 2, 1 / max(os_, 0.15) ** 2
        vx = (wb * bx + wo * ox) / (wb + wo)
        vy = (wb * by + wo * oy) / (wb + wo)
        se = 1 / math.sqrt(wb + wo)
    else:
        vx, vy, se = next(iter(est.values()))
    vx -= bias[0]
    vy -= bias[1]
    speed = math.hypot(vx, vy)
    heading = math.degrees(math.atan2(vx, vy)) % 360
    out.update({"vx_ms": round(vx, 3), "vy_ms": round(vy, 3),
                "speed_kn": round(speed * MS_TO_KN, 2), "speed_se_kn": round(se * MS_TO_KN, 2),
                "heading_deg": round(heading, 1)})
    if orient is not None and orient == orient:
        d = abs((heading % 180) - float(orient)) % 180
        out["axis_diff_deg"] = round(min(d, 180 - d), 1)
    moving = speed > max(sigma * se, min_kn / MS_TO_KN)
    out["moving"] = int(moving)
    if moving and "nir" in sub:
        dist_px = np.linspace(wake_l, 3 * wake_l, 25) / res
        behind, ahead = wake_strips(sub["nir"], r0, c0, vx / speed, vy / speed, dist_px)
        out["wake_behind"], out["wake_ahead"] = round(behind, 5), round(ahead, 5)
    out["kin_status"] = "ok"
    return out


def sql_value(v):
    if v is None or (isinstance(v, float) and v != v):
        return "NULL"
    if isinstance(v, str):
        return "'{}'".format(v.replace("'", "''"))
    return repr(v)


def main():
    options, flags = gs.parser()
    import os
    import tempfile

    import numpy as np  # noqa: F401  (fails early when numpy is missing)

    try:
        import scipy  # noqa: F401
    except ImportError:
        gs.fatal(_("The scipy Python package is required."))

    bands = {k: options[k] for k in BAND_KEYS}
    lags = [float(v) for v in options["time_lags"].split(",")]
    if len(lags) != 4:
        gs.fatal(_("<time_lags> needs four values (blue, green, red, nir), got {}.").format(len(lags)))
    times = dict(zip(BAND_KEYS, lags))
    bias = [0.0, 0.0]
    if options["bias"]:
        bias = [float(v) for v in options["bias"].split(",")]
        if len(bias) != 2:
            gs.fatal(_("<bias> needs two values, vx and vy in m/s."))
    sigma, min_kn = float(options["sigma"]), float(options["min_speed"])
    inp, layer, out = options["input"], options["layer"], options["output"]
    lcol, ocol = options["length_column"], options["orient_column"]

    if gs.vector_info_topo(inp)["points"] == 0:
        gs.fatal(_("Vector map <{}> has no points.").format(inp))
    cols = gs.vector_columns(inp, layer=layer)
    for c in (lcol, ocol):
        if c and c not in cols:
            gs.fatal(_("Column <{}> not found in vector map <{}>.").format(c, inp))
    for c, _t in NEW_COLUMNS:
        if c in cols:
            gs.fatal(_("Vector map <{}> already has a column <{}>.").format(inp, c))
    key = gs.vector_db(inp)[int(layer)]["key"]

    coords = {}
    for line in gs.read_command("v.out.ascii", input=inp, layer=layer, format="point",
                                separator="pipe").splitlines():
        p = line.split("|")
        coords[int(p[2])] = (float(p[0]), float(p[1]))
    attrs = {}
    if lcol or ocol:
        sel = gs.vector_db_select(inp, layer=int(layer),
                                  columns=",".join([key] + [c for c in (lcol, ocol) if c]))
        names = sel["columns"]
        for vals in sel["values"].values():
            rec = dict(zip(names, vals))
            attrs[int(rec[key])] = rec

    def num(rec, c):
        try:
            return float(rec.get(c)) if c and rec.get(c) not in (None, "") else None
        except ValueError:
            return None

    gs.run_command("g.copy", vector=(inp, out), quiet=True)
    gs.run_command("v.db.addcolumn", map=out, layer=layer,
                   columns=",".join(f"{c} {t}" for c, t in NEW_COLUMNS), quiet=True)
    env = os.environ.copy()
    gs.message(_("Measuring {} vessels...").format(len(coords)))
    stmts = []
    n_moving = 0
    for i, (cat, (x, y)) in enumerate(sorted(coords.items())):
        gs.percent(i, len(coords), 2)
        rec = attrs.get(cat, {})
        r = measure(bands, times, x, y, num(rec, lcol), num(rec, ocol), bias, sigma, min_kn, env,
                    options["method"])
        n_moving += r.get("moving", 0) or 0
        sets = ", ".join(f"{c} = {sql_value(r.get(c))}" for c, _t in NEW_COLUMNS)
        stmts.append(f"UPDATE {out} SET {sets} WHERE {key} = {cat};")
    gs.percent(1, 1, 1)
    # The attribute table may carry the layer name; read it from the map.
    table = gs.vector_db(out)[int(layer)]["table"]
    stmts = [s.replace(f"UPDATE {out} ", f"UPDATE {table} ", 1) for s in stmts]
    with tempfile.NamedTemporaryFile("w", suffix=".sql", delete=False) as f:
        f.write("BEGIN;\n" + "\n".join(stmts) + "\nCOMMIT;\n")
        sql = f.name
    try:
        gs.run_command("db.execute", input=sql,
                       database=gs.vector_db(out)[int(layer)]["database"],
                       driver=gs.vector_db(out)[int(layer)]["driver"])
    finally:
        os.remove(sql)
    gs.vector_history(out)
    gs.message(_("{} of {} vessels moving.").format(n_moving, len(coords)))
    return 0


if __name__ == "__main__":
    sys.exit(main())
