# CS-hybrid model: SLR dune-toe implementation plan

Nov 12, 2026 · @Filipe

## 1. Scope

Both changes reduce to one new operation, a "base-plane shift", called once per step right after `morpho_update`. The plan adds two yaml flags, both off by default, so every existing run stays bit-identical.

| Flag | Change | What it does |
| --- | --- | --- |
| `slr_toe` | 1 | Dune toe and seaward flat of the hard structure rise with sea level (`slr * dt` per step). The toe slides back along the ramp face, the dike toe point slides up the dike face. |
| `dynamic_toe` | 2 | Berm point and maximum toe level, as you proposed. The toe erodes up the beach slope and then landward on a flat at `dfoot`; it progrades on the flat and then down the beach slope. The beach budget is exact, so the phantom beach volume is gone. |

What stays untouched: the logic of `morpho_update` (two small robustness edits, section 5 step 1), `deposition_dune`, vegetation, bar and nourishment routines, and the lee side behind the dike.

**Status on 13 November.** The sandbox holds the full implementation of both flags and was run on your test16d case (79,056 steps, three storms).

- **Blocking issue:** the `functions.py` you uploaded crashes on test16d at step 6243 with all flags off, and your GIF shows behaviour it does not have. It is not the version you run. I used a marked workaround to get a baseline; please send the current file (section 8, item 1).
- **Regression:** flags off, and `slr_toe` on with `slr = 0`, are bit-identical to that baseline on test16d.
- **Change 1:** 0.5 m of rise on test16d runs through dike exposure and detachment; toe and substrate stay consistent at every step and the SLR terms add no balance error.
- **Change 2:** the berm option runs test16d with and without SLR; the beach budget closes to about 1e-12 m³/m per step with no phantom correction. After each storm a flat 54 to 60 m wide is left at `dfoot` and the dune regrows over it.
- **Not yet shown:** a toe far down the beach slope. In test16d the dune never passes the berm point, and in the synthetic case the toe went only 0.15 m below `dfoot`.

Line numbers refer to the files you uploaded on 12 November 2026.

## 2. Current behaviour

The toe elevation already exists in two places that never talk to each other: `main.py` uses the constant `p.dfoot`, while `morpho_update` reads `dfoot = dune_z['ramp_toe']` (functions.py line 1463) every step. That is why the change is small.

- **Seaward flat = toe level.** `morpho_update` uses `s_hs + dfoot` as the dike crest elevation in 12 places (lines 1484, 1786, 1799, 1801, 2036, 2063, 2069, 2155, 2159, 3089, 3092, 3110). This is only true while `substrate_z[4] == substrate_z[5] == dune_z['ramp_toe']`. Raising the toe and the flat by the same amount keeps it true.
- **Dike face slope.** Ramp-against-dike formulas (e.g. lines 2050–2051, 2143–2144) assume the face slope is `betas_hs`. If only `substrate_z[4]` were raised with `substrate_x[4]` fixed, the face would get flatter and these formulas would break. Point 4 must slide along the face.
- **SLR today.** Only the Bruun sink exists: `qslr = slr * bactive` in `subaq`, subtracted from the beach at main.py line 1838. Neither conservation routine knows about `qslr`.
- **Phantom volume.** `vbeachtot = ramp_toe * (dclos + dfoot) + width * (dclos + db)`. When the toe moves by Δ at a fixed elevation, this changes by Δ·(dfoot − db) with no sediment flux. Lines 1833–1835 accept it (`vbeachtot_phantom`), and `check_mass_conservation` subtracts it afterwards (`phantom_correction`, mass\_conservation.py lines 376–396).
- **Lee foot on the seaward flat.** For a dune in front of the dike, the lee foot is also set to `dfoot` (functions.py lines 2930, 2989, 3011). It rests on the same plane, so it has to move with it.
- **Stored substrate.** `results['substrate_x']` and `['substrate_z']` are stored once. The post-run conservation check and your plot scripts assume they never change.

Pre-existing quirks I noticed and will not touch (they behave the same before and after):

- Some erosion branches set `dune_z['ramp_toe']` to a level on the dike (lines 1652, 1683). The base-plane shift skips those states.
- A few ramp branches mix absolute elevation and height above the toe (lines 2097–2098, 2283, 2339–2341).
- `np.interp` over `dune_x.values()` (lines 1529, 1577) includes the extra key `ramp_crest_temp` once an erosion branch has created it.

## 3. Design of change 1: toe and seaward flat follow sea level

The model keeps its fixed vertical datum; mean sea level `msl` becomes a state variable that rises by `slr * dt` per step, and the base plane `z_plane` (toe level = `substrate_z[4]` = `substrate_z[5]`) rises with it. The toe height above mean sea level, `h_beach = z_plane - msl`, therefore stays equal to `p.dfoot`, as in the original CS-model.

**New state variables in `main.py`**

| Variable | Meaning | Initial value |
| --- | --- | --- |
| `msl` | cumulative sea-level rise (m, model datum) | 0.0 |
| `z_plane` | toe level and seaward flat of the substrate (m, model datum) | `substrate_z[4]` |
| `h_beach` | `z_plane - msl`, toe height above mean sea level | `p.dfoot` |
| `m_face` | slope of the seaward dike face, from the initial substrate, never recomputed | `(z[3]-z[4]) / (x[4]-x[3])` |

**What the base-plane shift does** (function `shift_base_plane`, section 5 step 1)

1. Substrate: `z[4] = z[5] = z_new`, and `x[4] = x[3] + (z[3] - z_new) / m_face`. The dike is not moved; its toe point slides up the face, so the face slope stays `betas_hs`.
2. Seaward pivots resting on the plane are re-anchored, by state:

| State of the dune front | What happens when the plane rises by dz |
| --- | --- |
| Ramp (full or partial), toe seaward of the dike toe | Toe and `ramp_toe_temp` slide landward along the ramp face by dz / betas. Ramp slope is preserved. |
| Vertical step at the toe (`ramp_toe_temp` above `ramp_toe`) | Toe rises in place. If dz exceeds the step, the step disappears and the toe continues up the ramp face. |
| Scarp with empty ramp (`ramp_crest`, `ramp_toe_temp`, `ramp_toe` at the scarp foot) | All three rise in place and stay equal. |
| Toe at the dike toe (dike face exposed) | Toe follows substrate point 4 up the face. Pivots on the dike above the plane are not touched. |
| Toe drawn at plane level under the dike face (scarp landward of the dike toe) | Toe keeps its x and only follows the level. Found with test16d; without this case the toe level froze during the first storm and the run crashed 1700 steps later. |
| Toe not on the plane (erosion states that set `ramp_toe` on the dike) | Only the substrate changes. |

3. Lee foot: only if it rests on the seaward plane (dune in front of the dike). It slides along the lee slope to the new plane, or onto the dike face if it meets the face first. `back_ref` follows when it also sits on the seaward plane.
4. The function measures the dune volume before and after with `partial_volume` and returns the difference, `toe_relabel`. This is the slab at the base of the dune that is now below the toe level.

**Where the slab goes.** It stays physically under the dune, so it is not given to the beach wedge and the shoreline does not advance because of it. In the model's bookkeeping it replaces the bottom slab of the under-dune block, which leaves the control volume because the depth of closure rises with sea level. It is a loss term for the dune compartment and for the total in the conservation check, together with the Bruun sink `qslr * dt`. Decision 2 in section 8 covers the alternative.

**Water levels.** With a fixed datum, thresholds compare absolute levels, and the beach formulas use heights above mean sea level:

- absolute toe level (`dfoot_abs = z_plane`): erosion threshold and impact height, dry-beach test, the aeolian reduction with dune height;
- height above mean sea level (`h_beach`, with `swl - msl`): runup friction term X, `rprim`, dry-beach width `yr`, all beach-volume formulas.

`dune_z['shoreline']` becomes `msl` instead of 0, so the foreshore slope stays (toe − shoreline elevation) / width.

**Input SWL.** Sea-level rise is added in the simulation: the model uses `swl + msl`. This follows your answer that the SWL files do not contain the trend. A flag `swl_includes_slr` (default false) remains for forcing files that already include it.

**Order in the loop.** The shift runs right after `morpho_update` and before the vegetation and beach blocks, not at the very end. The beach update must see the new toe position, otherwise the small toe retreat (dz / betas) would be one step out of sync with the phantom term.

## 4. Design of change 2: berm point, maximum toe level, no phantom volume

This section follows your berm suggestion and replaces my earlier symmetric toe rule. The beach profile gets two sides: a flat at the maximum toe level, ending at a berm point, and the beach face from the berm point to the shoreline. The dune sits on that profile, and its toe is wherever its front meets it.

**4.1 Rules**

| Element | Rule |
| --- | --- |
| Maximum toe level | `z_max = msl + dfoot`. It moves only with sea level, never with erosion or accretion. |
| Berm point | `x_berm`, the seaward end of the flat at `z_max`. It starts at the initial toe. It does not move when the dune erodes or regrows. |
| Erosion, toe below `z_max` | The toe retreats up the beach slope (the line from the previous toe to the berm point) until it reaches the berm point. |
| Erosion, toe at `z_max` | The toe retreats landward on the flat. The flat between toe and berm point is exposed; the eroded sand goes to the beach face, so the shoreline advances. |
| Accretion on the flat | The toe advances at `z_max` until it reaches the berm point. No sand is needed to raise the bed, because the flat is already there. |
| Accretion past the berm point | The toe moves down the beach face (the line from the previous toe to the shoreline). Each step lowers the erosion threshold, so waves reach the toe more often. That is the hydrodynamic limit you describe. |
| Floor | `dfoot_min`, a lowest toe height above mean sea level, as a numerical safeguard only. |
| Beach-face slope limits | The face slope is `dfoot` / face width. Steeper than `beach_slope_max` (erosive conditions): the berm point retreats and the flat narrows, down to the toe. Flatter than `beach_slope_min` (accretive conditions): the berm point advances and the flat widens. In between, the berm point stays where it is and the shoreline moves. |

**Dune or berm first in erosive conditions.** The two are driven by different things, so the model does not have to choose. The berm point reacts to the beach volume: it retreats while the face is at `beach_slope_max` and sand keeps leaving the beach. The dune toe reacts to water level: it retreats when runup plus SWL exceeds the toe level. A slow sediment deficit therefore eats the flat first and reaches the dune only when the berm point arrives at the toe; a storm erodes the dune directly and, because the eroded sand flattens the face, usually leaves the berm point in place or pushes it seaward.

**Sea-level rise.** The flat and the berm point rise with the toe: `z_max = msl + dfoot`, the berm point keeps its x, and the beach face keeps its shape relative to mean sea level. That is the Bruun assumption (the active profile moves up with the sea), and the sand it takes is the existing sink `qslr = slr * bactive`, which makes the shoreline retreat. Once the face reaches `beach_slope_max`, the berm point retreats as well, so the profile moves up and landward. One caveat: `bactive` is a constant from the config, so it does not grow when the flat widens.

**4.2 How it is computed.** `morpho_update` still works on a flat plane at the current toe level and is not changed. After it, the shift from section 3 is called with the level at which the toe lies on the beach profile (`find_toe_level`, section 5 step 1): `z_max` if the dune front meets the flat, otherwise the root of a one-dimensional equation on the beach slope. This is a geometric condition solved every step, not an increment, so the 20 m toe jumps of the narrow-beach branches need no special handling.

Dike interaction always happens with the plane at `z_max`, because the berm point is seaward of the dike toe. The dike branches of `morpho_update` therefore see the same geometry as today.

**4.3 Volumes: your question.** Yes. If the dune moves along the beach profile, the dune volume can be measured each step above the profile that was there before the dune, and the sand below it stays beach. Two things implement that:

- **Dune volume.** New output `vdune_true`: volume above the beach base profile (flat at `z_max` up to the berm point, beach slope from there to the toe). `vdune` remains the volume above the flat plane at toe level that `morpho_update` works with; the two are equal whenever the toe is at `z_max`.
- **Beach budget.** The sand seaward of the toe is carried as a state variable, `seaward_vol`, and changes only by real fluxes and exact relabels:

```latex
seaward\_vol_1 = seaward\_vol_0 + Q_{beach} + toe\_relabel - dblock
```

Q\_beach is (−qwnet − qslr + qnour − qbar)·dt + vnourb + eroded volume, as today. `toe_relabel` is the dune volume relabelled by the shift, and `dblock` the sand under the footprint the dune gained or lost (`block_volume`). The shoreline is then solved from `seaward_vol`:

- toe on the flat: `seaward_vol = (x_berm - toe) * (dfoot + dclos) + W_face * dclos + (W_face - bl) / al`;
- toe on the slope: `seaward_vol = W * dclos + (W - bl) / al_eff`, with `al_eff = al * dfoot / h_toe`.

`vbeachtot_phantom` and `phantom_correction` are not used in this mode.

**4.4 Hydrodynamics with a berm**

- Erosion threshold and impact height use the toe level, as in change 1.
- Foreshore slope for runup: the beach face (berm point to shoreline), not the average from toe to shoreline.
- Dissipation width X: beach face plus the flat in front of the toe.
- Dry-beach fetch: waterline on the beach face; the flat is dry while the total water level is below the toe.

**4.5 Pros and cons**

|  | Pros | Cons and risks |
| --- | --- | --- |
| Mass balance | No phantom volume in either direction. On test16d the beach budget closes to about 1e-12 m³/m per step through three storms with dike exposure. | The base under a prograded dune is one straight line from berm point to toe. When the beach slope changes, small volumes are relabelled between dune and beach. The budget stays exact, but `vdune_true` can change slightly without transport. |
| Erosion | The storm profile is the expected one: a flat at dune-foot level and a wider beach face. Today the whole beach pivots around the retreated toe. | While the toe is below `z_max`, the retreat is computed on a flat plane at toe level and re-levelled every step. The error is second order in the step size, but it is an approximation. |
| Calibration | The toe is never above `dfoot`, so the erosion threshold is never higher than the one you calibrated. | Runup slope, dissipation width and fetch change as soon as a flat exists. Results differ from the current model after the first erosion event even with the toe at 3 m (test16d: final shoreline 518.2 m against 517.2 m, eroded volume 425.2 against 426.0 m³/m). |
| Progradation | Self-limiting through the hydrodynamics; no upper clamp or per-step limit to tune. | Not demonstrated yet: in my tests the toe went at most 0.15 m below `z_max`. The floor `dfoot_min` is still needed and its value is a guess. |
| Berm | One extra state variable and one pivot. | Two more parameters to choose (`beach_slope_min`, `beach_slope_max`). Between the limits the berm point is fixed, and the flat left by a storm stays until the dune regrows over it or the slope limit moves it. |
| Beach shape | Identical to today while the toe is on the flat. | On the slope the `al`, `bl` relation is scaled with toe height. Your `bl` is 38.9 m, not zero, so this is an approximation there. |
| Code | `morpho_update` logic untouched; the dike branches keep their geometry. | A toe level that is not a round number exposes float comparisons in `morpho_update`. One lost crest deposits until I fixed it (section 7); others may exist. The new `berm` key appears in `dune_x` and `dune_z`, so scripts that loop over all pivots see one more point. |

**Alternative I considered:** keep the toe at `dfoot` always and only make the budget exact, charging the beach for the bed raised under an advancing toe. It removes the phantom volume with less code, but the dune cannot limit its own progradation, which is the point of your suggestion. I did not pursue it.

## 5. Step-by-step implementation

Steps 0 to 10, in this order, with a test gate after steps 3, 6 and 8 (section 6). Steps 0–8 deliver change 1; the berm option adds to steps 1 to 8. The sandbox patch (`cs_hybrid_slr_toe_prototype.diff`, regenerated on 13 November) contains steps 1–7. It also contains one sandbox-only workaround in `morpho_update`, marked as such, that is not part of this plan (section 7, first row).

### Step 0. Freeze a baseline

1. Commit the current code and tag it.
2. Run and keep the pickles of three reference cases with the unmodified code: Katwijk storm case, Oosteroever 3-year case, and one year of the 30-year case with `slr > 0`.
3. Write `compare_runs.py`: loads two pickles and checks `np.array_equal` on `vdune`, `vbeach`, `vramp`, `eroded_volume`, `overwash_volume`, `qwnet`, `runup`, and dictionary equality of every `dune_x` and `dune_z` entry. It prints the first step that differs.

### Step 1. Add five helpers at the end of `functions.py`

Five functions, tested in the sandbox on the synthetic cases and on test16d. `_reanchor_pivots` holds the geometry, `shift_base_plane` adds the volume measurement, and `beach_base_level` with `find_toe_level` serve the berm option.

```python
def block_volume(x_toe, z_plane, z_bot, x3, z3, m_face):
    """
    Bookkeeping volume of the sand below the base plane, landward of the dune
    toe and seaward of the seaward dike face (taken as a straight line through
    the seaward dike crest), down to z_bot. Signed and exact:
        d/dx_toe = z_plane - z_bot ;  d/dz_plane = x_toe - x4(z_plane)
    Only differences taken with the same z_bot are meaningful.
    """
    z_mid = 0.5 * (z_plane + z_bot)
    return (z_plane - z_bot) * (x_toe - (x3 + (z3 - z_mid) / m_face))


def _reanchor_pivots(dune_x, dune_z, substrate_x, substrate_z, z_new, m_face, tol_z=1e-9, tol_x=1e-9):
    """
    Geometry of the base-plane shift (no volumes): substrate points 4 and 5 go
    to z_new (point 4 slides along the seaward dike face) and the dune pivots
    resting on the seaward plane are re-anchored. Returns new dune_x, dune_z,
    substrate_x, substrate_z and info. Raises ValueError if the plane would
    bury the dune front or the lee crest.
    """
    z_old = float(substrate_z[4])
    z_new = float(z_new)
    info = {'case': 'none', 'lee': 'none'}
    dx, dz = dict(dune_x), dict(dune_z)
    sx, sz = [float(v) for v in substrate_x], [float(v) for v in substrate_z]
    if z_new == z_old:
        return dx, dz, sx, sz, info
    x3, z3 = sx[3], sz[3]
    x4_old = sx[4]
    if z_new >= z3 - 0.05:
        raise ValueError(f"base plane ({z_new:.3f} m) reaches the dike crest ({z3:.3f} m): not supported")
    x4_new = x3 + (z3 - z_new) / m_face
    sx[4], sz[4], sz[5] = x4_new, z_new, z_new
    sx[5] = max(sx[5], x4_new + 1.0)

    def on_plane(k):
        return (dx.get(k) is not None and dz.get(k) is not None
                and abs(dz[k] - z_old) <= tol_z and dx[k] >= x4_old - tol_x)

    # ---------------- seaward pivots ----------------
    walk = [k for k in ('ramp_toe', 'ramp_toe_temp', 'ramp_crest', 'ramp_dike_temp', 'seaward_crest')
            if dx.get(k) is not None and dz.get(k) is not None]   # seaward -> landward
    toe_level = (dx.get('ramp_toe') is not None and abs(dz['ramp_toe'] - z_old) <= tol_z)
    if toe_level and dx['ramp_toe'] < x4_old - tol_x:
        # toe drawn at plane level under the dike face (scarp landward of the dike toe): it keeps its
        # position and only follows the level; partial_volume clips it against the dike
        info['case'] = 'under_face'
        for k in walk:
            if abs(dz[k] - z_old) <= tol_z and dx[k] < x4_old + tol_x:
                dz[k] = z_new
    elif on_plane('ramp_toe'):
        xt = dx['ramp_toe']
        if abs(xt - x4_old) <= tol_x:
            # toe at the dike toe: follows point 4
            info['case'] = 'dike_toe'
            for k in walk:
                coincident = abs(dx[k] - x4_old) <= tol_x and abs(dz[k] - z_old) <= tol_z
                buried = (z_new > z_old and dz[k] < z_new and x3 - tol_x <= dx[k] <= x4_old + tol_x)
                if coincident or buried:
                    dx[k], dz[k] = x4_new, z_new
        elif z_new > z_old:
            info['case'] = 'rise'
            a = next((i for i, k in enumerate(walk) if dz[k] > z_new + tol_z), None)
            if a is None or a == 0:
                raise ValueError("base plane rose above the dune front: dune buried (not handled)")
            ka, kp = walk[a], walk[a - 1]
            if abs(dx[ka] - dx[kp]) <= tol_x:
                xc = dx[kp]
            else:
                xc = dx[kp] + (z_new - dz[kp]) * (dx[ka] - dx[kp]) / (dz[ka] - dz[kp])
            xc = max(xc, x4_new)
            for k in walk[:a]:
                dx[k], dz[k] = xc, z_new
        else:
            info['case'] = 'drop'
            group = [k for k in walk if abs(dx[k] - xt) <= tol_x and abs(dz[k] - z_old) <= tol_z]
            q = next((k for k in walk if k not in group), None)
            xc = xt
            if q is not None and (xt - dx[q]) > tol_x and (dz[q] - z_old) > tol_z:
                xc = xt + (z_old - z_new) * (xt - dx[q]) / (dz[q] - z_old)
            for k in walk:
                if k in group:
                    dx[k], dz[k] = xc, z_new
                elif on_plane(k):
                    dz[k] = z_new
    # ---------------- lee foot resting on the seaward plane ----------------
    if on_plane('lee_foot'):
        xl, xlc, zlc = dx['lee_foot'], dx['lee_crest'], dz['lee_crest']
        follow_ref = dx['back_ref'] >= x4_old - tol_x     # back_ref also sits on the seaward plane
        if (xlc - xl) > tol_x and (zlc - z_old) > tol_z:
            sl = (zlc - z_old) / (xlc - xl)
            xn, zn = xl + (z_new - z_old) / sl, z_new
            info['lee'] = 'slide'
            if xn < x4_new:      # lee line meets the dike face before the plane
                xn = (z3 + m_face * x3 - zlc + sl * xlc) / (sl + m_face)
                zn = z3 - m_face * (xn - x3)
                info['lee'] = 'on_face'
            if zn >= zlc - tol_z:
                raise ValueError("base plane rose above the lee crest: dune buried (not handled)")
        else:
            xn, zn = max(xl, x4_new), z_new
            info['lee'] = 'vertical'
        dx['lee_foot'], dz['lee_foot'] = xn, zn
        if follow_ref:
            dx['back_ref'], dz['back_ref'] = xn, zn
    return dx, dz, sx, sz, info


def shift_base_plane(dune_x, dune_z, substrate_x, substrate_z, z_new, m_face,
                     tol_z=1e-9, tol_x=1e-9):
    """
    Move the seaward base plane (substrate points 4 and 5 = dune toe level) to
    z_new and re-anchor the dune pivots that rest on it.
    Returns dune_x, dune_z, substrate_x, substrate_z (new objects), tau, info
    tau > 0: dune volume relabelled to below-plane (beach side) sand.
    """
    if float(z_new) == float(substrate_z[4]):
        return dune_x, dune_z, substrate_x, substrate_z, 0.0, {'case': 'none', 'lee': 'none'}
    v_before = partial_volume(dune_x, dune_z, substrate_x, substrate_z,
                              dune_x['lee_foot'], dune_x['ramp_toe'])
    dx, dz, sx, sz, info = _reanchor_pivots(dune_x, dune_z, substrate_x, substrate_z, z_new, m_face, tol_z, tol_x)
    v_after = partial_volume(dx, dz, sx, sz, dx['lee_foot'], dx['ramp_toe'])
    return dx, dz, sx, sz, v_before - v_after, info


def beach_base_level(x, x_berm, z_max, x_toe0, z_toe0, x_shore, z_msl, z_floor, slope_min=0.0):
    """
    Elevation of the beach base profile at x: flat at z_max landward of the berm
    point, then the line berm -> previous toe (if that toe was seaward of the
    berm and below z_max) -> shoreline. Not lower than z_floor.
    """
    if x <= x_berm:
        return z_max
    xs, zs = [x_berm], [z_max]
    if x_toe0 > x_berm + 1e-9 and z_toe0 < z_max - 1e-12 and x_toe0 < x_shore - 1e-9:
        xs.append(x_toe0); zs.append(z_toe0)
    xs.append(max(x_shore, xs[-1] + 1e-6)); zs.append(z_msl)
    zb = max(float(np.interp(x, xs, zs)), z_floor)
    if slope_min > 0:   # the beach in front of a low toe is not flatter than slope_min: it builds up instead
        zb = max(zb, min(z_max, z_msl + slope_min * (x_shore - x)))
    return zb


def find_toe_level(dune_x, dune_z, substrate_x, substrate_z, m_face, base_level, z_floor, z_max):
    """
    Level of the base plane for which the dune toe lies on the beach base
    profile (base_level(x)). z_max if the toe is on the flat part; between
    z_floor and z_max if the dune front stands on the beach slope.
    """
    def g(z):
        try:
            dx = _reanchor_pivots(dune_x, dune_z, substrate_x, substrate_z, z, m_face)[0]
        except ValueError:
            return 1.0          # plane above the dune front: the toe must be lower
        return z - base_level(dx['ramp_toe'])
    if g(z_max) <= 1e-12:
        return z_max
    if g(z_floor) >= 0.0:
        return z_floor
    return float(brentq(g, z_floor, z_max, xtol=1e-12, rtol=1e-14))
```

Two existing places in `morpho_update` also change. Neither alters the output with the flags off (checked bit for bit on test16d and the synthetic cases).

**Line 2997.** The condition lacks the x test that its twin at line 3337 has. With a rising plane, `back_ref` ends up seaward of the lee foot and `partial_volume` sees a notch.

```python
# functions.py line 2997, old
if dune_z['lee_foot'] > dune_z['back_ref']:
# new
if dune_x['back_ref'] > dune_x['lee_foot'] or dune_z['lee_foot'] > dune_z['back_ref']:
```

**Before `if qs_crest > 0:` (line 2470).** The crest branch tests compare elevations with `<` and `==`. With a toe level that is not a round number, `(z - dfoot) + dfoot` can differ from `z` by one bit, the detached-crest branch is taken by mistake, and the crest deposit of that step is lost. Insert:

```python
for _k in ('lee_crest', 'crest_temp', 'ramp_dike_temp', 'ramp_crest'):
    if dune_z.get(_k) is not None and dune_z[_k] != dune_z['seaward_crest'] and abs(dune_z[_k] - dune_z['seaward_crest']) <= 1e-9:
        dune_z[_k] = dune_z['seaward_crest']
```

### Step 2. State and flags in `main.py`

Insert before `# Initialize results dictionary` (line 1281):

```python
slr_toe = _config_flag(getattr(p, 'slr_toe', False))
dynamic_toe = _config_flag(getattr(p, 'dynamic_toe', False))
swl_includes_slr = _config_flag(getattr(p, 'swl_includes_slr', False))
msl = 0.0
z_plane = float(substrate_z[4])
h_ref = float(p.dfoot)
m_face = (substrate_z[3] - substrate_z[4]) / (substrate_x[4] - substrate_x[3])
face_x3, face_z3 = float(substrate_x[3]), float(substrate_z[3])
substrate_x0, substrate_z0 = list(substrate_x), list(substrate_z)
if slr_toe or dynamic_toe:
    if abs(substrate_z[4] - p.dfoot) > 1e-12 or abs(substrate_z[5] - p.dfoot) > 1e-12:
        raise ValueError("slr_toe/dynamic_toe need base_elevation_front == dfoot")
    if abs(m_face - p.betas_hs) > 1e-6 * max(1.0, abs(p.betas_hs)):
        log(f"Warning: betas_hs ({p.betas_hs:g}) differs from the dike face slope of the geometry ({m_face:g})")
    log(f"Sea-level rise over the run: {p.slr * p.dt * time_length_dt:.4f} m "
        f"({p.slr * 365.25 * 86400 * 1000:.2f} mm/yr)")
h_min = float(getattr(p, 'dfoot_min', 0.3 * h_ref))      # lowest toe height above mean sea level
# beach-face slope limits (toe height above sea level / face width): steeper than slope_max -> the berm point
# retreats; flatter than slope_min -> the berm point advances (0 = no lower limit)
slope_max = float(getattr(p, 'beach_slope_max', h_ref / max(p.bl + 1.0, 30.0)))
slope_min = float(getattr(p, 'beach_slope_min', 0.0))
w_face_min = max(h_ref / slope_max, p.bl + 1.0)
w_face_max = h_ref / slope_min if slope_min > 0 else float('inf')
x_berm = float(dune_x['ramp_toe'])                         # berm point: seaward end of the flat at dfoot
block_cum = 0.0
if dynamic_toe:
    # berm pivot between ramp_toe and shoreline (visible berm crest; equals the toe when the dune covers it)
    _xs, _zs = dune_x.pop('shoreline'), dune_z.pop('shoreline')
    dune_x['berm'], dune_z['berm'] = x_berm, dune_z['ramp_toe']
    dune_x['shoreline'], dune_z['shoreline'] = _xs, _zs
# sand seaward of the toe down to the depth of closure (berm block + beach face), carried as state
seaward_vol = (dune_x['shoreline'] - dune_x['ramp_toe']) * p.dclos + ((dune_x['shoreline'] - dune_x['ramp_toe']) - p.bl) / p.al
vbeach_x = None
```

Three points about this block:

- The log line with the total rise is the unit check for `slr` (it must be in m/s).
- The `berm` key must sit between `ramp_toe` and `shoreline` in the dictionary. `morpho_update` interpolates over all pivot values in insertion order (lines 1529, 1577), so the order has to stay monotonic in x. For the same reason the stored pivot is the visible berm crest, `max(x_berm, toe)`; the true berm position is kept in `x_berm` and in `results['x_berm']`.
- Add `'berm'` to `_profile_order` in `_save_profile_snapshot` (line 190).

### Step 3. Hydrodynamics and aeolian block (`main.py` lines 1568–1737)

At line 1568, replace the SWL line by:

```python
swl = swl_time_series[i] + (0.0 if swl_includes_slr else msl)
h_beach = z_plane - msl
dfoot_abs = z_plane
al_eff = p.al * h_ref / h_beach if dynamic_toe else p.al
```

Then edit each `p.dfoot` as follows. Every line gives the same number as today when `msl = 0` and the plane has not moved.

| Line | Now | Becomes |
| --- | --- | --- |
| 1613, 1615 | `dune_z['ramp_toe'] / ...` | `(dune_z['ramp_toe'] - dune_z['shoreline']) / ...` |
| 1618 | `2 * vbeach / p.dfoot * (1 - swl / p.dfoot)` | `2 * vbeach / h_beach * (1 - (swl - msl) / h_beach)` |
| 1623, 1625 | `p.dfoot` in `rprim` and its test | `h_beach` |
| 1644–1658 | `p.dfoot` (11 times: threshold, impact height, crest height) | `dfoot_abs` |
| 1722 | `(1 - (r + swl) / p.dfoot)` | `(1 - (r + swl - msl) / h_beach)` |
| 1726 | `r + swl >= p.dfoot` | `r + swl >= dfoot_abs` |
| 1735, 1737 | `... - p.dfoot` | `... - dfoot_abs` |

For the berm option, four more edits in the same block. All are inside `if dynamic_toe`, so they do not affect gate A.

```python
# after the al_eff line (line 1568 block)
toe_on_flat = (h_beach >= h_ref - 1e-9)

# foreshore slope (lines 1612-1615): put this branch first, keep the two existing ones as elif / else
if dynamic_toe:
    slope_foreshore = (dune_z['berm'] - dune_z['shoreline']) / max(dune_x['shoreline'] - dune_x['berm'], 30.0)

# dissipation width (line 1618)
X = 2 * (vbeach_x if (dynamic_toe and vbeach_x is not None) else vbeach) / h_beach * (1 - (swl - msl) / h_beach)

# dry-beach limit (line 1722)
if dynamic_toe and toe_on_flat:
    yr = dune_x['berm'] + (1 - (r + swl - msl) / h_beach) * (dune_x['shoreline'] - dune_x['berm'])
else:
    yr = dune_x['ramp_toe'] + (1 - (r + swl - msl) / h_beach) * (dune_x['shoreline'] - dune_x['ramp_toe'])
```

**Gate A:** run T1 (section 6) now. It must be bit-identical.

### Step 4. Save the start-of-step beach state (`main.py` lines 1766–1768)

Replace `p.al` by `al_eff` and `p.dfoot` by `h_beach` in the three lines, then add:

```python
x_toe_0 = dune_x['ramp_toe']; w_0 = dune_x['shoreline'] - dune_x['ramp_toe']
wedge_0 = w_0 * p.dclos + vbeach_old
z_bot = msl - p.dclos
x_shore_0 = dune_x['shoreline']
if dynamic_toe and coupling_mode:
    # shoreline imposed from outside: rebuild the seaward volume from the geometry
    if toe_on_flat:
        seaward_vol = (x_berm - x_toe_0) * (h_ref + p.dclos) + (x_shore_0 - x_berm) * p.dclos + ((x_shore_0 - x_berm) - p.bl) / p.al
    else:
        seaward_vol = w_0 * p.dclos + (w_0 - p.bl) / al_eff
```

### Step 5. The shift, right after `morpho_update` (insert before line 1781)

```python
toe_relabel = 0.0
shift_info = {'case': 'none', 'lee': 'none'}
z_plane_0, msl_0 = z_plane, msl
if slr_toe or dynamic_toe:
    dz_slr = p.slr * p.dt if slr_toe else 0.0
    msl = msl + dz_slr
    z_max = msl + h_ref                       # highest toe level: rises with sea level only
    if dynamic_toe:
        # the toe must lie on the beach base profile: flat at z_max up to the berm point, beach slope seaward of it
        base = lambda x: beach_base_level(x, x_berm, z_max, x_toe_0, z_plane_0 + dz_slr, x_shore_0, msl, msl + h_min, slope_min)
        z_target = find_toe_level(dune_x, dune_z, substrate_x, substrate_z, m_face, base, msl + h_min, z_max)
    else:
        z_target = z_plane + dz_slr
    if z_target != z_plane:
        dune_x, dune_z, substrate_x, substrate_z, toe_relabel, shift_info = shift_base_plane(
            dune_x, dune_z, substrate_x, substrate_z, z_target, m_face)
        z_plane = float(substrate_z[4])
    dune_z['shoreline'] = msl
h_beach_1 = z_plane - msl
al_eff_1 = p.al * h_ref / h_beach_1 if dynamic_toe else p.al
toe_on_flat_1 = (h_beach_1 >= h_ref - 1e-9)
```

`substrate_x` and `substrate_z` are rebound to new lists here, so `results['substrate_x']` keeps the initial geometry. `vramp` and `vdune` (lines 1824–1825) are computed after this block and need no change.

### Step 6. Beach update (`main.py` lines 1833–1847)

6a (change 1) keeps the current equations with `h_beach_1` and `al_eff_1`. 6b (berm option) adds the exact branch. Replace the block by:

```python
q_beach = (-qwnet - qslr + qnour - qbar) * p.dt + vnourb[i] + deltavramp_ero
dblock = (block_volume(dune_x['ramp_toe'], z_plane, z_bot, face_x3, face_z3, m_face)
          - block_volume(x_toe_0, z_plane_0, z_bot, face_x3, face_z3, m_face))
if dynamic_toe:
    seaward_vol = seaward_vol + q_beach + toe_relabel - dblock
    xt = dune_x['ramp_toe']
    if toe_on_flat_1:
        x_berm = max(x_berm, xt)
        if not coupling_mode:
            w_face = (p.al * (seaward_vol - (x_berm - xt) * (h_ref + p.dclos)) + p.bl) / (1 + p.al * p.dclos)
            if (w_face < w_face_min and x_berm > xt) or w_face > w_face_max:
                # face steeper than slope_max: the berm point erodes back (not past the toe);
                # face flatter than slope_min: the berm point builds out and the flat widens
                w_lim = w_face_min if w_face < w_face_min else w_face_max
                x_berm = max(xt, xt + (seaward_vol - w_lim * p.dclos - (w_lim - p.bl) / p.al) / (h_ref + p.dclos))
                w_face = (p.al * (seaward_vol - (x_berm - xt) * (h_ref + p.dclos)) + p.bl) / (1 + p.al * p.dclos)
            dune_x['shoreline'] = x_berm + w_face
        w_face = dune_x['shoreline'] - x_berm
        vbeach = (x_berm - xt) * h_ref + (w_face - p.bl) / p.al            # sand above mean sea level seaward of the toe
        vbeach_x = 0.5 * (x_berm - xt) * h_ref + (w_face - p.bl) / p.al    # equivalent volume for the dissipation width X
        dune_x['berm'], dune_z['berm'] = x_berm, z_plane
    else:
        if not coupling_mode:
            w_1 = (al_eff_1 * seaward_vol + p.bl) / (1 + al_eff_1 * p.dclos)
            dune_x['shoreline'] = xt + w_1
        vbeach = ((dune_x['shoreline'] - xt) - p.bl) / al_eff_1
        vbeach_x = vbeach
        dune_x['berm'], dune_z['berm'] = xt, z_plane                        # berm point is under the dune
else:
    vbeachtot_phantom = dune_x['ramp_toe'] * (p.dclos + h_beach_1) + (dune_x['shoreline'] - dune_x['ramp_toe']) * (p.dclos + db_old)
    vbeachtotnew = vbeachtot_phantom + (-qwnet - qslr + qnour - qbar) * p.dt + vnourb[i] + deltavramp_ero
    if not coupling_mode:
        dune_x['shoreline'] = dune_x['ramp_toe'] + (al_eff_1 * (vbeachtotnew - dune_x['ramp_toe'] * (h_beach_1 + p.dclos)) + p.bl) / (1 + al_eff_1 * p.dclos)

if not dynamic_toe:
    vbeach = ((dune_x['shoreline'] - dune_x['ramp_toe']) - p.bl) / al_eff_1
block_cum = block_cum + dblock
```

In the `else` branch keep the sum written out in the original order. Writing it as `vbeachtot_phantom + q_beach` changes the rounding and breaks the bit-identical regression (it did in the sandbox).

**Gate B:** run T1 to T4.

### Step 7. Per-step conservation and outputs

`mass_conservation.py`:

- `calculate_compartment_volumes(..., dfoot, block=None)`: when `block` is given, `beach = block + (shoreline - ramp_toe) * dclos + vbeach`; otherwise the current formula.
- `calculate_timestep_conservation`: read `qslr`, `toe_relabel` and `exact_beach` from `fluxes` and use

```python
qslr_vol = fluxes.get('qslr', 0.0) * dt
relabel = fluxes.get('toe_relabel', 0.0)
exact = fluxes.get('exact_beach', False)
expected_dune_change  = total_aeolian - eroded + nour_dune - relabel
expected_beach_change = -qwnet_vol - qbar_vol + eroded + nour_beach - qslr_vol + (relabel if exact else 0.0)
expected_bar_change   = qbar_vol + nour_bar
expected_total_change = nour_dune + nour_beach + nour_bar - qslr_vol - (0.0 if exact else relabel)
```

`main.py`:

- Both calls to `calculate_compartment_volumes` (lines 1400 and 1857) take a `block` argument. Initial call: `block=(block_cum if dynamic_toe else None)`. In the loop, pass `h_beach_1` instead of `p.dfoot` and `block=(block_cum + seaward_vol - ((dune_x['shoreline'] - dune_x['ramp_toe']) * p.dclos + vbeach) if dynamic_toe else None)`, so that the beach compartment equals `block_cum + seaward_vol`. Forgetting the initial call gives a single beach error of about 1500 m³/m at step 0.
- `fluxes` (line 1862): add `'qslr': qslr, 'toe_relabel': toe_relabel, 'exact_beach': dynamic_toe`.
- New entries in `results`, initialised with the start value where one exists and appended every step: `dfoot` (z\_plane), `msl`, `substrate_x4`, `x_berm`, `block` (block\_cum), `qslr`, `toe_relabel`, `dblock`, `shift_case` (`info['case'] + '/' + info['lee']`), `vdune_ref` (`partial_volume` with `substrate_x0`, `substrate_z0`), and the flags `dynamic_toe`, `slr_toe`.
- `vdune_true` (berm option): `partial_volume` over the beach base profile. Copy the current substrate, set point 4 to `z_max` on the dike face, and when the toe is seaward of the berm point append `(x_berm, z_max)` and `(toe, z_plane)`. With the option off, store `vdune`.

### Step 8. Post-run check (`check_mass_conservation`, lines 269–396)

1. Build the substrate of a given step from the results instead of using one substrate for both ends: copy `results['substrate_x']` and `['substrate_z']`, then set `x[4] = results['substrate_x4'][k]` and `z[4] = z[5] = results['dfoot'][k]`. Use k = 0 (or 1 with `skip_first_step`) for the initial dune volume and k = −1 for the final one.
2. 8a (change 1): `expected_total_change = phantom_correction + total_nourishment - sum(qslr) * dt - sum(toe_relabel)`. The phantom loop stays as it is, because the toe height above mean sea level equals the config `dfoot`.
3. 8b (berm option): when `results['dynamic_toe']` is true, take both beach volumes from the per-step conservation record (`results['conservation']['beach_volume']`), skip the phantom loop, and use `expected_total_change = total_nourishment - sum(qslr) * dt`.
4. Add `sum(qslr) * dt` to `total_transport` and print it in the report.

This step is not in the sandbox patch. The balances in section 6 were checked with the per-step record and a separate script.

**Gate C:** run T5 and T6.

### Step 9. Plotting and external scripts

- `_save_profile_snapshot` draws the sea at 0.0 (lines 638, 664, 688). Pass `msl` through `bc_info` and use it there. The substrate passed in the loop is already the current one. With `'berm'` in `_profile_order` the flat and the beach face are drawn without further changes.
- `plot_oosteroever_decadal.py` and the Katwijk scripts: the dune-toe line and the storm threshold should read `results['dfoot']` per step; dune volume comparisons between runs should use `vdune_true` or `vdune_ref`.
- `run_example_test16d.py` works unchanged: `extract_dune_section` keeps only pivots between lee foot and toe, and the profile plot simply gains the berm point.

### Step 10. Yaml and README

| Key | Default | Meaning |
| --- | --- | --- |
| `slr_toe` | false | toe and seaward flat rise with `slr` |
| `swl_includes_slr` | false | set true only if the SWL forcing already contains the sea-level trend |
| `dynamic_toe` | false | berm point and maximum toe level; exact beach budget |
| `dfoot_min` | 0.3 × `dfoot` | lowest toe height above mean sea level (numerical floor) |
| `beach_slope_max` | `dfoot` / max(`bl` + 1, 30 m) | steepest beach face; beyond it the berm point retreats |
| `beach_slope_min` | 0 (no limit) | flattest beach face; beyond it the berm point advances |

`dfoot` itself becomes the maximum toe height above mean sea level. The keys `dfoot_max` and `dz_toe_max_step` of the first version of this plan are gone.

## 6. Verification tests

Nine tests, each with a pass criterion that needs no judgement. The sandbox column covers two synthetic geometries and your test16d case (79,056 steps, three storms, dike exposure and detachment). Test16d ran on the uploaded code plus one workaround (section 7, first row), so every row has to be repeated on your current code.

| Test | Setup | Pass criterion | Sandbox |
| --- | --- | --- | --- |
| T1 Regression | Flags off | `compare_runs.py` reports bit-identical output | Passed, synthetic and test16d |
| T2 Neutral switch | `slr_toe: true`, `slr: 0` | Bit-identical to baseline | Passed, synthetic and test16d |
| T3 Geometry invariants | `slr_toe: true`, `slr` > 0, assert every step | `dune_z['ramp_toe'] == substrate_z[4] == substrate_z[5]` when the toe is on the plane; `substrate_x[4] == x3 + (z3 - z_plane) / m_face`; `z_plane - msl == dfoot` within 1e-9 | Passed on test16d with 0.5 m of rise, zero violations |
| T4 SLR balance | Same run | Final total change − phantom correction + Σqslr·dt + Σtoe\_relabel differs from the baseline value only by the difference in cumulative `dune_error` | Passed on test16d: both differ by −0.021 m³/m (dune error −0.059 against −0.038) |
| T5 Exact budget | `dynamic_toe: true`, with and without SLR | Per-step `beach_error` below 1e-9 m³/m; total error equals dune error | Passed on test16d: largest step error 1e-12, total −0.039 (no SLR) and −0.060 (0.5 m rise), equal to the dune error |
| T6 Dike contact | test16d with each option | No exception; `shift_case` contains `dike_toe` and `under_face`; T3 to T5 hold through exposure, detachment and recovery | Passed: 357 steps at the dike toe, 9 under the dike face, 138 detached steps |
| T7 Berm behaviour on erosion | test16d, `dynamic_toe: true` | Toe stays at `dfoot` and retreats on the flat; berm point fixed; `ramp_toe <= berm <= shoreline` every step | Passed: flat 54 to 60 m wide after each storm, 13 m at the end; shoreline 1 to 3 m seaward of the baseline |
| T8 Toe on the beach slope | A case where the dune progrades past the berm point, then a storm | Toe level falls on advance, returns to `dfoot` on retreat, never below `dfoot_min`; T5 holds | Partly: synthetic case only, toe 0.15 m below `dfoot`, budget exact. In test16d the toe never passes the berm point |
| T9 Sanity of magnitudes | Oosteroever 30-year case, 4 mm/yr | `dfoot` ends at start + 0.12 m; toe retreat from the shift alone about 0.12 / betas; shoreline retreat close to the Bruun-only run | Not tested |

Two extra checks that cost minutes and catch most mistakes:

- **Operator unit test.** Take 20 saved states from an existing pickle (`results['dune_x'][k]`, including storm steps), apply `shift_base_plane` with dz = ±0.001 and ±0.05 m, and check: `tau` has the sign of dz, |tau| ≤ 1.05 × dz × (toe − dike toe distance) + 0.01, and one following `morpho_update` call with a small `qs_ramp` runs without exception.
- **`shift_case` histogram.** Print `collections.Counter(results['shift_case'])` after every test run. It shows at once which states were exercised and which were not.

## 7. Troubleshooting

Start every debugging session with two outputs: the first step that differs (`compare_runs.py`) and the `shift_case` of that step.

**Found on test16d and in the berm tests (13 November)**

| Symptom | Cause | Fix |
| --- | --- | --- |
| The uploaded `functions.py` crashes on `config_test16d.yaml` at step 6243 with flags off: `brentq`, f(a) and f(b) must have different signs, line 2416 | After the first storm the ramp is rebuilt on the dike face under a sand scarp at x ≈ 295 m. When it reaches the scarp foot it keeps climbing the face at about 1.2 mm per step, because the 1 cm tolerance at line 2034 keeps it in the dike branch, and 9 steps later it falls out of that tolerance. Your GIF shows the ramp continuing against the scarp, and a different lee side after the third storm, so the uploaded file is not the version that made the GIF | Use your current `functions.py`. The sandbox uses a workaround (flag `ramp_at_scarp` before line 2034, marked in the diff) only to get a baseline; do not merge it |
| `brentq` crash in the ramp routine about 1700 steps after the first storm, only with `slr_toe` | A scarp landward of the dike toe leaves `ramp_toe` at plane level under the dike face. The first operator did not move it, so the toe level froze while the plane kept rising | Case `under_face` in `_reanchor_pivots` (step 1) |
| Crest deposit of single steps lost with a moving toe level (−0.68 m³/m in 14,000 steps of the synthetic berm test) | `ramp_dike_temp` one bit above `seaward_crest` after `(z - dfoot) + dfoot`; the detached-crest branch at line 2476 is taken | The level snap before `if qs_crest > 0:` (step 1). Residual afterwards +0.006 m³/m |
| `np.interp` in the erosion branch returns nonsense after adding the `berm` pivot | `berm` appended after `shoreline`, or stored landward of the toe, so the pivot x values are not monotonic | Insert the key between `ramp_toe` and `shoreline` and store `max(x_berm, toe)` (step 2) |

**General**

The first four rows are problems I hit in the first sandbox round.

| Symptom | Likely cause | Fix |
| --- | --- | --- |
| T1 fails with differences of 1e-12 that start at the first storm | Rounding only: the beach sum was regrouped, or `al_eff` was computed as `al * h_ref / h_beach` with the flags off | Keep the original expression order in the non-dynamic branch; `al_eff = p.al` unless `dynamic_toe` |
| One beach error of about `ramp_toe * (dclos + dfoot)` at step 0 with `dynamic_toe` | Initial `calculate_compartment_volumes` call still uses the old block formula | Pass `block=block_cum` in the call at line 1400 too |
| Dune errors of ±1e-3 m³/m that cancel in pairs, dune in front of the dike, rising plane | `back_ref` seaward of `lee_foot` gives a notch in `partial_volume` | The one-line fix at functions.py line 2997 (step 1) |
| Crash or dune error in an erosion step that also occurs with the flags off | Existing behaviour, not caused by the change. In the sandbox the original code reported eroded volume without changing the profile once a scarp reached the lee crest, and crashed in `brentq` when erosion reached the dike | Reproduce with flags off first. Fix separately, not in this branch |
| `brentq`: f(a) and f(b) must have different signs, only with flags on | A pivot is a hair off the plane after the shift, so a bracket built from `dfoot` has zero or negative width | Check that the plane-resting pivots received the same float as `substrate_z[4]`. Never recompute `z_plane` by another route than `float(substrate_z[4])` |
| A branch that used to trigger no longer does (e.g. dike-exposed recovery) | Float equality on elevations (functions.py lines 1285, 1795, 2513, 2651, 3321) or on x (lines 1592, 1606, 1805, 1836, 2064) fails once `dfoot` is not a round number, because `(z3 - z4) + z4` can differ from `z3` by one bit | Contingency: define `z_dike = substrate_z[3]` at the top of `morpho_update` and replace the 12 `s_hs + dfoot` by it; turn line 1795 into `abs(dune_z['seaward_crest'] - z_dike) < 1e-9`. Run T1 after it |
| `shift_base_plane` raises "reaches the dike crest" | Plane within 5 cm of the dike crest: low dike or wrong `slr` unit | Check the logged total rise. A buried dike is outside this plan |
| `shift_base_plane` raises "dune buried" | Plane rose above the dune front or lee crest (very small dune, large rise) | Start from a higher incipient dune, or decide how a buried dune should reset before handling it |
| Total rise in the log is absurd | `slr` given per year or per step | `slr` is in m/s (4 mm/yr = 1.27e-10) |
| Erosion almost disappears with `slr_toe` | `swl_includes_slr` was set to true while the SWL forcing has no trend: the toe rises but the water does not | Leave `swl_includes_slr` at its default (false) |
| Erosion roughly doubles with `slr_toe` | Trend counted twice: it is in the forcing and the model adds it as well (the default) | Set the flag to true |
| Toe level goes down to `dfoot_min` (`dynamic_toe`) | Waves never reach the advancing toe (calm forcing, low runup), so nothing stops the progradation before the floor | Check runup and SWL in the forcing; set `dfoot_min` to the level where you expect the limit. At the floor the beach is charged for the toe advance, so the budget stays exact |
| Shoreline drifts compared with the baseline under `dynamic_toe` | Expected: runup slope, dissipation width and fetch use the beach face and the flat (section 4.4). Larger drifts when the toe sits at `dfoot_min` | Expected within a few metres; if larger, check `al * dfoot` against 2 and `bl` |
| Post-run check fails although per-step totals are fine | Check still uses a single substrate, or `qslr` and `toe_relabel` are missing in the expected change | Step 8 |
| Dune-toe line or storm shading wrong in the paper figures | Plot script uses the config `dfoot` | Read `results['dfoot']` (step 9) |
| `vdune` jumps without transport | The plane moved; the volume is measured above the plane | Use `vdune_ref` for comparisons |

Rules that prevent most of the above:

- One source of truth for the toe level: `z_plane`, set only from `substrate_z[4]` after a shift. `p.dfoot` is used only as `h_ref` and in the non-dynamic beach relation.
- Do not edit `morpho_update` in this branch, apart from line 2997, the level snap before the crest block, and the contingency above.
- Commit after each gate, so `git bisect` stays useful.

## 8. Assumptions and decisions to confirm before coding

Two points are settled by your answers; six remain. Item 1 blocks further testing.

**Settled**

- Sea-level rise is not in the SWL series; the model adds it (`swl_includes_slr` defaults to false).
- The lee side does not rise. Back base, lee foot behind the dike and `dfoot_back` stay fixed; only a lee foot resting on the seaward flat moves with the plane.

**Open**

- [ ] **1. Send the `functions.py` you actually run.**

  The uploaded one crashes on `config_test16d.yaml` at step 6243 with all flags off (section 7, first row). All test16d results here rest on a sandbox workaround.
- [ ] **2. The dune slab below the rising toe stays under the dune (not added to the beach face).**

  This conserves real sediment. On test16d with 0.5 m of rise the slab is 19.5 m³/m in total, against 500 m³/m taken by the Bruun term. Please check which one the original CS-model code does; switching is one line.
- [ ] **3. `base_elevation_front == dfoot` in every config you will run.**

  True in test16d. The plan requires it (step 2 raises an error otherwise).
- [ ] **4. `slr` is constant and in m/s.**

  A time-varying rate would need `slr` read as a forcing series. Not included.
- [ ] **5. Values for `beach_slope_min` and `beach_slope_max`.**

  Implemented as you suggested: the berm point advances at the minimum slope and retreats at the maximum. The defaults keep the berm point fixed on accretion (no minimum) and let it retreat only on a very narrow face. I tested 0.018 and 0.025 on test16d, values I picked around its initial slope of 0.021; they need to come from your sites.
- [ ] **6. Floor `dfoot_min` and the beach shape on the slope.**

  The default 0.3 × `dfoot` is a placeholder. On the slope the `al`, `bl` relation is scaled with toe height, which is approximate for your `bl` = 38.9 m.

Two limits of what I verified: no test has yet taken the toe far down the beach slope, so the hydrodynamic limit on progradation is untested; and the post-run `check_mass_conservation` (step 8) is not updated in the sandbox.
