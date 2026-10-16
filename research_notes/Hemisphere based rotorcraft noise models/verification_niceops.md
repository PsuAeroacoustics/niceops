# Adversarial verification of `niceops_inventory.md` against the code

Compiled 2026-10-03. Each claim was checked by trying to refute it from the source.

**Code base.** All `file:line` citations are to the worktree
`/Users/evg5332/Code/D/NICEOPS/.claude/worktrees/hemisphere-review` at HEAD `fb9fad5`. Paths are relative to it.

The inventory was written against the main checkout at `e472b6f`, which is 10 commits *ahead* of this HEAD. Between the two, only `niceops.cpp`, `trajectory.*`, `ground_effect.*`, `refraction.*` and `stratified.*` differ. `sphere_database.cpp`, `interpolator.cpp`, `util/sel.hpp`, `weighting.cpp` and `propagation.hpp` are identical, so every claim below except the `niceops.cpp` guard details holds for both commits. The inventory's own line numbers are for `e472b6f` and are off by tens of lines here.

## Verdict summary

| # | Claim | Verdict |
| --- | --- | --- |
| 1 | Spectral is the default; `--no-spectral` selects broadband | **CONFIRMED**, with one caveat: the default *demo* run is still broadband |
| 2 | phi uses \|up\|, so the upper hemisphere is never sampled | **CONFIRMED**. It matters for receivers above the source, for TPP tilt in manoeuvres, and for upward-launched refracted rays |
| 3 | SEL window runs first-to-last within 10 dB, not a walk-out from the peak | **CONFIRMED**; CLAUDE.md and README are stale |
| 4 | `-D` is not normalised by 86 400 s | **CONFIRMED** |
| 5 | Queries are clamped, the hull is extended 10%, nearest condition outside the hull | **CONFIRMED**, with clarifications |
| 6 | +20 log10(LF); clamp at 5 g on some databases, unbounded on fixed_load_factor | **PARTLY**: the correction is 20 log10(CT/CT_ref), and the clamp is the database's CT maximum |
| 7 | No Doppler shift, no convective amplification, reception weight (1 + dr/dt / c) | **CONFIRMED** |
| 8 | Condition blend in dB; within-sphere bilinear in dB | **CONFIRMED**: both are in dB, and energy appears only in the band sum |
| 9 | Refraction / ground effect on the broadband path | **CONFIRMED refused**: every route is an input error |
| 10 | No LAmax, EPNL/PNL or multiple-source support | **CONFIRMED** |
| A | README: "Databases are not tracked" | **REFUTED** by `git ls-files` |
| B | CLAUDE.md: "No database in `databases/` carries spectra" | **REFUTED** (new; not in the inventory) |
| C | CLAUDE.md/README: the windowed path and `--db_cutoff` are unavailable with `--spectral` | **REFUTED**: both run on spectral now (new) |
| D | Asserts are compiled out in Release despite the docs | **CONFIRMED** |
| E | The projection is transverse Mercator (docs) or polyconic (code review) | **Adjudicated: polyconic.** The docs are wrong |

---

## 1. Spectral propagation is the default — CONFIRMED (with a caveat)

**Evidence**
- `src/niceops.cpp:546-552`: `RunOptions` carries `bool spectral = true;`, commented "On by default".
- `src/niceops.cpp:2082`: CLI default `bool spectral = true`.
- `src/niceops.cpp:2186-2190`: the flag is `add_flag("--spectral,!--no-spectral", ...)`, with help text "(default on; falls back to broadband when the database carries no spectra). --no-spectral propagates one broadband level...".
- Fallback, `src/niceops.cpp:983-1007`:
  - if the database has no spectra, the run switches to broadband with a log line;
  - but it is an `input_error` if `--spectral` was named explicitly (`spectral_explicit`, `:2517`), or if `--ground_effect`, `--atmosphere` or `--propagation_table` was given.
- Doc side: README `:168` describes `--spectral` as an opt-in flag, and README `:284` heads the broadband column "default". CLAUDE.md `:992` says "`--spectral` is opt-in".

**Explanation.** The code defaults to spectral, and the docs are stale (the change came in `9d6e986`, which is in HEAD's history).

The caveat matters for the review: the default demo database, `databases/S-76D_M3.nod`, has **no** `frequency`/`amplitude` variables (an `ncdump -h` grep finds none). So the default demo run falls back to broadband. That is why `golden/197421_LAE.csv` is still a broadband result despite the default.

Any shipped `*_fixed_load_factor.nod` database does carry spectra (B407: `frequency = 31`), so a run on one of those is spectral by default. "What the default does" therefore depends on the database chosen, which a review should state rather than say "spectral by default" alone.

## 2. phi from |up|: the upper hemisphere is never sampled — CONFIRMED

**Evidence**
- `src/sphere_database.cpp:711-723`, `phi_degrees(up, across)`: `std::atan(across / std::fabs(up))`. The range is [-90, 90], and the 0/0 case returns 90.
- Every lookup goes through it:
  - `lookup_angles` (`:726-732`), with `up = frame.normal.dot(direction)`;
  - `Database::level_from_geometry` (`:2062-2077`);
  - the spectral path's `lookup_angles(state.frame, path.direction)` (`:1570`), and `niceops.cpp:1515` (band-history angles).
- A grep finds no other phi computation in `src/`.

**Explanation.** A propagation direction with `normal·d > 0` (above the tip-path plane) and one with `normal·d < 0` and the same `along` and `across` components get identical (phi, theta). The lookup is a reflection through the TPP.

The sphere grid itself spans phi -180 to 170. The only time a cell beyond ±90 is touched is at phi = ±90 exactly, where `locate` may pick the [90, 100] cell with zero weight on the 100 knot. So the stored upper half carries no weight in any level.

For panam-built databases this loses nothing relative to the file. panam fills the upper half by mirroring through the horizontal plane (inventory §1), which in steady flight is the same reflection as NICEOPS's fold. A database carrying a *genuine* upper hemisphere (a computed sphere, say from a first-principles code) would have that half silently ignored.

**When `up > 0` actually arises.** The TPP normal is the load-factor unit vector (`trajectory.cpp:461-465`; `niceops.cpp:765`). In steady unaccelerated flight that is vertical, so "above the TPP" means above the aircraft's horizontal plane. Four configurations reach it:
- **Receiver above the source.** The microphone is the observer point raised by its height. Scattered observers carry their own z (`x,y,z[,height]`), and the default point-receiver height is 1.2 m. The mirror is read for:
  - the last frames of an approach to touchdown, or a ground-run or hover-taxi below the microphone height;
  - any rooftop or hillside observer in a `-b` file above a low aircraft.

  The ground model still assumes the flat plane z = 0 (`propagation.hpp:76-79`), so "terrain" exists only through observer z, and only in the direct path.
- **TPP tilt in manoeuvres.** A deceleration of a·g tilts the normal aft by about atan(a). Forward observers at depression angles smaller than that tilt are above the TPP. For example, 0.1 g gives 5.7°, which covers the 2-10° grazing band where the validation residuals concentrate. A banked turn puts outside-of-turn observers at depressions under the bank angle above the TPP. They read the mirrored (lower) level about the *tilted* plane. That is the QSAM convention, but it means no upper-hemisphere directivity (for example, BVI's upward lobe, or the rotor's in-plane thickness noise above the plane) can enter.
- **Refraction.** Under `--refraction`, the sphere is keyed on the ray's *launch* direction (`propagation.hpp:9-14`; `sphere_database.cpp:1570` uses `path.direction`). Downward refraction (downwind, inversions, ducts) launches long-range rays *above* the horizontal; `refraction.cpp:134-135` handles the "ray launched upward ... passes over its top" branch. Those rays read the mirror image below.

## 3. SEL window — CONFIRMED (code is first-to-last; CLAUDE.md is stale)

**Evidence**
- `src/util/sel.hpp:22-60`, `ten_db_down_interval`: scans forward from index 0 for the first sample `>= threshold`, then backward from the end for the last. Dips between them are included. The comment says it "replaces the D original's walk out from the peak, which stopped at the first sample below the threshold".
- `peak_range` (`:62-75`) and `sound_exposure_level` (`:154-166`) call it.
- The peak-search path `windowed_sound_exposure_level` uses the same function (`src/niceops.cpp:356-377`), with a `certainly_below` bound that skips states provably below the threshold.
- SEL = 10 log10(T Σ 10^(L/10) w_i) (`sel.hpp:139-149`).
- Doc side: CLAUDE.md `:289` and `:332` still say `peak_range` "walks out from the peak and stops at the **first** sample below peak - 10". CLAUDE.md also names `expand_window()`, which no longer exists anywhere in `src/` (grep).

**Explanation.** The current code is the standard 10 dB-down duration (the 14 CFR 36 EPNL duration convention). CLAUDE.md's sections "The SEL integral is oversampled..." and "Walking out from the peak..." describe the superseded D rule, including the "~0.7% of observers depend on one sample" fragility argument. A review should cite the code, not those sections.

## 4. `-D` is not normalised by 86 400 s — CONFIRMED

**Evidence**
- `src/util/sel.hpp:175-183`, `day_night_level`: adds 10 dB to night samples, then returns `sound_exposure_level` of the penalised history. That is again a 10 dB-down *windowed* SEL, located on the penalised history.
- `src/niceops.cpp:1801-1810`: per trajectory, `dnl_energy[index] += 10^(0.1·day_night_level(...))`.
- `src/niceops.cpp:2695-2702`: `dnl[i] = 10.0 * std::log10(dnl_energy[i])`, written directly.
- A grep for `86400` finds it only in `is_night` (`sel.hpp:129`, a modulo for time of day).
- README `:219`: "With one daytime operation the DNL equals the SEL exactly."

**Explanation.** The output is a night-penalised cumulative SEL (L_AE,dn, the sum of exposures), not DNL in the ANSI/FAA sense. A true DNL subtracts 10 log10(86 400) = 49.37 dB. A second departure: each operation contributes only its 10 dB-down window, not the whole-event exposure. The time step is still the full trajectory step.

## 5. Envelope clamping, 10% extension, nearest outside — CONFIRMED (with clarifications)

**Evidence**
- `src/sphere_database.cpp:1170-1203`: each axis is normalised by the database's min/max.
- `:1215-1218` and `:1247-1248`: queries are clamped with `clamp01`.
- `src/interpolator.cpp:54-78`: every condition with *any* coordinate within 10 ε of 0 or 1 is duplicated, with those coordinates moved to -0.1 or 1.1 (`extend_domain_by`, default 0.1, `interpolator.hpp:51`). The copy carries the same data index.
- `interpolator.cpp:113-114`: the query is clamped to [0, 1]^d again.
- `:128-136`: a barycentric weight of 1 or more snaps to that sphere.
- `:145-148`: outside the mesh, weight 1 on `nearest(clamped)` (Euclidean in normalised coordinates).
- `:87-93`: fewer than d+1 vertices means nearest only.

**Explanation.** "Clamped with a 10% extension" is right, but what is extended is the *triangulation*, not the query. Queries never leave the unit box. The duplicated boundary vertices carry the same sphere, so the blend near the envelope's faces is flat (constant extrapolation) rather than a true extension.

The unit box is per-axis min/max. A non-rectangular envelope (for example, high mu only at moderate gamma) leaves box corners outside the convex hull even after the extension. Queries there take the nearest condition in normalised space, which weights a 0.1 change in mu the same as a 0.1 change in gamma.

## 6. Load-factor treatment — PARTLY CONFIRMED

**Evidence**
- `fixed_load_factor` databases:
  - `src/sphere_database.cpp:1162-1169`: CT must be constant across spheres, and the triangulation drops to 2-D.
  - `:1221-1243`: `level_correction_db = CT > CT_ref ? 20 log10(CT/CT_ref) : 0`, with no cap. The comment says "there is no upper bound".
  - Applied in `eval` (`:1353`) and per band, uniformly, in `eval_spectrum` (`:1446-1447`).
- Other databases:
  - `:1244-1254`: CT is a third, normalised and clamped axis.
  - NICEOPS adds no LF term itself; the 20 log10(LF) lives in panam's materialised copies.
  - `S-76D_M3.nod` holds CT at {0, 1, 1.1, 1.2, 1.5, 2, 3, 5} × 0.005863 (164 conditions each, from `ncdump -v thrust_coefficient`). That is a clamp at exactly 5 g, plus a CT = 0 copy, so LF < 1 interpolates toward the 1 g level.
- CT itself is LF·W/(rho V_tip² pi R²) from the vehicle JSON (`trajectory.cpp:467`).

**Explanation.** The correction is 20 log10(CT/CT_ref), with CT_ref read from the database and CT computed from the JSON. It equals 20 log10(LF) only when the JSON's weight, density, tip speed and radius reproduce the database's reference. For B407 they do: the JSON gives CT = 0.0043348 against the database's 0.0043407, a ratio of 0.9986, so no correction at 1 g. A heavier JSON weight, or a lower density (a hot day), on a fixed_load_factor database adds level at *every* state through this term, silently.

"Clamp at 5 g" is the database's own CT maximum. It is 5 g for panam's LF list, but in general it is whatever the file holds. Neither policy ever reduces level below 1 g: fixed has no correction below CT_ref, and the CT-axis databases' 0-g copy carries the 1 g level.

## 7. No Doppler shift, no convective amplification; reception-time weight — CONFIRMED

**Evidence**
- A grep of `src/` for `doppler|convective|mach` finds only:
  - `sphere_database.cpp:2341-2342, 2384`: the AAM *export* writer, which writes `DOPPLER_SHIFT_REMOVED = 0`;
  - a tip-Mach comment.
- NICEOPS never *reads* a Doppler flag from a `.nod`. It would not notice a database whose shift had been removed.
- `src/util/sel.hpp:85-97`, `reception_factor`: max(1 + (dr/dt)/c, 1e-3), with dr/dt by central differences. It is applied as the per-sample weight in every SEL path (`sel.hpp:159-165`; `niceops.cpp:378-380`).
- Ranges are the straight-line range (`niceops.cpp:1779-1780`), or c·travel_time under refraction (`:1717, 1727, 1758`).

**Explanation.** No frequency is shifted before absorption, ground or A-weighting, and no (1 - M cos θ)^-n amplitude factor is applied. Both are taken as baked into each sphere at its measured speed. The only kinematic term is the reception-time Jacobian (Doppler time compression of the exposure).

`-i` time histories are on the emission-time index (`write_time_history_file(levels, integration_step, ...)`, `niceops.cpp:1812-1817`). `--band_history` adds a reception-time column (`:1820-1824`).

## 8. Interpolation domain — CONFIRMED: dB everywhere, energy only in the band sum

**Evidence**
- Broadband `eval` (`src/sphere_database.cpp:1281-1320`):
  - `blend_conditions` takes a weighted sum of dBA (and EAA) values across the four conditions (`:424-440`);
  - the four corners are then summed with the bilinear stencil weights, as a weighted sum of dB.
- Spectral `eval_spectrum` (`:1418-1448`): the same, per band, on `spectra_`. These are stored dB values: the netCDF `amplitude` is read as is (`:1066-1068`), and non-finite values are floored to `no_energy_floor_db = -1000`.
- `blend_present` (`:389-419`): renormalises over non-empty entries, so a gated cell drops out instead of averaging toward -1000 dB.
- Energy appears only in `weighted_level` (`src/weighting.cpp:45-61`), which takes Σ g_b·10^(L_b/10) over bands, and in the SEL sum.

**Explanation.** Both the barycentric blend across (mu, gamma, CT) and the within-sphere bilinear (phi, theta) interpolation are linear in dB (log-pressure), per band on the spectral path. This is not energy (p²) interpolation. Between dissimilar spheres, a dB average sits below the energy average (it is the geometric mean of the energies). Each band is interpolated independently, with no interpolation of spectral shape or frequency.

## 9. Refraction / ground effect on the broadband path — CONFIRMED refused

**Evidence** (all in `src/niceops.cpp`)
- `:2518-2523`: `--ground_effect` without spectral is an input error.
- `:2544-2549`: `--atmosphere` without spectral is an error.
- `:2550-2554`: `--refraction` other than `none` requires `--atmosphere`, so refraction is refused on broadband by transitivity.
- `:2529-2533`: `--propagation_table` needs spectral.
- `:2534-2538`: `--plate_table` needs `--ground_effect`.
- `:985-1001`: on a spectra-less database, any of these is an error rather than a silent fallback.

**Explanation.** The broadband path is spreading plus EAA·(r - R)/1000 only (`sphere_database.cpp:1343-1353`). Every advanced propagation option needs a spectral database, which closes the inventory's §5 gap.

Two error messages at `:2521` and `:2527` still say "the broadband path (--no-spectral or --db_cutoff)". That is stale: `--db_cutoff` no longer forces broadband (see C).

## 10. LAmax, EPNL/PNL, multiple sources — CONFIRMED absent

**Evidence**
- A grep of `src/` for `lamax|epnl|pnl|noy|tone_corr|loudness|c_weight|time.above|leq` finds only a comment citing the EPNL *duration* rule (`sel.hpp:25`).
- `-m` (`niceops.cpp:237-257`) writes the arithmetic mean of per-observer SEL in dB (not an energy mean), the maximum *SEL* over observers, and the fraction of observers ≥ 65 dB SEL.
- One source per aircraft: `database.source_state(Vec3(state.x, state.y, state.z), ...)` at `niceops.cpp:765`. A grep for `tail.rotor|source_offset` finds nothing. The JSON's `Tail Rotor` block (for example `databases/B407.json`) is never read (`trajectory.cpp:385-424` reads only Main Rotor radius and tip speed, Vehicle weight, and Atmosphere).

**Explanation.** Only A-weighted SEL (and the penalised sum of item 4) is output. LAmax can be post-processed from `-i` histories, but is not built in.

A crude two-source run is possible: two `-t` tracks with their own `-d` databases, summed with `-D`. It is not true source superposition, because each contribution is windowed to its own 10 dB-down duration and the output is the night-penalised sum.

---

## Further doc/code conflicts found

**A. README `:256` "Databases are not tracked in this repository" — REFUTED.** `git ls-files databases` lists 24 files: all the `.nod` databases and their JSON files.

**B. CLAUDE.md `:960` "No database in `databases/` carries spectra" — REFUTED.** `B407_fixed_load_factor.nod` has `frequency = 31`, `fixed_load_factor` and `shared_grid_and_frequency`; the other `*_fixed_load_factor.nod` databases are the spectral variants (inventory §1). The spectral path no longer needs the synthetic database from `test_spectral_propagation.cpp` to run.

**C. CLAUDE.md `:365`/`:992-996` and README `:1081` — REFUTED.** They say the peak-search path is skipped with `--spectral`, and that `--spectral` is rejected with `--db_cutoff`.
- `niceops.cpp:1445`: `expand_from_peak = direction_bounds.usable() && !culling && !needs_full_history`, with no spectral term.
- The comment at `:1431-1440` lists only `-a`, `-i`, `--band_history` and `-D` as reasons to march.
- `:1318-1321` and `:1396-1403`: `--db_cutoff` on spectral culls per direction (`cull_by_direction = culling && spectral`) on `spectral_direction_bounds`.
- README `:1065` agrees with the code. The inventory did not flag this.

**D. Asserts are compiled out in Release — CONFIRMED.** `CMakeLists.txt:16-23` defaults the build type to Release and only *comments* "do NOT add -DNDEBUG". CMake's own default `CMAKE_CXX_FLAGS_RELEASE` is `-O3 -DNDEBUG`, verified by configuring an empty project in scratch, and nothing in `CMakeLists.txt` removes it. So the `assert(std::isfinite(...))` guards (for example `sphere_database.cpp:1346-1349, 2075-2077`) are inactive in the documented build. This contradicts CLAUDE.md's "asserts live" port contract.

**E. Projection naming — the code is the American polyconic.** `geodetic.hpp:8`, CLAUDE.md `:997` and README `:1181` call it a transverse Mercator. But `geodetic.cpp:68-85` computes:
- ms = N cos φ / sin φ = N cot φ;
- E = N cot φ · sin(Δλ sin φ);
- N = M(φ) - M(φ0) + N cot φ (1 - cos(Δλ sin φ)).

That is the ellipsoidal polyconic (Snyder eqs. 18-2/18-3 form). The code review is right and the docs are wrong. It is immaterial over a few km, but it should not be cited as TM.

**Not re-verified here:** the evaluated-state fraction ("about 94%" against the docs' "70%, 1.2x"), which would need a built binary and a timed run, and the inventory's validation numbers (§10), which come from reports rather than code.

**Inward spreading has no floor — confirmed.** The inventory flags this as a known issue, not a doc conflict: `sphere_database.cpp:1344` (broadband) and `:1415` (spectral) take 20 log10(r/R) with no floor. Inside the sphere, `excess_factor·(r - R)` also turns negative, adding level on the broadband path.
