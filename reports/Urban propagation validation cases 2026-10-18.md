# Urban propagation validation cases, 2026-10-18

> **What this is.** A search for data to validate the urban receivers and the point-source model,
> building on `research_notes/Urban canyon propagation methods/validation_data_and_inputs.md` and the
> plan's validation section. Four research passes ran in parallel: full-scale aircraft among
> buildings; drones and eVTOL; scale models and controlled sources; benchmarks, the turbulence
> floor and impulsive sources. Facts marked (U) come from abstracts or snippets only; everything
> else was read in the source.
>
> **Since then.** EDGE is in the harness (`urban/edge/`) and scored in
> [FIDELITY](../docs/FIDELITY.md); its UE06 and UE07 measurements may be swapped against their
> geometry ([Next Steps](NICEOPS%20Urban%20Propagation%20Next%20Steps.md#validation-still-open)).
> Of the data requests at the end, NLR's is under way.

## Bottom line

- No dataset has a rotorcraft at 30-1500 m heard in a real street canyon with band data,
  geometry and meteorology. NASA/TP-20250011527 (May 2026, Sec. 3.3.2) names this as a gap and
  recommends canonical geometries with simple sources. The proposed campaign in the plan remains
  the only route to that case.
- Two new downloadable sets are worth starting on now: RWTH's EDGE (fully resolved below) and
  Wang et al. 2026, a 1:40 scale model with a moving elevated source over houses.
- The best real-aircraft data (Schiphol container courtyards, Paris heliport rooftops, NLR
  DigiCity with drone flights, PoliMi rotor over a box) are all "on request". Six emails would
  unlock most of the ladder (last section).
- Three entries in the plan's dataset table are wrong (see Corrections).

## Ladder by model component

| Component | Best case | Second case | Data |
| --- | --- | --- | --- |
| Roof-edge diffraction, shadow order | EDGE UE00-07 (terraced rows, six roof shapes, drone source) | PoliMi rotor over box | EDGE public; PoliMi on request |
| Vertical-corner diffraction | EDGE UE20-22 (crossing, dense block, street corner) | DigiCity corner mics 3, 7 and two dedicated corner mics | EDGE public; DigiCity on request |
| Double diffraction, 3-D | EDGE UE14-19 (detached house, 1/12 and 1/5) | EDB2 BTM references you build (GPL, out of process) | Public |
| Facade images, canyon | EDGE UE12/13 (canyon, one side raised) | DigiCity reflection-only mics 2, 8, 10; FAA TC-18/15 wall | Public |
| Diffuse facade scattering s(f) | EDGE UE23-26, UE28/30 (periodic MDF scatterers) + Heimes interlab coefficients for the same scatterers | SonicBOBS: diffuse reflection did not help for an isolated building | Public |
| Turbulence floor | Forssén et al. 2014 model, fully specified | Forssén & Ögren 2002 container barrier with measured C_v^2, C_T^2 | In papers (open PDFs) |
| Low-band canyon modes | Dragna et al. 2022 2-D Euler canyons (numerical) | Richoux et al. LAUM scale street (U) | In papers |
| Point-source mode, drone | DigiCity hovers and flyovers (8.5 and 20 m) | Wang et al. 2026 moving source | DigiCity on request; Wang public |
| Moving elevated source, whole chain | Lugten Schiphol courtyards (2383 flyovers) | MOTUS Paris heliport (15,329 helicopter events) | On request |

## Downloadable now

**EDGE, RWTH Aachen** (Palenda & Vorländer, Zenodo 10.5281/zenodo.21264202, CC BY-SA 4.0,
2026-08-19). Score 5.
- Hemi-anechoic room, rigid floor, lower limit 100 Hz. Scales 1/12 (UE00-07, 12-15, 20-26),
  1/5 (UE08-09, 16-19), 1/3 (UE10-11), 1/20 (UE27-30).
- Scenes: terraced rows with flat, gable, half-gable, pent, mansard and shed roofs; a canyon
  (UE12/13); a detached house; crossing streets; a dense block; a 23 m deep canyon with periodic
  facade scatterers; a model of the IHTA Park buildings (9.6 and 24.7 m full scale).
- Elevated "drone" source in UE00, UE12-14, UE16, UE20 and UE27, at 15-30 m full scale. It sits
  about 3x roof height but at 20-24 deg elevation from the receivers (arithmetic from the SOFA
  positions). No source is overhead. UE21-23 have ground sources only.
- Calibrated IRs in SOFA, 96 kHz, with temperature and humidity; loudspeaker directivity on a
  1 deg grid; geometry in FreeCAD, STEP and Blender with per-face materials; carpet impedance fits.
- Usable band is not stated. Inferred from the drivers: about 170 Hz-3.5 kHz full scale at 1/12,
  100 Hz-2 kHz at 1/20, 40 Hz-4 kHz at 1/5 (U).
- Run NICEOPS at model scale, with rigid walls and the measured directivity.
- Goldens built from EDGE geometry would carry CC BY-SA on those files.

**Wang, Anai, Yano & Sakamoto 2026** (Acoustical Science and Technology 47(3):209-220,
10.1250/ast.e25.73, CC BY-ND; data Zenodo 10.5281/zenodo.15833542, CC BY 4.0). Score 4.
- 1/40 scale, University of Tokyo hemi-anechoic room. An air-jet point source moves along a
  100 m full-scale line at 3, 5, 10 and 15 m height, 85 positions per pass.
- Detached houses 8x8 and 8x16 m, five layouts (16.8-33.6 % density), heights 4, 7 and 10 m;
  12 ground receivers per layout.
- Building insertion loss per source position, octave bands 125 Hz-2 kHz full scale; about 4300
  receiver-condition pairs.
- Layouts are only in Fig. A1 and must be digitized. Elevation angles are low.

**Heimes et al. interlaboratory scattering database** (Acta Acustica 10, 65 (2026),
10.1051/aacus/2026060; Zenodo 17660108 coefficients, 17660089 raw IRs, 17660055 BEM; CC BY 4.0).
Score 3, for s(f) only.
- Measured surfaces include a 1:40 facade bay projection; BEM covers the IHTA Park facades, and
  probably the EDGE scatterers (U). No measured balcony, window or corrugated container facade
  exists in any open database.

**Turbulence floor, open papers.**
- Forssén et al. 2014 (Acta Acustica 100(5):810-815; author PDF at Chalmers 202744). The full
  regression is given: b1 = -52.8 + 10 lg C_v^2 (b2 = 11.3, b3 = -17.1) and
  b1 = -49.6 + 10 lg C_T^2 (b2 = 11.5, b3 = -13.1), the (10/3) lg(f/1 kHz) term, canyon terms
  gamma1 = 7 or 14 dB and gamma2 = 2H/W. Check the sign of gamma1 against the published page.
  Score 5 as a specification.
- Forssén & Ögren 2002 (Applied Acoustics 63(2):173-187; Chalmers 23206). Two containers as a
  barrier 2.55 m high, source 8 m in front, receivers 4-18 m behind, sonic anemometer. Calm
  C_v^2 about 0.4 and windy about 15 m^(4/3)/s^2 (units to confirm). Turbulence raises the
  shadow level about 5 dB at 4 kHz and 10 dB at 8 kHz, little below 2 kHz. Score 4.

**Reference builders.** Svensson's EDB2 (GPL v3, exact BTM to 6th order for IRs) builds 3-D
finite-building references with roof and corner edges together. This extends the 2-D BEM check
in FIDELITY. Svensson, Brick & Forssén, Forum Acusticum 2014 gives cube and building cases (full
text not found).

## Available on request

| Case | What it adds | Score | Contact |
| --- | --- | --- | --- |
| Lugten et al. 2024, Building and Environment 255:111400, Schiphol container courtyards | 10 facade mics (0.2 m off, 1.5 or 3.9 m up), 1/3-octave every 0.125 s plus WAV, radar tracks, KNMI weather, 2383 flyovers; 1:50 scale model of 32 no-wind events (FA2023) | 4 | m.c.lugten@tudelft.nl; data statement: "on request" |
| MOTUS, Paris-Issy heliport (Caillet et al. ERF 2024; CEAS 2026, 10.1007/s13272-026-00957-z) | 15 rooftop stations 14-40 m up with exact coordinates; LAeq and LCeq every 0.1 s; arrival azimuth and elevation every 0.1 s; DGAC radar; 15,329 events of EC135, EC145, AS355, A109, AW139, R44 | 4 | Bruitparif (Sineau, Hellot); Airbus Helicopters (Caillet, Dieumegard) |
| NLR DigiCity (ERF2025-518 and -009) | 24 PCB mics at 1.2 m; M600 flyovers at 8.5 and 20 m, 1 and 10 m/s; hovers at 7 spots with rotor RPM and attitude logged; omni at spots 2 and 5. NLR's own model is 10-15 dB low at mic 7 | 4 | Kylie.Knepper@nlr.nl; project lead Timo Nijman |
| PoliMi rotor over a box (Gibertini et al., Aerospace 10:701, 2023) | 1/100 rotor at 5 descent points over the GARTEUR DLR box (45 m building full scale); floor grid 5-50 cm from the box; 20-35 dB shielding | 3-4 | giuseppe.gibertini@polimi.it, zanotti@aero.polimi.it ("public on request") |
| NASA SonicBOBS 2009 and Edwards houses 2006/07 | F/A-18 booms on an office building and houses; Riegel & Sparrow JASA 145:2317 (2019) compare ray tracing and radiosity | 3-4 | Loubeau, Klos (NASA LaRC) |
| Empa Zurich medical drone routes (Affolter et al., Quiet Drones 2026) | Six months, 6 stations at 4 m, three on buildings; sonAIR with reflections within about 0.6 dB at quiet stations | 3 | claudio.affolter@empa.ch |
| CRREL training village (Albert & Liu 2010; Mehra et al. 2014) | 15 buildings, 14 receivers, ground-level C4, ARD within about 3 dB to 450 Hz | 3 | D. Albert (ERDC-CRREL) |
| Heli-STAR Vol. 5 (DOT/FAA/ND-97/13, DTIC ADA370872) | S-76 approaches to a ground helipad at Norcross, GA; mics on a roof, at a building's base and on an open slope; GPS tracks | 3 | Retry DTIC (403 today); Ahuja, Funk (GTRI) for raw data |
| Van Renterghem et al. 2013 Ghent shielded facades | Long-term 1/3-octave at exposed and shielded facades; up to 10 dB between weak and strong turbulence classes above 1 kHz | 3 | Timothy Van Renterghem (UGent) |

## Corrections to the plan's dataset table

- **Hao & Kang 2014** (Applied Acoustics 84:73-82) is a CadnaA simulation, not a measurement.
  Drop it as data.
- **ACRP 02-79 "566 SEL events at four airports"** are the project's own fixed-wing measurements
  at LAX and LGB (A321, B737, A320), 1 s LAeq and LAmax only. The "NYC Urban Helicopter" data is
  Volpe's July 2000 AS350 campaign for the FAA 2004 Report to Congress, mostly over water at
  Liberty State Park, with one microphone near a facade (AAM 1.8-3.3 dB low there). Its tables
  are in a separately posted Appendix G. Score 2.
- **Heli-STAR Vol. 5** was an S-76 at a ground-level helipad (Vol. 1, Sec. 8.2.4), not a rooftop
  landing as Gibertini et al. state.
- **Donavan**: the MIT thesis on DSpace is the 1976 ScD "Sound propagation in urban spaces"
  (1721.1/27489); the V/STOL title is a 1974 JASA paper with Lyon (U).
- **Lugten 2024 Table 2**: its LAmax rows for 125 Hz-1 kHz are identical, an apparent copy error.

## Gaps nobody fills

- A rotorcraft at 30-1500 m over real canyons, in bands, with tracks and weather.
- Measured low-band canyon modes with an elevated source.
- Measured scattering coefficients of real facades (balconies, windows, corrugated containers).
- Turbulence scattering into shadows from an elevated source. All turbulence data use ground or
  road sources.

## Suggested order

1. EDGE UE00-07 and UE12/13 at model scale: roof diffraction and canyon images, with no requests
   needed. Then UE20-22 for corners and UE23-26 for s(f).
2. Wang et al. 2026, after digitizing Fig. A1: moving elevated source, insertion loss by band.
3. Implement Forssén 2014 as a reference for the turbulence floor, and check it on Forssén &
   Ögren 2002.
4. Send the requests below now, since replies take weeks.
5. Lugten, MOTUS and DigiCity drone flights as whole-chain tests when the data arrives.

## Requests to send

1. Martijn Lugten (TU Delft): field-lab 1/3-octave events, tracks, container layout, and the
   1:50 model data.
2. Kylie Knepper (NLR): DigiCity layout and mic coordinates, omni calibration, drone flight data.
   NLR and TU/e are recruiting a PhD on urban sound modelling validated at DigiCity, a possible
   joint channel.
3. Bruitparif and Airbus Helicopters: MOTUS event tables, arrival angles and radar tracks.
4. Gibertini and Zanotti (PoliMi): the full floor-grid database.
5. Loubeau and Klos (NASA LaRC): SonicBOBS waveforms and building geometry.
6. Juliet Page (Volpe): FAA 2004 Appendix G and the NYC microphone positions.
