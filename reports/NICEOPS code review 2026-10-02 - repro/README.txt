Reproduction scripts for "NICEOPS code review 2026-10-02.md" (section 6).

They were run against a scratch copy of HEAD 5d761b8, never inside the repository:

  git archive HEAD | tar -x -C <scratch>/tree          # plus symlinks to databases/*.nod if untracked
  cmake -S <scratch>/tree -B <scratch>/tree/build -DCMAKE_BUILD_TYPE=Release -DNICEOPS_PLOTTING=OFF
  cmake --build <scratch>/tree/build -j8
  # a run directory laid out like the repo root: symlinks to tracks/ databases/ atmospheres/ plates/
  # arrays/ golden/ and tests/data, plus writable state/ integrated/ log/
  cd <scratch>/run && bash exp1.sh   # paths in the scripts assume ../tree/build/niceops

exp1.sh  golden_compare blind to NaN; -i overwrite across trajectories; exception at -j 1 vs -j 4
exp2.sh  flat/descending time, NaN coordinate, --trim past the end, grid count under --units, invisible errors
exp3.sh  one NaN x in 6551 rows; grid columns; DNL on a clockless geodetic track
pole2.cpp  DirectionBounds near the nose and tail axes (link against build/libniceops_core.a; see the header of the report)
anti.cpp   LocalFrame across the antimeridian
