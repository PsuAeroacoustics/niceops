# The PE tools' thread count changes the wall time and nothing else: tables
# built with the default (NICEOPS_THREADS unset), with --threads 2 and with
# NICEOPS_THREADS=3 hold the same bits, still and turbulent, PE and ground
# model, and niceops_pe_compare writes the same CSV.  A NICEOPS_THREADS that is
# not a positive count is refused rather than ignored.
file(MAKE_DIRECTORY "${OUT}")
file(WRITE "${OUT}/windy.csv"
  "z_m,T_K,RH,wind_speed_mps,wind_direction\n0,285,30,0,270\n10,285,30,4,270\n1000,280,30,8,270\n")
set(common --method pe --units m --bands 125 250 --atmosphere "${OUT}/windy.csv" --bearings 2
  --source_heights 2 20 --receiver_heights 0 1.2 --max_range 200 --step 50 --samples 2)
set(unset_env ${CMAKE_COMMAND} -E env --unset=NICEOPS_THREADS)

macro(build name)
  execute_process(COMMAND ${ARGN} RESULT_VARIABLE status ERROR_VARIABLE summary)
  if(NOT status EQUAL 0 OR NOT summary MATCHES "Wrote")
    message(FATAL_ERROR "${name} failed (${status}): ${summary}")
  endif()
endmacro()
macro(expect_identical a b)
  execute_process(COMMAND "${COMPARE}" "${OUT}/${a}" "${OUT}/${b}" RESULT_VARIABLE status
    OUTPUT_VARIABLE said ERROR_VARIABLE said)
  if(NOT status EQUAL 0)
    message(FATAL_ERROR "${a} and ${b} differ: ${said}")
  endif()
  message(STATUS "${a} = ${b}: ${said}")
endmacro()

build(default ${unset_env} "${TOOL}" ${common} -o "${OUT}/default.nc")
build(threads2 ${unset_env} "${TOOL}" ${common} --threads 2 -o "${OUT}/threads2.nc")
if(NOT summary MATCHES " on 2 threads: ")
  message(FATAL_ERROR "--threads 2 did not reach the PE: ${summary}")
endif()
build(environment ${CMAKE_COMMAND} -E env NICEOPS_THREADS=3 "${TOOL}" ${common} -o "${OUT}/environment.nc")
if(NOT summary MATCHES " on 3 threads: ")
  message(FATAL_ERROR "NICEOPS_THREADS=3 did not reach the PE: ${summary}")
endif()
expect_identical(default.nc threads2.nc)
expect_identical(default.nc environment.nc)

set(turbulence --gamma_t 3e-6 --realizations 2 --seed 5)
build(turbulent ${unset_env} "${TOOL}" ${common} ${turbulence} -o "${OUT}/turbulent.nc")
build(turbulent1 ${unset_env} "${TOOL}" ${common} ${turbulence} -j 1 -o "${OUT}/turbulent1.nc")
expect_identical(turbulent.nc turbulent1.nc)

set(ground --method ground_model --units m --bands 125 250 --source_heights 2 20 --receiver_heights 0 1.2
  --max_range 200 --step 50)
build(ground ${unset_env} "${TOOL}" ${ground} -o "${OUT}/ground.nc")
build(ground2 ${unset_env} "${TOOL}" ${ground} --threads 2 -o "${OUT}/ground2.nc")
expect_identical(ground.nc ground2.nc)

# The environment's value is checked, and --threads wins over it.
execute_process(COMMAND ${CMAKE_COMMAND} -E env NICEOPS_THREADS=many "${TOOL}" ${common} -o "${OUT}/refused.nc"
  RESULT_VARIABLE status ERROR_VARIABLE summary)
if(status EQUAL 0 OR NOT summary MATCHES "NICEOPS_THREADS='many' is not a positive thread count")
  message(FATAL_ERROR "expected NICEOPS_THREADS=many refused, got ${status}: ${summary}")
endif()
build(override ${CMAKE_COMMAND} -E env NICEOPS_THREADS=many "${TOOL}" ${common} --threads 2 -o "${OUT}/override.nc")
expect_identical(default.nc override.nc)

# niceops_pe_compare: the PE's bands on one thread or on the default.
set(path --atmosphere "${OUT}/windy.csv" --units m --toward 90 --source_height 20 --receiver_height 1.2
  --bands 125 250 --max_range 300 --step 50 --samples 2 --shadow_model resolved)
foreach(run default single)
  if(run STREQUAL "single")
    set(extra --threads 1)
  else()
    set(extra)
  endif()
  execute_process(COMMAND ${unset_env} "${COMPARE_TOOL}" ${path} ${extra} RESULT_VARIABLE status
    OUTPUT_VARIABLE csv_${run} ERROR_VARIABLE summary)
  if(NOT status EQUAL 0 OR NOT csv_${run} MATCHES "pe_125")
    message(FATAL_ERROR "niceops_pe_compare (${run}) failed (${status}): ${summary}")
  endif()
endforeach()
if(NOT csv_default STREQUAL csv_single)
  message(FATAL_ERROR "niceops_pe_compare differs with --threads 1:\n${csv_default}\n---\n${csv_single}")
endif()
