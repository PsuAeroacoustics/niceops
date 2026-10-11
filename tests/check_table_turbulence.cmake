# niceops_propagation_table's reciprocity and turbulence options, end to end:
# both marches and a turbulent table with a height profile build, and the
# combinations that make no sense are refused with their reason.
file(MAKE_DIRECTORY "${OUT}")
file(WRITE "${OUT}/uniform.csv" "z_m,T_K,RH,pressure_kpa\n0,280,20,101.325\n1000,280,20,101.325\n")
file(WRITE "${OUT}/gamma.csv" "z_m,gamma\n0,1e-5\n50,1e-6\n")
set(common --method pe --units m --bands 250 --atmosphere "${OUT}/uniform.csv"
  --source_heights 2 20 --receiver_heights 0 1.2 --max_range 200 --step 50 --samples 1)
foreach(march receiver source)
  execute_process(COMMAND "${TOOL}" ${common} --march ${march} -o "${OUT}/${march}.nc"
    RESULT_VARIABLE status ERROR_VARIABLE summary)
  if(NOT status EQUAL 0 OR NOT summary MATCHES "Wrote")
    message(FATAL_ERROR "--march ${march} failed (${status}): ${summary}")
  endif()
endforeach()
execute_process(COMMAND "${TOOL}" ${common} --gamma_profile "${OUT}/gamma.csv" --realizations 2 --seed 3
  --outer_scale 5 -o "${OUT}/turbulent.nc" RESULT_VARIABLE status ERROR_VARIABLE summary)
if(NOT status EQUAL 0 OR NOT summary MATCHES "Wrote")
  message(FATAL_ERROR "turbulent table failed (${status}): ${summary}")
endif()
macro(expect_refusal pattern)
  execute_process(COMMAND "${TOOL}" ${ARGN} RESULT_VARIABLE status ERROR_VARIABLE summary)
  if(status EQUAL 0 OR NOT summary MATCHES "${pattern}")
    message(FATAL_ERROR "expected a refusal matching '${pattern}', got ${status}: ${summary}")
  endif()
endmacro()
# Under ground_model, --gamma_t is the HARMONOISE coherence; the phase screens' own options are PE's,
# and Chessell's --turbulence and --roughness are the ground model's.
expect_refusal("apply to --method pe only" --method ground_model --units m --bands 250 --source_heights 2 20
  --receiver_heights 1.2 --max_range 100 --step 50 --gamma_t 1e-6 --realizations 2 -o "${OUT}/refused.nc")
expect_refusal("apply to --method ground_model only" ${common} --turbulence 4e-4 -o "${OUT}/refused.nc")
expect_refusal("checks a still table" ${common} --gamma_t 1e-6 --validation_csv "${OUT}/v.csv" -o "${OUT}/refused.nc")
file(WRITE "${OUT}/zero.csv" "z_m,gamma\n0,0\n50,0\n")
expect_refusal("zero everywhere" ${common} --gamma_profile "${OUT}/zero.csv" -o "${OUT}/refused.nc")
