# Exercise both acceptance and rejection, and check that the CSV contains
# actual between-node errors rather than comparisons at the table's nodes.
execute_process(COMMAND "${TOOL}" --method ground_model --units m --bands 250
  --source_heights 10 20 --receiver_heights 1.2 --max_range 100 --step 10
  --validation_csv "${OUT}/validation.csv" --validation_max_error 100
  -o "${OUT}/validation.nc" RESULT_VARIABLE status ERROR_VARIABLE summary)
if(NOT status EQUAL 0)
  message(FATAL_ERROR "Validation unexpectedly failed: ${summary}")
endif()
file(READ "${OUT}/validation.csv" rows)
if(NOT rows MATCHES "0,15,1\\.2[0-9]*,15,250,.*covered")
  message(FATAL_ERROR "No between-node source/range validation row: ${rows}")
endif()
execute_process(COMMAND "${TOOL}" --method ground_model --units m --bands 250
  --source_heights 10 20 --receiver_heights 1.2 --max_range 100 --step 10
  --validation_csv "${OUT}/rejected.csv" --validation_max_error 0
  -o "${OUT}/rejected.nc" RESULT_VARIABLE status ERROR_VARIABLE summary)
if(NOT status EQUAL 2)
  message(FATAL_ERROR "Expected tolerance rejection (exit 2), got ${status}: ${summary}")
endif()
file(WRITE "${OUT}/uniform.csv" "z_m,T_K,RH,pressure_kpa\n0,280,20,101.325\n1000,280,20,101.325\n")
execute_process(COMMAND "${TOOL}" --method pe --units m --bands 125
  --atmosphere "${OUT}/uniform.csv" --source_heights 2 2.2 --receiver_heights 1.2
  --max_range 100 --step 1 --samples 21 --validation_samples 41
  --validation_csv "${OUT}/pe.csv" --validation_max_error 2
  -o "${OUT}/pe.nc" RESULT_VARIABLE status ERROR_VARIABLE summary)
if(NOT status EQUAL 0)
  message(FATAL_ERROR "PE midpoint validation failed: ${summary}")
endif()
file(READ "${OUT}/pe.csv" rows)
if(NOT rows MATCHES "covered")
  message(FATAL_ERROR "PE reference had no covered rows")
endif()
