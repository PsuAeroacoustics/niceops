# Checks that files are byte for byte the same, pair by pair, for CTest cases
# that hold several outputs of two runs equal.  cmake -DPAIRS="a1;b1;a2;b2" -P
# tests/compare_pairs.cmake, from the source root; it names the first pair
# that differs.
list(LENGTH PAIRS count)
math(EXPR odd "${count} % 2")
if(count EQUAL 0 OR odd)
  message(FATAL_ERROR "PAIRS needs an even, nonzero number of files, got ${count}")
endif()
math(EXPR last "${count} - 1")
foreach(i RANGE 0 ${last} 2)
  math(EXPR j "${i} + 1")
  list(GET PAIRS ${i} a)
  list(GET PAIRS ${j} b)
  execute_process(COMMAND ${CMAKE_COMMAND} -E compare_files "${a}" "${b}" RESULT_VARIABLE differs)
  if(differs)
    message(FATAL_ERROR "${a} and ${b} differ")
  endif()
endforeach()
message(STATUS "${count} files agree pair by pair")
