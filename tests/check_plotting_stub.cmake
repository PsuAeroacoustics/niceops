# Every plt:: function plotting.cpp defines, plotting_disabled.cpp defines too.
#
# Only the -DNICEOPS_PLOTTING=OFF build links the stub, and CI builds with plotting, so a plt:: entry
# point added to plotting.cpp alone would break that build's link without a red check anywhere.  The
# plotting build compiles the stub as an object library, and this compares the two objects' defined
# external symbols.
#
#   cmake -DNM=<nm> -DCORE=<libniceops_core.a> -DSTUB=<plotting_disabled object> -P check_plotting_stub.cmake

function(plt_symbols out file)
  execute_process(COMMAND ${NM} -C --defined-only ${file}
    OUTPUT_VARIABLE listing RESULT_VARIABLE status ERROR_VARIABLE err)
  if(NOT status EQUAL 0)
    message(FATAL_ERROR "check_plotting_stub: ${NM} failed on ${file}: ${err}")
  endif()
  # An archive's listing heads each member with "member:" (GNU) or "<archive>(member):" (llvm); keep
  # the plotting.cpp member's lines only.
  string(REPLACE ";" "\\;" listing "${listing}")
  string(REPLACE "\n" ";" lines "${listing}")
  set(in_member TRUE)
  if(ARGN)
    set(in_member FALSE)
  endif()
  set(symbols)
  foreach(line IN LISTS lines)
    if(ARGN AND line MATCHES "^(.*\\()?([^()/]+)\\)?:$")
      if(CMAKE_MATCH_2 STREQUAL ARGN)
        set(in_member TRUE)
      else()
        set(in_member FALSE)
      endif()
    elseif(in_member AND line MATCHES " T _?(plt::.*)$")
      list(APPEND symbols "${CMAKE_MATCH_1}")
    endif()
  endforeach()
  list(REMOVE_DUPLICATES symbols)
  set(${out} "${symbols}" PARENT_SCOPE)
endfunction()

plt_symbols(real ${CORE} plotting.cpp.o)
plt_symbols(stub ${STUB})
if(NOT real)
  message(FATAL_ERROR "check_plotting_stub: no plt:: symbols found for plotting.cpp.o in ${CORE}")
endif()
set(missing)
foreach(symbol IN LISTS real)
  if(NOT symbol IN_LIST stub)
    list(APPEND missing "${symbol}")
  endif()
endforeach()
if(missing)
  list(JOIN missing "\n  " missing)
  message(FATAL_ERROR "check_plotting_stub: plotting_disabled.cpp lacks\n  ${missing}")
endif()
list(LENGTH real count)
message(STATUS "check_plotting_stub: plotting_disabled.cpp defines all ${count} plt:: functions")
