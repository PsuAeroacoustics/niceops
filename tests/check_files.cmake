# Checks that files exist or do not, for CTest cases whose subject is what a run
# wrote.  cmake -DFILES="a;b" -DABSENT="c;d" -P tests/check_files.cmake, from the
# source root.
foreach(file IN LISTS FILES)
  if(NOT EXISTS "${file}")
    message(FATAL_ERROR "expected ${file} to exist")
  endif()
endforeach()
foreach(file IN LISTS ABSENT)
  if(EXISTS "${file}")
    message(FATAL_ERROR "expected ${file} not to exist")
  endif()
endforeach()
