# Writes OUTPUT, the header src/build_info.cpp includes, with the source tree's
# git state.  Run at every build (the niceops_build_info target), and rewrites
# the header only when its text changes, so an unchanged tree recompiles
# nothing.  A tree with no git, such as an unpacked archive, records "unknown".
execute_process(
  COMMAND git describe --always --dirty --tags
  WORKING_DIRECTORY "${SOURCE_DIR}"
  OUTPUT_VARIABLE NICEOPS_GIT_DESCRIBE
  OUTPUT_STRIP_TRAILING_WHITESPACE
  ERROR_QUIET
  RESULT_VARIABLE NICEOPS_GIT_RESULT)
if(NOT NICEOPS_GIT_RESULT EQUAL 0 OR NICEOPS_GIT_DESCRIBE STREQUAL "")
  set(NICEOPS_GIT_DESCRIBE "unknown")
endif()
set(NICEOPS_BUILD_INFO_TEXT "#pragma once
#define NICEOPS_GIT_DESCRIBE \"${NICEOPS_GIT_DESCRIBE}\"
")
if(EXISTS "${OUTPUT}")
  file(READ "${OUTPUT}" NICEOPS_OLD_TEXT)
else()
  set(NICEOPS_OLD_TEXT "")
endif()
if(NOT NICEOPS_OLD_TEXT STREQUAL NICEOPS_BUILD_INFO_TEXT)
  file(WRITE "${OUTPUT}" "${NICEOPS_BUILD_INFO_TEXT}")
endif()
