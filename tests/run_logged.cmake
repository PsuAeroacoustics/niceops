# Run a command with a fresh log file and print the log: cmake -DLOG=<file>
# "-DCOMMAND=<program;arg;...>" -P run_logged.cmake.  The program is given
# -L info -l <file> by the caller.  The logger appends, so the file is removed
# first; info lines reach only the file, so a CTest that matches on them
# (PASS/FAIL_REGULAR_EXPRESSION) needs them printed.  A failing command fails.
file(REMOVE "${LOG}")
execute_process(COMMAND ${COMMAND} RESULT_VARIABLE result)
if(EXISTS "${LOG}")
    file(READ "${LOG}" contents)
    message("${contents}")
endif()
if(NOT result EQUAL 0)
    message(FATAL_ERROR "exited with ${result}")
endif()
