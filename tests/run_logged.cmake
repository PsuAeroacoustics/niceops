# Run a command with a fresh log file and print the log: cmake -DLOG=<file>
# "-DCOMMAND=<program;arg;...>" -P run_logged.cmake.  The program is given
# -L info -l <file> by the caller.  The logger appends, so the file is removed
# first; info lines reach only the file, so a CTest that matches on them
# (PASS/FAIL_REGULAR_EXPRESSION) needs them printed.  CTest ignores the exit
# status of a test with PASS_REGULAR_EXPRESSION, so a failing command prints a
# line that the caller's FAIL_REGULAR_EXPRESSION matches (LOGGED_RUN_FAILED in
# CMakeLists.txt) as well as failing.
file(REMOVE "${LOG}")
execute_process(COMMAND ${COMMAND} RESULT_VARIABLE result)
if(EXISTS "${LOG}")
    file(READ "${LOG}" contents)
    message("${contents}")
endif()
if(NOT result EQUAL 0)
    message(FATAL_ERROR "run_logged: the command exited with ${result}")
endif()
