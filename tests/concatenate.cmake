# Concatenate files: cmake -DOUT=<file> "-DINPUTS=<a;b;...>" -P concatenate.cmake
# CTest uses it to assemble the expected output of a mixed-receiver run from
# runs over each of its subsets.
file(WRITE "${OUT}" "")
foreach(input ${INPUTS})
    file(READ "${input}" contents)
    file(APPEND "${OUT}" "${contents}")
endforeach()
