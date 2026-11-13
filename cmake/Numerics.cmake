# Floating-point policy for the production library and numerical header tests.
# Apply explicitly only to the targets already using these flags. Research and
# portable-integer definitions remain directory-wide in the root CMakeLists.
function(certpolysolve_strict_floating_point target)
    if(MSVC)
        target_compile_options(${target} PRIVATE /fp:strict)
    else()
        # Contract a multiply/add only where the source explicitly uses std::fma.
        target_compile_options(${target} PRIVATE -fno-fast-math -ffp-contract=off)
    endif()
endfunction()
