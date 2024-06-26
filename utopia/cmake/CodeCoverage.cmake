add_library(coverage_config INTERFACE)
install(TARGETS coverage_config EXPORT UtopiaTargets)

if(UTOPIA_ENABLE_CODE_COVERAGE AND CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
    message(STATUS "[STATUS] code coverage enabled")
    # Add required flags (GCC & LLVM/Clang)

    if(UTOPIA_ENABLE_CODE_COVERAGE AND CMAKE_CXX_COMPILER_ID MATCHES "Clang")
        target_compile_options(
            coverage_config
            INTERFACE -O0 # no optimization
                      -g # generate debug info
                      --coverage # sets all required flags
                      -fprofile-instr-generate
                      -fcoverage-mapping
                      -ftest-coverage)
    else()
        target_compile_options(
            coverage_config
            INTERFACE -O0 # no optimization
                      -g # generate debug info
                      --coverage # sets all required flags
                      -ftest-coverage)
    endif()

    if(CMAKE_VERSION VERSION_GREATER_EQUAL 3.13)
        target_link_options(coverage_config INTERFACE --coverage)
    else()
        target_link_libraries(coverage_config INTERFACE --coverage)
    endif()

    set(UTOPIA_ENABLE_CODE_COVERAGE ON)
endif()

macro(utopia_link_default_targets target)
    target_link_libraries(${target} PUBLIC coverage_config)
endmacro()

# 1) remove previous gcda files find . -name "*.gcda" -print0 | xargs -0 rm

# install lcov for post-processing the coverage see
# https://clang.llvm.org/docs/SourceBasedCodeCoverage.html#the-code-coverage-workflow
# LLVM_PROFILE_FILE=$PWD/utopia_cov.proraw ./utopia_test -test blas -test views
# -test expr -test qp_solver -test solvers -test utilities llvm-profdata merge
# -sparse utopia_cov.proraw -o utopia_cov.profdata llvm-cov show ./utopia_test
# -instr-profile=utopia_cov.profdata llvm-cov export -format=lcov -instr-profile
# utopia_cov.profdata ./utopia_test  > cov.info genhtml cov.info
# --output-directory out
