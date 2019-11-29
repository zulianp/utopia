alias llvm-cov=llvm-cov-mp-7.0
alias llvm-profdata=llvm-profdata-mp-7.0

LLVM_PROFILE_FILE=$PWD/utopia_cov.proraw ./utopia_test
llvm-profdata merge -sparse utopia_cov.proraw -o utopia_cov.profdata
# llvm-cov show ./utopia_test -instr-profile=utopia_cov.profdata
llvm-cov export -format=lcov -instr-profile  utopia_cov.profdata ./utopia_test  > cov.info
# or
# lcov --no-external --capture --directory . --output-file cov.info
genhtml cov.info  --output-directory out