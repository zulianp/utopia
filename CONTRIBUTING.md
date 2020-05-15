### Pull Request Checklist

Before a PR can be merged, it should satisfy the following:

- [ ] Synchronize with target branch.
- [ ] Resolve all conflicts with target branch.
- [ ] Code builds using `make complete`.
- [ ] Code is formatted using clang-format (see [Code Formatting](https://bitbucket.org/zulianp/utopia/wiki/Code%20formatting))
- [ ] Run `make utest` to make sure all unit tests pass.
- [ ] Run `make bench` and make sure all benchmarks run through the end.
- [ ] Update `wiki`:
    - [ ] Is this a new feature users need to be aware of? New or updated example or app? If yes add a new short wiki entry for explaining usage.
- [ ] New capability:
   - [ ] All significant new classes, methods and functions have Doxygen-style documentation in source comments.
   - [ ] Consider adding new sample app in existing apps to highlight the new capability.
   - [ ] Consider saving cool simulation pictures with the new capability and add it to the `wiki` example page.
