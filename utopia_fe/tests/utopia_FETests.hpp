#ifndef UTOPIA_FE_FE_TESTS_HPP
#define UTOPIA_FE_FE_TESTS_HPP

#include "utopia_FETest.hpp"

#include <ostream>
#include <string>
#include "libmesh/parallel_mesh.h"

namespace utopia {

    class FETests {
    public:
        FETests();
        void run(libMesh::Parallel::Communicator &comm, int argc, char *argv[]);
        int add_test(const std::string &command, std::unique_ptr<FETest> &&app);
        void print_usage(std::ostream &os) const;

    private:
        std::map<std::string, std::unique_ptr<FETest>> tests_;
    };
}  // namespace utopia

#endif  // UTOPIA_FE_FE_TESTS_HPP
