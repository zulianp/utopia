#ifndef UTOPIA_FE_FE_TESTS_HPP
#define UTOPIA_FE_FE_TESTS_HPP

#include "utopia_FETest.hpp"

#include "libmesh/parallel_mesh.h"
#include <string>
#include <ostream>


namespace utopia {

    class FETests {
    public:
        FETests();
        void run(libMesh::Parallel::Communicator &comm, int argc, char * argv[]);
        int add_test(const std::string &command, std::unique_ptr<FETest> &&app);
        void print_usage(std::ostream &os) const;

    private:
        std::map<std::string, std::unique_ptr<FETest>> tests_;


    };
}

#endif //UTOPIA_FE_FE_TESTS_HPP
