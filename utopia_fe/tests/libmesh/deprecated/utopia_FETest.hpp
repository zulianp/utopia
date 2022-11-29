#ifndef UTOPIA_FE_TEST_HPP
#define UTOPIA_FE_TEST_HPP

#include "utopia_Testing.hpp"

#include "libmesh/parallel_mesh.h"

#include "utopia_Utils.hpp"

namespace utopia {

    class FETest {
    public:
        virtual ~FETest() {}

        inline void init(libMesh::Parallel::Communicator &comm) { comm_ = utopia::make_ref(comm); }
        virtual void run(Input &in) = 0;

        inline libMesh::Parallel::Communicator &comm() { return *comm_; }

    private:
        std::shared_ptr<libMesh::Parallel::Communicator> comm_;
    };
}  // namespace utopia

#endif  // UTOPIA_FE_TEST_HPP