#ifndef UTOPIA_FE_APP_HPP
#define UTOPIA_FE_APP_HPP

#include "utopia_App.hpp"
#include "utopia_Utils.hpp"
#include "utopia_ui.hpp"

#include <string>
#include "libmesh/parallel_mesh.h"

namespace utopia {

    class FEApp : public App {
    public:
        virtual ~FEApp() {}

        inline void init(libMesh::Parallel::Communicator &comm) { comm_ = utopia::make_ref(comm); }

        inline libMesh::Parallel::Communicator &comm() { return *comm_; }

    private:
        std::shared_ptr<libMesh::Parallel::Communicator> comm_;
    };
}  // namespace utopia

#endif  // UTOPIA_FE_APP_HPP
