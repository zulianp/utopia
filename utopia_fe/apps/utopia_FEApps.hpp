#ifndef UTOPIA_FE_APPS_HPP
#define UTOPIA_FE_APPS_HPP

#include "utopia_FEApp.hpp"

#include <string>
#include "libmesh/parallel_mesh.h"

namespace utopia {

    class FEApps {
    public:
        FEApps();
        void run(libMesh::Parallel::Communicator &comm, int argc, char *argv[]);
        void add_app(const std::string &command, std::unique_ptr<FEApp> &&app);
        void print_usage(std::ostream &os) const;

    private:
        std::map<std::string, std::unique_ptr<FEApp>> apps_;
    };
}  // namespace utopia

#endif  // UTOPIA_FE_APPS_HPP
