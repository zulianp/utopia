#ifndef UTOPIA_LIBMESH_LIBRARY_HPP
#define UTOPIA_LIBMESH_LIBRARY_HPP

#include "utopia_Library.hpp"

#include <memory>
#include <string>

namespace utopia {

    class LibMeshLibrary final : public Library {
    public:
        LibMeshLibrary();
        ~LibMeshLibrary();
        void init(int argc, char *argv[]);
        int finalize();
        std::string name() const;
        std::string version_string() const;

    private:
        class Impl;
        std::unique_ptr<Impl> impl_;
    };

}  // namespace utopia

#endif  // UTOPIA_LIBMESH_LIBRARY_HPP
