#ifndef UTOPIA_PETSC_LIBRARY_HPP
#define UTOPIA_PETSC_LIBRARY_HPP

#include "utopia_Library.hpp"

namespace utopia {
    class PetscLibrary : public Library {
    public:
        void init(int argc, char *argv[]) override;
        int finalize() override;

        std::string name() const override;
        std::string version_string() const override;
    };
}  // namespace utopia

#endif  // UTOPIA_PETSC_LIBRARY_HPP
