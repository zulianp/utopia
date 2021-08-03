#ifndef UTOPIA_PETSC_IO_HPP
#define UTOPIA_PETSC_IO_HPP

#include <memory>

namespace utopia {
    // FIXME move this to file
    class Path;

    class PetscVector;
    class PetscCommunicator;

    namespace petsc {

        class DMBase;

        class PetscIO {
        public:
            class Wrapper;

            bool open(const PetscCommunicator &comm, const Path &path);

            bool write(const DMBase &dm);
            bool write(const PetscVector &x);

            void close();
            PetscIO();
            ~PetscIO();

        private:
            std::unique_ptr<Wrapper> wrapper_;
        };

    }  // namespace petsc
}  // namespace utopia

#endif  // UTOPIA_PETSC_IO_HPP
