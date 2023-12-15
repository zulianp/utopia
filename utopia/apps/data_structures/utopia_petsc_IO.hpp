#ifndef UTOPIA_PETSC_IO_HPP
#define UTOPIA_PETSC_IO_HPP

#include <memory>
#include <string>

namespace utopia {

    // FIXME move this to file
    class Path;
    class PetscDMBase;
    class PetscVector;
    class PetscCommunicator;

    class PetscIO {
    public:
        class Wrapper;

        bool open(const PetscCommunicator &comm, const Path &path, const std::string &mode = "w");

        bool write(const PetscDMBase &dm);
        bool write(const PetscVector &x);
        bool read(PetscDMBase &dm);
        bool read(PetscVector &x);

        void close();
        PetscIO();
        ~PetscIO();

    private:
        std::unique_ptr<Wrapper> wrapper_;
    };

}  // namespace utopia

#endif  // UTOPIA_PETSC_IO_HPP
