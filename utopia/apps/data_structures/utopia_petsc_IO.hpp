#ifndef UTOPIA_PETSC_IO_HPP
#define UTOPIA_PETSC_IO_HPP

#include <memory>

namespace utopia {

    //FIXME move this to file
    class Path;
    class PetscDMBase;
    class PetscVector;
    class PetscCommunicator;

    class PetscIO {
    public:

        class Wrapper;

        bool open(const PetscCommunicator &comm, const Path &path);

        bool write(const PetscDMBase &dm);
        bool write(const PetscVector &x);

        void close();
        PetscIO();
        ~PetscIO();

    private:
        std::unique_ptr<Wrapper> wrapper_;
    };

}

#endif //UTOPIA_PETSC_IO_HPP
