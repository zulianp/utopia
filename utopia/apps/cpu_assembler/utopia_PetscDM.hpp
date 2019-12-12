#ifndef UTOPIA_PETSC_DM_HPP
#define UTOPIA_PETSC_DM_HPP

#include <memory>
#include "utopia_petsc_Vector.hpp"
#include "utopia_petsc_Communicator.hpp"
#include "utopia_petsc_ForwardDeclarations.hpp"

namespace utopia {

    class PetscDM final {
    public:
        using SizeType = PetscInt;
        using Scalar   = PetscScalar;

        class Impl;
        class Elements;

        void create_matrix(PetscMatrix &mat);

        PetscDM(
            const PetscCommunicator     &comm,
            const std::vector<SizeType> &dims,
            const std::vector<Scalar>   &box_min,
            const std::vector<Scalar>   &box_max
        );

        PetscDM();
        ~PetscDM();

        void describe() const;

        class Elem {
        public:

            Elem(const Elements &e, const SizeType &idx) : e_(e), idx_(idx) {}
            SizeType n_nodes() const;
            SizeType node_id(const SizeType k) const;

        private:
            const Elements &e_;
            SizeType idx_;
        };

        void each_element(const std::function<void(const SizeType &, const Elem &e)> &f);

    private:
        std::unique_ptr<Impl> impl_;
    };

}

#endif //UTOPIA_PETSC_DM_HPP
