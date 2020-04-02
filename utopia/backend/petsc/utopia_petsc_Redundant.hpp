#ifndef UTOPIA_PETSC_REDUNDANT_HPP
#define UTOPIA_PETSC_REDUNDANT_HPP

#include "utopia_petsc_Types.hpp"

namespace utopia {

    template<class M, class V>
    class Redundant {};

    class PetscVecScatter;

    template<>
    class Redundant<PetscMatrix, PetscVector> : public virtual Clonable {
    public:
        using Layout = PetscTraits::Layout;

        Redundant();
        ~Redundant();

        Redundant * clone() const override;
        void init(const Layout &lo);

        void create_sub_vector(
            const PetscVector &vec,
            PetscVector &vec_sub
            PetscVecScatter &scatter_to_sub,
            PetscVecScatter &scatter_to_super
        );
        void create_sub_matrix(const PetscMatrix &mat, PetscMatrix &mat_sub);

        void super_to_sub(const PetscMatrix &mat, PetscMatrix &mat_sub);

        // void super_to_sub(const PetscVector &vec,     const PetscVecScatter &scatter, PetscVector &vec_sub);
        // void sub_to_super(const PetscVector &vec_sub, const PetscVecScatter &scatter, PetscVector &vec);

    private:
        PetscSubcomm   psubcomm;
        PetscMatrix    pmats;  // check for update
        PetscVector    sol_sub, rhs_sub, sol_dup, rhs_dup;
        VecScatter     scatterin, scatterout;
        int n_sub_comm_;
    };

}
#endif //UTOPIA_PETSC_REDUNDANT_HPP

