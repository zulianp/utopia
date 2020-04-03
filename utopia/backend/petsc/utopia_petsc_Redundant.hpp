#ifndef UTOPIA_PETSC_REDUNDANT_HPP
#define UTOPIA_PETSC_REDUNDANT_HPP

#include "utopia_petsc_Types.hpp"
#include "utopia_PetscIS.hpp"

namespace utopia {

    template<class M, class V>
    class Redundant {};


    class PetscVecScatter {
    public:
        void create(
            const PetscVector &from, const PetscIS &from_is,
            const PetscVector &to, const PetscIS &to_is
        );

        void apply(const PetscVector &from, PetscVector &to) const;
        void begin(const PetscVector &from, PetscVector &to) const;
        void end(const PetscVector &from, PetscVector &to) const;

        PetscVecScatter();
        ~PetscVecScatter();

    private:
        class Wrapper;
        std::shared_ptr<Wrapper> wrapper_;
    };

    template<>
    class Redundant<PetscMatrix, PetscVector> /*: public virtual Clonable*/ {
    public:
        using Layout = PetscTraits::Layout;

        Redundant();
        ~Redundant();

        // Redundant * clone() const override;
        void init(const Layout &lo, const SizeType n_sub_comm);

        void create_sub_vector(
            const PetscVector &vec,
            PetscVector &vec_sub,
            PetscVecScatter &scatter_to_sub,
            PetscVecScatter &scatter_to_super
        );

        void create_sub_matrix(const PetscMatrix &mat, PetscMatrix &mat_sub);
        void super_to_sub(const PetscMatrix &mat, PetscMatrix &mat_sub);

        // void super_to_sub(const PetscMatrix &mat, PetscMatrix &mat_sub);

        void super_to_sub(const PetscVector &vec,     const PetscVecScatter &scatter, PetscVector &vec_sub);
        void sub_to_super(const PetscVector &vec_sub, const PetscVecScatter &scatter, PetscVector &vec);

    private:
        PetscSubcomm   psubcomm;
        // PetscMatrix    pmats;  // check for update
        // PetscVector    sol_sub, rhs_sub, sol_dup, rhs_dup;
        // VecScatter     scatterin, scatterout;
        int n_sub_comm_;
        Layout sub_layout_;

        std::unique_ptr<PetscIS> is_super_to_sub_from, is_super_to_sub_to;
        std::unique_ptr<PetscIS> is_sub_to_super_from, is_sub_to_super_to;
    };

}
#endif //UTOPIA_PETSC_REDUNDANT_HPP

