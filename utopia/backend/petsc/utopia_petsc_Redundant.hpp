#ifndef UTOPIA_PETSC_REDUNDANT_HPP
#define UTOPIA_PETSC_REDUNDANT_HPP

#include "utopia_PetscIS.hpp"
#include "utopia_petsc_Types.hpp"

namespace utopia {

    template <class M, class V>
    class Redundant {};

    class PetscVecScatter {
    public:
        void create(const PetscVector &from, const PetscIS &from_is, const PetscVector &to, const PetscIS &to_is);

        void apply(const PetscVector &from, PetscVector &to) const;
        void begin(const PetscVector &from, PetscVector &to) const;
        void end(const PetscVector &from, PetscVector &to) const;

        PetscVecScatter();
        ~PetscVecScatter();

    private:
        class Wrapper;
        std::shared_ptr<Wrapper> wrapper_;
    };

    template <>
    class Redundant<PetscMatrix, PetscVector> /*: public virtual Clonable*/ {
    public:
        using Layout = PetscTraits::Layout;

        Redundant();
        ~Redundant();

        // Redundant * clone() const override;
        void init(const Layout &lo, const SizeType n_sub_comm);

        inline void init(const Layout &lo) { init(lo, n_sub_comm_); }

        inline void n_sub_comm(const int n_sub_comm) { n_sub_comm_ = n_sub_comm; }

        inline int n_sub_comm() const { return n_sub_comm_; }

        void create_sub_vector(PetscVector &vec_sub);

        void create_sub_matrix(const PetscMatrix &mat, PetscMatrix &mat_sub);
        void super_to_sub(const PetscMatrix &mat, PetscMatrix &mat_sub);

        // void super_to_sub(const PetscMatrix &mat, PetscMatrix &mat_sub);

        void super_to_sub(const PetscVector &vec, PetscVector &vec_sub);
        void sub_to_super(const PetscVector &vec_sub, PetscVector &vec);

        bool empty() const { return !static_cast<bool>(psubcomm); }

    private:
        PetscSubcomm psubcomm{nullptr};
        PetscVector buff_, empty_;

        int n_sub_comm_{2};
        Layout sub_layout_, child_layout_;

        std::unique_ptr<PetscIS> is_super_to_sub_from, is_super_to_sub_to;
        std::unique_ptr<PetscIS> is_sub_to_super_from, is_sub_to_super_to;
        PetscVecScatter scatter_to_sub;
        PetscVecScatter scatter_to_super;

        inline Redundant(const Redundant &) { assert(false); }

        inline Redundant &operator=(const Redundant &) {
            assert(false);
            return *this;
        }
    };
}  // namespace utopia

#endif  // UTOPIA_PETSC_REDUNDANT_HPP
