#ifndef UTOPIA_PETSC_INTREPID2_HPP
#define UTOPIA_PETSC_INTREPID2_HPP

#include "utopia_Base.hpp"
#include "utopia_Field.hpp"

#include "utopia_intrepid2_Base.hpp"

#include "utopia_petsc_intrepid2_OmniAssembler.hpp"
// #include "utopia_petsc_intrepid2_Transport.hpp"

#ifdef UTOPIA_WITH_PETSC
#include "utopia_petsc_ForwardDeclarations.hpp"
#endif

#include "utopia_CreateFE.hpp"
#include "utopia_LocalToGlobal.hpp"

#include "utopia_intrepid2_ForwardDeclarations.hpp"
#include "utopia_kokkos_Field.hpp"
#include "utopia_petsc_ForwardDeclarations.hpp"
// #include "utopia_petsc_intrepid2_GradientField.hpp"

#include "utopia_petsc_FunctionSpace.hpp"

#include <Kokkos_DynRankView.hpp>

namespace utopia {

    using PetscScalar_t = utopia::Traits<utopia::petsc::FunctionSpace>::Scalar;

    template <typename T>
    using PetscViewDevice_t = utopia::intrepid2::ViewDevice<T>;
    using PetscIntViewDevice_t = utopia::intrepid2::ViewDevice<int>;

    template <typename Scalar>
    using Intrepid2Field = utopia::kokkos::Field<utopia::intrepid2::FE<Scalar>>;

    template <typename Scalar>
    class CreateFE<utopia::petsc::FunctionSpace, utopia::intrepid2::FE<Scalar>> {
    public:
        static void apply(const utopia::petsc::FunctionSpace &space,
                          utopia::intrepid2::FE<Scalar> &fe,
                          const int degree = 0);
    };

    template <typename Scalar>
    class CreateFEOnBoundary<utopia::petsc::FunctionSpace, utopia::intrepid2::FE<Scalar>> {
    public:
        static void apply(const utopia::petsc::FunctionSpace &space,
                          utopia::intrepid2::FE<Scalar> &fe,
                          const int degree = 2);

        static void apply(const utopia::petsc::FunctionSpace &space,
                          utopia::intrepid2::FE<Scalar> &fe,
                          const std::string &part_name,
                          const int degree = 2);
    };

    template <typename Scalar>
    class ConvertField<Field<utopia::petsc::FunctionSpace>, Intrepid2Field<Scalar>> {
    public:
        using Vector = utopia::Traits<utopia::petsc::FunctionSpace>::Vector;
        using SizeType = utopia::Traits<utopia::petsc::FunctionSpace>::SizeType;

        static void apply(const Field<utopia::petsc::FunctionSpace> &from, Intrepid2Field<Scalar> &to);
        static void apply(const Field<utopia::petsc::FunctionSpace> &from, Intrepid2Field<Scalar> &to, int var);
    };

#ifdef UTOPIA_WITH_PETSC
    template <typename Scalar>
    class LocalToGlobal<utopia::petsc::FunctionSpace, PetscViewDevice_t<Scalar>, PetscMatrix> {
    public:
        static void apply(const utopia::petsc::FunctionSpace &space,
                          const PetscViewDevice_t<Scalar> &element_matrices,
                          AssemblyMode mode,
                          PetscMatrix &matrix);
    };

    template <typename Scalar>
    class LocalToGlobal<utopia::petsc::FunctionSpace, PetscViewDevice_t<Scalar>, PetscVector> {
    public:
        static void apply(const utopia::petsc::FunctionSpace &space,
                          const PetscViewDevice_t<Scalar> &element_vectors,
                          AssemblyMode mode,
                          PetscVector &vector);

        static void apply(const utopia::petsc::FunctionSpace &space,
                          const PetscViewDevice_t<Scalar> &element_vectors,
                          AssemblyMode mode,
                          PetscVector &vector,
                          const int n_var);

        static void side_apply(const utopia::petsc::FunctionSpace &space,
                               const PetscViewDevice_t<Scalar> &element_vectors,
                               AssemblyMode mode,
                               PetscVector &vector,
                               const std::string &part_name);
    };

#endif  // UTOPIA_WITH_PETSC

    template <typename Scalar>
    class GlobalToLocal<utopia::petsc::FunctionSpace,
                        Traits<utopia::petsc::FunctionSpace>::Vector,
                        PetscViewDevice_t<Scalar>> {
    public:
        using Vector = utopia::Traits<utopia::petsc::FunctionSpace>::Vector;
        using SizeType = utopia::Traits<utopia::petsc::FunctionSpace>::SizeType;

        static void apply(const utopia::petsc::FunctionSpace &space,
                          const Vector &vector,
                          PetscViewDevice_t<Scalar> &element_vectors,
                          const int n_comp = 1);
    };

    void l2_norm(const Field<utopia::petsc::FunctionSpace> &field, std::vector<PetscScalar_t> &norms);

}  // namespace utopia

#endif  // UTOPIA_PETSC_INTREPID2_HPP