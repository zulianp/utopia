#ifndef UTOPIA_STK_INTREPID2_HPP
#define UTOPIA_STK_INTREPID2_HPP

#include "utopia_Base.hpp"
#include "utopia_Field.hpp"

#include "utopia_intrepid2_Base.hpp"

#include "utopia_stk_intrepid2_OmniAssembler.hpp"
#include "utopia_stk_intrepid2_Transport.hpp"

#ifdef UTOPIA_ENABLE_PETSC
#include "utopia_petsc_ForwardDeclarations.hpp"
#endif

#include "utopia_CreateFE.hpp"
#include "utopia_LocalToGlobal.hpp"

#include "utopia_intrepid2_ForwardDeclarations.hpp"
#include "utopia_kokkos_Field.hpp"
#include "utopia_stk_ForwardDeclarations.hpp"
#include "utopia_stk_intrepid2_GradientField.hpp"

#include "utopia_stk_FunctionSpace.hpp"

#include <Kokkos_DynRankView.hpp>

namespace utopia {

    using StkScalar_t = utopia::Traits<utopia::stk::FunctionSpace>::Scalar;

    template <typename T>
    using StkViewDevice_t = utopia::intrepid2::ViewDevice<T>;
    using StkIntViewDevice_t = utopia::intrepid2::ViewDevice<int>;

    template <typename Scalar>
    using Intrepid2Field = utopia::kokkos::Field<utopia::intrepid2::FE<Scalar>>;

    template <typename Scalar>
    class CreateFE<utopia::stk::FunctionSpace, utopia::intrepid2::FE<Scalar>> {
    public:
        static void apply(const utopia::stk::FunctionSpace &space,
                          utopia::intrepid2::FE<Scalar> &fe,
                          const int degree = 0);
    };

    template <typename Scalar>
    class CreateFEOnBoundary<utopia::stk::FunctionSpace, utopia::intrepid2::FE<Scalar>> {
    public:
        static void apply(const utopia::stk::FunctionSpace &space,
                          utopia::intrepid2::FE<Scalar> &fe,
                          const int degree = 2);

        static void apply(const utopia::stk::FunctionSpace &space,
                          utopia::intrepid2::FE<Scalar> &fe,
                          const std::string &part_name,
                          const int degree = 2);
    };

    template <typename Scalar>
    class ConvertField<Field<utopia::stk::FunctionSpace>, Intrepid2Field<Scalar>> {
    public:
        using Vector = utopia::Traits<utopia::stk::FunctionSpace>::Vector;
        using SizeType = utopia::Traits<utopia::stk::FunctionSpace>::SizeType;

        static void apply(const Field<utopia::stk::FunctionSpace> &from, Intrepid2Field<Scalar> &to);
        static void apply(const Field<utopia::stk::FunctionSpace> &from, Intrepid2Field<Scalar> &to, int var);
    };

    template <typename Scalar>
    class ConvertField<Intrepid2Field<Scalar>, Field<utopia::stk::FunctionSpace>> {
    public:
        using Vector = utopia::Traits<utopia::stk::FunctionSpace>::Vector;
        using SizeType = utopia::Traits<utopia::stk::FunctionSpace>::SizeType;

        static void apply(const Intrepid2Field<Scalar> &from, Field<utopia::stk::FunctionSpace> &to);
    };

#ifdef UTOPIA_ENABLE_PETSC
    template <typename Scalar>
    class LocalToGlobal<utopia::stk::FunctionSpace, StkViewDevice_t<Scalar>, PetscMatrix> {
    public:
        static void apply(const utopia::stk::FunctionSpace &space,
                          const StkViewDevice_t<Scalar> &element_matrices,
                          AssemblyMode mode,
                          PetscMatrix &matrix);
    };

    template <typename Scalar>
    class LocalToGlobal<utopia::stk::FunctionSpace, StkViewDevice_t<Scalar>, PetscVector> {
    public:
        static void apply(const utopia::stk::FunctionSpace &space,
                          const StkViewDevice_t<Scalar> &element_vectors,
                          AssemblyMode mode,
                          PetscVector &vector);

        static void apply(const utopia::stk::FunctionSpace &space,
                          const StkViewDevice_t<Scalar> &element_vectors,
                          AssemblyMode mode,
                          PetscVector &vector,
                          const int n_var);

        static void side_apply(const utopia::stk::FunctionSpace &space,
                               const StkViewDevice_t<Scalar> &element_vectors,
                               AssemblyMode mode,
                               PetscVector &vector,
                               const std::string &part_name);
    };

#endif  // UTOPIA_ENABLE_PETSC

    template <typename Scalar>
    class GlobalToLocal<utopia::stk::FunctionSpace,
                        Traits<utopia::stk::FunctionSpace>::Vector,
                        StkViewDevice_t<Scalar>> {
    public:
        using Vector = utopia::Traits<utopia::stk::FunctionSpace>::Vector;
        using SizeType = utopia::Traits<utopia::stk::FunctionSpace>::SizeType;

        static void apply(const utopia::stk::FunctionSpace &space,
                          const Vector &vector,
                          StkViewDevice_t<Scalar> &element_vectors,
                          const int n_comp = 1);
    };

    void l2_norm(const Field<utopia::stk::FunctionSpace> &field, std::vector<StkScalar_t> &norms);

}  // namespace utopia

#endif  // UTOPIA_STK_INTREPID2_HPP