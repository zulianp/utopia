#ifndef UTOPIA_LIBMESH_INTREPID2_HPP
#define UTOPIA_LIBMESH_INTREPID2_HPP

#include "utopia_Base.hpp"
#include "utopia_Traits.hpp"

#include "utopia_Field.hpp"

#include "utopia_CreateFE.hpp"
#include "utopia_LocalToGlobal.hpp"
#include "utopia_petsc_ForwardDeclarations.hpp"

#include "utopia_intrepid2_Field.hpp"
#include "utopia_intrepid2_ForwardDeclarations.hpp"
#include "utopia_libmesh_ForwardDeclarations.hpp"

#include "utopia_libmesh_FunctionSpace_new.hpp"

#include <Kokkos_DynRankView.hpp>

namespace utopia {

    template <typename Scalar>
    class CreateFE<utopia::libmesh::FunctionSpace, utopia::intrepid2::FE<Scalar>> {
    public:
        static void apply(const utopia::libmesh::FunctionSpace &space,
                          utopia::intrepid2::FE<Scalar> &fe,
                          const int degree = 0);
    };

    template <typename Scalar>
    class CreateFEOnBoundary<utopia::libmesh::FunctionSpace, utopia::intrepid2::FE<Scalar>> {
    public:
        static void apply(const utopia::libmesh::FunctionSpace &space,
                          utopia::intrepid2::FE<Scalar> &fe,
                          const int degree = 2);

        static void apply(const utopia::libmesh::FunctionSpace &space,
                          utopia::intrepid2::FE<Scalar> &fe,
                          const std::string &part_name,
                          const int degree = 2);
    };

    template <typename Scalar>
    class ConvertField<Field<utopia::libmesh::FunctionSpace>, utopia::intrepid2::Field<Scalar>> {
    public:
        using DynRankView = ::Kokkos::DynRankView<Scalar>;
        using Vector = utopia::Traits<utopia::libmesh::FunctionSpace>::Vector;
        using SizeType = utopia::Traits<utopia::libmesh::FunctionSpace>::SizeType;

        static void apply(const Field<utopia::libmesh::FunctionSpace> &from, utopia::intrepid2::Field<Scalar> &to);
    };

    template <typename Scalar>
    class LocalToGlobal<utopia::libmesh::FunctionSpace, ::Kokkos::DynRankView<Scalar>, PetscMatrix> {
    public:
        using DynRankView = ::Kokkos::DynRankView<Scalar>;
        static void apply(const utopia::libmesh::FunctionSpace &space,
                          const DynRankView &element_matrices,
                          AssemblyMode mode,
                          PetscMatrix &matrix);
    };

    template <typename Scalar>
    class LocalToGlobal<utopia::libmesh::FunctionSpace, ::Kokkos::DynRankView<Scalar>, PetscVector> {
    public:
        using DynRankView = ::Kokkos::DynRankView<Scalar>;
        static void apply(const utopia::libmesh::FunctionSpace &space,
                          const DynRankView &element_vectors,
                          AssemblyMode mode,
                          PetscVector &vector);

        static void side_apply(const utopia::libmesh::FunctionSpace &space,
                               const DynRankView &element_vectors,
                               AssemblyMode mode,
                               PetscVector &vector,
                               const std::string &part_name);
    };

    template <class FunctionSpace, class ElementTensors, class GlobalTensor, typename... Args>
    void side_local_to_global(const FunctionSpace &space,
                              const ElementTensors &element_matrices,
                              AssemblyMode mode,
                              GlobalTensor &tensor,
                              Args &&... args) {
        LocalToGlobal<FunctionSpace, ElementTensors, GlobalTensor>::side_apply(
            space, element_matrices, mode, tensor, args...);
    }

    template <typename Scalar>
    class GlobalToLocal<utopia::libmesh::FunctionSpace,
                        Traits<utopia::libmesh::FunctionSpace>::Vector,
                        ::Kokkos::DynRankView<Scalar>> {
    public:
        using DynRankView = ::Kokkos::DynRankView<Scalar>;
        using Vector = utopia::Traits<utopia::libmesh::FunctionSpace>::Vector;
        using SizeType = utopia::Traits<utopia::libmesh::FunctionSpace>::SizeType;

        static void apply(const utopia::libmesh::FunctionSpace &space,
                          const Vector &vector,
                          DynRankView &element_vectors,
                          const int n_comp = 1);
    };

}  // namespace utopia

#endif  // UTOPIA_LIBMESH_INTREPID2_HPP