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

namespace utopia {

    using LibMeshScalar_t = utopia::Traits<utopia::libmesh::FunctionSpace>::Scalar;

    template <typename T>
    using LibMeshViewDevice_t = utopia::intrepid2::ViewDevice<T>;
    using LibMeshIntViewDevice_t = utopia::intrepid2::ViewDevice<int>;

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
        using Vector = utopia::Traits<utopia::libmesh::FunctionSpace>::Vector;
        using SizeType = utopia::Traits<utopia::libmesh::FunctionSpace>::SizeType;

        static void apply(const Field<utopia::libmesh::FunctionSpace> &from, utopia::intrepid2::Field<Scalar> &to);
    };

    template <typename Scalar>
    class LocalToGlobal<utopia::libmesh::FunctionSpace, LibMeshViewDevice_t<Scalar>, PetscMatrix> {
    public:
        static void apply(const utopia::libmesh::FunctionSpace &space,
                          const LibMeshViewDevice_t<Scalar> &element_matrices,
                          AssemblyMode mode,
                          PetscMatrix &matrix);
    };

    template <typename Scalar>
    class LocalToGlobal<utopia::libmesh::FunctionSpace, LibMeshViewDevice_t<Scalar>, PetscVector> {
    public:
        static void apply(const utopia::libmesh::FunctionSpace &space,
                          const LibMeshViewDevice_t<Scalar> &element_vectors,
                          AssemblyMode mode,
                          PetscVector &vector);

        static void side_apply(const utopia::libmesh::FunctionSpace &space,
                               const LibMeshViewDevice_t<Scalar> &element_vectors,
                               AssemblyMode mode,
                               PetscVector &vector,
                               const std::string &part_name);
    };

    template <typename Scalar>
    class GlobalToLocal<utopia::libmesh::FunctionSpace,
                        Traits<utopia::libmesh::FunctionSpace>::Vector,
                        LibMeshViewDevice_t<Scalar>> {
    public:
        using Vector = utopia::Traits<utopia::libmesh::FunctionSpace>::Vector;
        using SizeType = utopia::Traits<utopia::libmesh::FunctionSpace>::SizeType;

        static void apply(const utopia::libmesh::FunctionSpace &space,
                          const Vector &vector,
                          LibMeshViewDevice_t<Scalar> &element_vectors,
                          const int n_comp = 1);
    };

}  // namespace utopia

#endif  // UTOPIA_LIBMESH_INTREPID2_HPP