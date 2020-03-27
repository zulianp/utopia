#ifndef UTOPIA_PETSC_DM_OLD_HPP
#define UTOPIA_PETSC_DM_OLD_HPP

#include "utopia_petsc_Matrix.hpp"
#include "utopia_petsc_Vector.hpp"
#include "utopia_petsc_Communicator.hpp"
#include "utopia_petsc_ForwardDeclarations.hpp"
#include "utopia_Box.hpp"
#include "utopia_FunctionSpace.hpp"
#include "utopia_ArrayView.hpp"
#include "utopia_VectorView.hpp"
#include "utopia_Quadrature.hpp"
#include "utopia_UniformQuad4.hpp"
#include "utopia_UniformHex8.hpp"
#include "utopia_Temp.hpp"
#include "utopia_Writable.hpp"
#include "utopia_petsc_FE.hpp"
#include "utopia_DeviceView.hpp"
#include "utopia_SideSets.hpp"
#include "utopia_petsc_DMDA.hpp"

#include <array>
#include <memory>
#include <functional>

namespace utopia {
    template<int Dim>
    using PetscDM = utopia::PetscDMDA<StaticVector<PetscScalar, Dim>, ArrayView<PetscInt, Dim>>;

    // template<int Dim>
    // class DMDAMirror;

    // template<int Dim>
    // class PetscDMElements;

    // template<int Dim>
    // class DMDANodes;

    // template<class Space>
    // class DirichletBoundaryCondition;

    // // template<int Dim>
    // // class PetscNode {
    // // public:
    // //     PetscNode(const DMDANodes<Dim> &nodes, const SizeType &idx) : nodes_(nodes), idx_(idx) {}
    // //     SizeType idx() const { return idx_; }
    // //     bool is_ghost() const;

    // // private:
    // //     const DMDANodes<Dim> &nodes_;
    // //     SizeType idx_;
    // // };

    // template<int Dim>
    // class PetscDM final {
    // public:
    //     static const std::size_t UDim = Dim;
    //     using SizeType  = PetscInt;
    //     using Scalar    = PetscScalar;
    //     using Point     = utopia::StaticVector<Scalar, Dim>;
    //     // using Elem      = utopia::PetscElem<Dim>;
    //     // using Node      = utopia::PetscNode<Dim>;
    //     using Device    = utopia::Device<PETSC>;
    //     using NodeIndex = utopia::ArrayView<const SizeType>;
    //     using IntArray  = utopia::ArrayView<SizeType, UDim>;
    //     using ScalarArray = utopia::ArrayView<Scalar, UDim>;

    //     using SideSets = utopia::SideSets<Dim>;
    //     using Comm     = utopia::PetscCommunicator;

    //     class Impl;

    //     PetscDM(
    //         const PetscCommunicator     &comm,
    //         const std::array<SizeType, UDim> &dims,
    //         const std::array<Scalar, UDim>   &box_min,
    //         const std::array<Scalar, UDim>   &box_max,
    //         const SizeType &n_components = 1
    //     );

    //     void build(
    //         const PetscCommunicator     &comm,
    //         const std::array<SizeType, UDim> &dims,
    //         const std::array<Scalar, UDim>   &box_min,
    //         const std::array<Scalar, UDim>   &box_max,
    //         const SizeType &n_components = 1
    //     );

    //     void build_simplicial_complex(
    //         const PetscCommunicator     &comm,
    //         const std::array<SizeType, UDim> &dims,
    //         const std::array<Scalar, UDim>   &box_min,
    //         const std::array<Scalar, UDim>   &box_max,
    //         const SizeType &n_components = 1
    //     );

    //     PetscDM();
    //     ~PetscDM();

    //     PetscCommunicator &comm();
    //     const PetscCommunicator &comm() const;

    //     constexpr static typename SideSets::Sides sides()
    //     {
    //         return SideSets::sides();
    //     }

    //     void point(const SizeType &local_node_idx, Point &p) const;
    //     void cell_point(const SizeType &idx, Point &translation) const;
    //     void cell_size(const SizeType &idx, Point &cell_size) const;

    //     bool is_node_on_boundary(const SizeType &idx) const;
    //     bool is_node_on_boundary(const SizeType &idx, SideSet::BoundaryIdType b_id) const;
    //     void node(const SizeType &idx, Point &node) const;
    //     void elem(const SizeType &idx, Elem &e) const;
    //     void nodes(const SizeType &idx, NodeIndex &nodes) const;
    //     void nodes_local(const SizeType &idx, NodeIndex &nodes) const;

    //     void create_matrix(PetscMatrix &mat) const;
    //     void create_vector(PetscVector &vec) const;

    //     // void create_local_matrix(PetscMatrix &mat);
    //     void create_local_vector(PetscVector &vec) const;
    //     void local_to_global(const PetscVector &local,  PetscVector &global) const;
    //     void global_to_local(const PetscVector &global, PetscVector &local) const;

    //     void describe() const;

    //     inline static constexpr SizeType dim() { return Dim; }
    //     Range local_node_range() const;
    //     SizeType n_local_nodes_with_ghosts() const;
    //     Scalar min_spacing() const;

    //     Range element_range() const;

    //     const DMDAMirror<Dim> &mirror() const;

    //     // const IntArray &local_nodes_begin() const;
    //     // const IntArray &local_nodes_end() const;

    //     const IntArray &dims() const;
    //     const Point &box_min() const;
    //     const Point &box_max() const;

    //     inline Impl &impl()
    //     {
    //         return *impl_;
    //     }

    //     inline const Impl &impl() const
    //     {
    //         return *impl_;
    //     }

    //     SizeType n_nodes() const;
    //     SizeType n_local_nodes() const;

    //     SizeType n_elements() const;
    //     SizeType n_components() const;

    //     bool is_ghost(const SizeType &global_node_idx) const;
    //     bool is_boundary(const SizeType &global_node_idx) const;

    //     bool on_boundary(const SizeType &elem_idx) const;

    //     void set_field_name(const SizeType &nf, const std::string &name);
    //     void set_field_names(const std::vector<std::string> &names);

    //     std::unique_ptr<PetscDM> uniform_refine() const;

    //     void dmda_set_interpolation_type_Q0();
    //     void dmda_set_interpolation_type_Q1();
    //     void create_interpolation(const PetscDM &target, PetscMatrix &I) const;

    //     void update_mirror();

    //     std::unique_ptr<PetscDM> clone(const SizeType &n_components) const;
    //     std::unique_ptr<PetscDM> clone() const;
    // private:
    //     std::unique_ptr<Impl> impl_;
    // };
}

#endif //UTOPIA_PETSC_DM_OLD_HPP
