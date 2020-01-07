#ifndef UTOPIA_PETSC_DM_HPP
#define UTOPIA_PETSC_DM_HPP

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

#include <array>
#include <memory>
#include <functional>

namespace utopia {

    class SideSet {
    public:
        using BoundaryIdType = int;
        inline static constexpr BoundaryIdType left() { return 1; }
        inline static constexpr BoundaryIdType right() { return 2; }
        inline static constexpr BoundaryIdType bottom() { return 3; }
        inline static constexpr BoundaryIdType top() { return 4; }
        inline static constexpr BoundaryIdType front() { return 5; }
        inline static constexpr BoundaryIdType back() { return 6; }
        inline static constexpr BoundaryIdType invalid() { return 0; }
    };

    template<int Dim>
    class PetscDMElements;

    template<int Dim>
    class PetscDMNodes;

    template<class Space>
    class DirichletBoundaryCondition {};

    template<int Dim>
    class PetscNode {
    public:
        PetscNode(const PetscDMNodes<Dim> &nodes, const SizeType &idx) : nodes_(nodes), idx_(idx) {}
        SizeType idx() const { return idx_; }
        bool is_ghost() const;

    private:
        const PetscDMNodes<Dim> &nodes_;
        SizeType idx_;
    };

    template<int Dim>
    class PetscDM final {
    public:
        static const std::size_t UDim = Dim;
        using SizeType  = PetscInt;
        using Scalar    = PetscScalar;
        using Point     = utopia::StaticVector<Scalar, Dim>;
        using Elem      = utopia::PetscElem<Dim>;
        using Node      = utopia::PetscNode<Dim>;
        using Device    = utopia::Device<PETSC>;
        using NodeIndex = typename Elem::NodeIndex;

        class Impl;

        PetscDM(
            const PetscCommunicator     &comm,
            const std::array<SizeType, UDim> &dims,
            const std::array<Scalar, UDim>   &box_min,
            const std::array<Scalar, UDim>   &box_max
        );

        PetscDM();
        ~PetscDM();

        PetscCommunicator &comm();
        const PetscCommunicator &comm() const;

        void cell_point(const SizeType &idx, Point &translation);
        void cell_size(const SizeType &idx, Point &cell_size);

        bool is_local_node_on_boundary(const SizeType &idx, SideSet::BoundaryIdType b_id) const;
        bool is_node_on_boundary(const SizeType &idx, SideSet::BoundaryIdType b_id) const;
        void node(const SizeType &idx, Point &node) const;
        void elem(const SizeType &idx, Elem &e) const;
        void nodes(const SizeType &idx, NodeIndex &nodes) const;
        void nodes_local(const SizeType &idx, NodeIndex &nodes) const;

        void create_matrix(PetscMatrix &mat) const;
        void create_vector(PetscVector &vec) const;

        // void create_local_matrix(PetscMatrix &mat);
        void create_local_vector(PetscVector &vec) const;
        void local_to_global(const PetscVector &local,  PetscVector &global) const;
        void global_to_local(const PetscVector &global, PetscVector &local) const;

        void describe() const;

        void each_element(const std::function<void(const Elem &)> &f);
        void each_node(const std::function<void(const Node &)> &f);
        void each_node_with_ghosts(const std::function<void(const Node &)> &f);

        inline static constexpr SizeType dim() { return Dim; }
        Range local_node_range() const;
        SizeType n_local_nodes_with_ghosts() const;

        // void dims(SizeType *arr) const;
        void box(Scalar *min, Scalar *max) const;

        void local_node_ranges(SizeType *begin, SizeType *end) const;

        Range local_element_range() const;

#if UTOPIA_PETSC_VERSION_GREATER_EQUAL_THAN(3, 11, 0)
        void local_element_ranges(SizeType *begin, SizeType *end) const;

        template<class Array>
        void local_element_ranges(Array &begin, Array &end) const
        {
            SizeType temp_begin[3], temp_end[3];
            local_element_ranges(temp_begin, temp_end);

            SizeType d = dim();
            for(SizeType i = 0; i < d; ++i) {
                begin[i] = temp_begin[i];
                end[i] = temp_end[i];
            }
        }

#endif UTOPIA_PETSC_VERSION_GREATER_EQUAL_THAN(3, 11, 0)

        template<class Array>
        void local_node_ranges(Array &begin, Array &end) const
        {
            SizeType temp_begin[3], temp_end[3];
            local_node_ranges(temp_begin, temp_end);

            SizeType d = dim();
            for(SizeType i = 0; i < d; ++i) {
                begin[i] = temp_begin[i];
                end[i] = temp_end[i];
            }
        }

        template<class Array>
        void dims(Array &arr) const
        {
            std::array<SizeType, UDim> temp;
            dims(temp);
            SizeType d = dim();
            for(SizeType i = 0; i < d; ++i) {
                arr[i] = temp[i];
            }
        }

        void dims(std::array<SizeType, UDim> &arr) const;

        template<class Array>
        void box(Array &min, Array &max) const
        {
            Scalar temp_min[3], temp_max[3];
            box(temp_min, temp_max);

            SizeType d = dim();
            for(SizeType i = 0; i < d; ++i) {
                min[i] = temp_min[i];
                max[i] = temp_max[i];
            }
        }

        inline Impl &impl()
        {
            return *impl_;
        }

        inline const Impl &impl() const
        {
            return *impl_;
        }

        SizeType n_nodes() const;
        SizeType n_elements() const;

        bool is_ghost(const SizeType &global_node_idx) const;
        bool is_boundary(const SizeType &global_node_idx) const;
        SideSet::BoundaryIdType boundary_id(const SizeType &global_node_idx) const;

    private:
        std::unique_ptr<Impl> impl_;
    };
}

#endif //UTOPIA_PETSC_DM_HPP
