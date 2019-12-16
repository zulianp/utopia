#ifndef UTOPIA_PETSC_DM_HPP
#define UTOPIA_PETSC_DM_HPP


#include "utopia_petsc_Vector.hpp"
#include "utopia_petsc_Communicator.hpp"
#include "utopia_petsc_ForwardDeclarations.hpp"
#include "utopia_Box.hpp"
#include "utopia_FunctionSpace.hpp"
#include "utopia_ArrayView.hpp"
#include "utopia_VectorView.hpp"

#include <array>
#include <memory>

namespace utopia {

    template<int Dim>
    class PetscDMElements;

    template<int Dim>
    class PetscDMNodes;

    template<int Dim>
    class PetscElem {
    public:
        using SizeType = PetscInt;
        using Scalar   = PetscScalar;
        using Point    = utopia::StaticVector<Scalar, Dim>;
        using Grad     = utopia::StaticVector<Scalar, Dim>;

        virtual ~PetscElem() {}
        PetscElem()
        {}

        inline void set(const ArrayView<SizeType> &nodes)
        {
            nodes_ = nodes;
        }

        inline const ArrayView<SizeType> &nodes() const
        {
            return nodes_;
        }

        inline SizeType node_id(const SizeType k) const
        {
            return nodes_[k];
        }

        inline const SizeType &node(const SizeType &i) const
        {
            return nodes_[i];
        }

        inline SizeType n_nodes() const
        {
            return nodes_.size();
        }

        inline constexpr static SizeType dim() { return Dim; }
        virtual Scalar fun(const SizeType &i, const Point &p) const = 0;

        inline void idx(const SizeType &idx)
        {
            idx_ = idx;
        }

        inline SizeType idx() const
        {
            return idx_;
        }

    private:
        ArrayView<SizeType> nodes_;
        SizeType idx_;
    };

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
        using Point     = utopia::VectorView<ArrayView<Scalar>>;
        using Elem = utopia::PetscElem<Dim>;
        using Node = utopia::PetscNode<Dim>;

        class Impl;

        PetscDM(
            const PetscCommunicator     &comm,
            const std::array<SizeType, UDim> &dims,
            const std::array<Scalar, UDim>   &box_min,
            const std::array<Scalar, UDim>   &box_max
        );


        PetscDM();
        ~PetscDM();

        void create_matrix(PetscMatrix &mat);
        void create_vector(PetscVector &vec);

        void describe() const;

        void each_element(const std::function<void(const Elem &)> &f);
        void each_node(const std::function<void(const Node &)> &f);
        void each_node_with_ghosts(const std::function<void(const Node &)> &f);

        inline static constexpr SizeType dim() { return Dim; }
        Range local_node_range() const;
        SizeType n_local_nodes_with_ghosts() const;

        void dims(SizeType *arr) const;
        void box(Scalar *min, Scalar *max) const;

        void local_node_ranges(SizeType *begin, SizeType *end) const;
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

    private:
        std::unique_ptr<Impl> impl_;
    };

    template<int Dim>
    class FunctionSpace<PetscDM<Dim>, 1> {
    public:
        using Mesh = utopia::PetscDM<Dim>;

        FunctionSpace(const std::shared_ptr<Mesh> &mesh)
        : mesh_(mesh)
        {}

        template<class Fun>
        void each_element(Fun fun)
        {
            mesh_->each_element(fun);
        }

        FunctionSpace &view_device()
        {
            return *this;
        }

        // Range local_element_range() const
        // {
        //     return mesh_.local_element_range();
        // }

    private:
        std::shared_ptr<Mesh> mesh_;

        FunctionSpace(const FunctionSpace &other)
        : mesh_(other.mesh_)
        {}
    };

}

#endif //UTOPIA_PETSC_DM_HPP
