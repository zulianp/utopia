#ifndef UTOPIA_STRUCTURED_GRID_HPP
#define UTOPIA_STRUCTURED_GRID_HPP

#include "utopia_Base.hpp"
#include "utopia_Traits.hpp"
#include "utopia_Algorithms.hpp"
#include "utopia_SideSets.hpp"
#include "utopia_CppMacros.hpp"

#include <type_traits>
#include <utility>

namespace utopia {

#define UTOPIA_CONSTEXPR_VOID void

    template<typename Array, typename SizeType>
    UTOPIA_INLINE_FUNCTION constexpr auto prod(const Array &vals, const SizeType &n) -> typename std::remove_reference<decltype(vals[0])>::type
    {
        if(n == 1) {
            return vals[0];
        } else {
            return vals[n-1] * prod(vals, n-1);
        }
    }

    template<typename Expr>
    UTOPIA_INLINE_FUNCTION constexpr auto prod(const Expr &expr) -> typename Traits<Expr>::Scalar
    {
        using SizeType = typename Traits<Expr>::SizeType;
        typename Traits<Expr>::Scalar ret = expr(0);

        const auto n = expr.size();
        for(SizeType i = 1; i < n; ++i) {
            ret *= expr(i);
        }

        return ret;
    }

    //follows petsc dm layout
    template<class Point, class IntArray, typename...>
    class StructuredGrid {
    public:
        using SizeType = typename Traits<IntArray>::ValueType;
        using Scalar   = typename Traits<Point>::Scalar;

        //if 0 it means it is dynamic and not static
        static constexpr const int StaticDim = Traits<IntArray>::StaticSize;

        ////////////////////////////////////////////////////////////////
        //////////////////////// GETTERS //////////////////////////////
        ////////////////////////////////////////////////////////////////

        UTOPIA_INLINE_FUNCTION constexpr SizeType dim() const
        {
            return dims_.size();
        }

        UTOPIA_INLINE_FUNCTION constexpr SizeType n_components() const
        {
            return n_components_;
        }

        UTOPIA_INLINE_FUNCTION constexpr const Point &box_min() const
        {
            return box_min_;
        }

        UTOPIA_INLINE_FUNCTION constexpr const Point &box_max() const
        {
            return box_max_;
        }

        UTOPIA_INLINE_FUNCTION constexpr const IntArray &dims() const { return dims_; }
        UTOPIA_INLINE_FUNCTION constexpr const IntArray &corners_begin() const { return corners_begin_; }
        UTOPIA_INLINE_FUNCTION constexpr const IntArray &corners_extent() const { return corners_extent_; }

        UTOPIA_INLINE_FUNCTION constexpr const IntArray &ghost_corners_begin() const { return ghost_corners_begin_; }
        UTOPIA_INLINE_FUNCTION constexpr const IntArray &ghost_corners_extent() const { return ghost_corners_extent_; }

        UTOPIA_INLINE_FUNCTION constexpr SizeType elements_x_cell() const { return elements_x_cell_; }
        UTOPIA_INLINE_FUNCTION constexpr SizeType dof_range_begin() const { return dof_range_begin_; }
        UTOPIA_INLINE_FUNCTION constexpr SizeType dof_range_end() const { return dof_range_end_; }

        ////////////////////////////////////////////////////////////////
        //////////////////////// SETTERS //////////////////////////////
        ////////////////////////////////////////////////////////////////

        UTOPIA_INLINE_FUNCTION Point &box_min()
        {
            return box_min_;
        }

        UTOPIA_INLINE_FUNCTION Point &box_max()
        {
            return box_max_;
        }

        UTOPIA_INLINE_FUNCTION IntArray &dims() { return dims_; }
        UTOPIA_INLINE_FUNCTION IntArray &corners_begin() { return corners_begin_; }
        UTOPIA_INLINE_FUNCTION IntArray &corners_extent() { return corners_extent_; }

        UTOPIA_INLINE_FUNCTION IntArray &ghost_corners_begin() { return ghost_corners_begin_; }
        UTOPIA_INLINE_FUNCTION IntArray &ghost_corners_extent() { return ghost_corners_extent_; }

        UTOPIA_INLINE_FUNCTION void set_n_components(const SizeType &val) { n_components_ = val; }
        UTOPIA_INLINE_FUNCTION void set_elements_x_cell(const SizeType &val) { elements_x_cell_ = val; }
        UTOPIA_INLINE_FUNCTION void set_dof_range_begin(const SizeType &val) { dof_range_begin_ = val; }
        UTOPIA_INLINE_FUNCTION void set_dof_range_end(const SizeType &val) { dof_range_end_ = val; }

        ////////////////////////////////////////////////////////////////

        UTOPIA_INLINE_FUNCTION constexpr SizeType n_nodes() const
        {
            return prod(dims_, dims_.size());
        }

        UTOPIA_INLINE_FUNCTION constexpr SizeType n_elements() const
        {
            const SizeType n = dim();

            SizeType ret = dims_[0] - 1;

            for(SizeType i = 1; i < n; ++i) {
                ret *= (dims_[i] - 1);
            }

            return ret * elements_x_cell_;
        }

        UTOPIA_INLINE_FUNCTION constexpr SizeType n_local_dofs() const
        {
            return dof_range_end_ - dof_range_begin_;
        }

        UTOPIA_INLINE_FUNCTION constexpr Scalar measure() const
        {
            return prod(box_max_ - box_min_);
        }

        UTOPIA_INLINE_FUNCTION constexpr Scalar min_spacing() const
        {
            const SizeType n = dim();

            Scalar min_h = (box_max_[0] - box_min_[0])/dims_[0];
            for(SizeType i = 1; i < n; ++i) {
                const Scalar h = (box_max_[i] - box_min_[i])/dims_[i];
                min_h = device::min(min_h, h);
            }

            return min_h;
        }

        UTOPIA_INLINE_FUNCTION constexpr Scalar cell_size(const SizeType &, const SizeType &d) const
        {
            return (box_max_[d] - box_min_[d])/(dims_[d] - 1);
        }

        UTOPIA_INLINE_FUNCTION constexpr bool is_node_on_boundary(const SizeType &idx) const //_local_no_ghost
        {
            const SizeType n = dim();

            SizeType current = idx;

            for(SizeType i = 0; i < n; ++i) {
                const SizeType next = current / corners_extent_[i];
                const SizeType tensor_index = current - next * corners_extent_[i] + corners_begin_[i];

                if(tensor_index == SizeType(0) || tensor_index == SizeType(dims_[i] - 1)) {
                    return true;
                }

                current = next;
            }

            return false;
        }

        UTOPIA_INLINE_FUNCTION void copy(const StructuredGrid &other)
        {
            n_components_    = other.n_components_;
            elements_x_cell_ = other.elements_x_cell_;
            dof_range_begin_ = other.dof_range_begin_;
            dof_range_end_   = other.dof_range_end_;

            device::copy(other.dims_, dims_);
            device::copy(other.corners_begin_, corners_begin_);
            device::copy(other.corners_extent_, corners_extent_);
            device::copy(other.ghost_corners_begin_, ghost_corners_begin_);
            device::copy(other.ghost_corners_extent_, ghost_corners_extent_);
            device::copy(other.box_min_, box_min_);
        }

        UTOPIA_INLINE_FUNCTION constexpr StructuredGrid(
            const IntArray &dims,
            const IntArray &corners_begin,
            const IntArray &corners_extent,
            const IntArray &ghost_corners_begin,
            const IntArray &ghost_corners_extent,
            const Point &box_min,
            const Point &box_max,
            const SizeType &n_components = 1,
            const SizeType &elements_x_cell = 1)
        : dims_(dims),
          corners_begin_(corners_begin),
          corners_extent_(corners_extent),
          ghost_corners_begin_(ghost_corners_begin),
          ghost_corners_extent_(ghost_corners_extent),
          n_components_(n_components),
          elements_x_cell_(elements_x_cell),
          dof_range_begin_(0),
          dof_range_end_(prod(dims, dims.size()) * n_components),
          box_min_(box_min),
          box_max_(box_max)
        {}

        ///////////////////////////////// void-constexpr methods (which do not work)////////////////////////

        template<class Array>
        UTOPIA_INLINE_FUNCTION UTOPIA_CONSTEXPR_VOID node_to_grid_coord(
            const SizeType &idx,
            Array &tensor_index
        ) const //_local_no_ghost
        {
            const SizeType n = dim();

            SizeType current = idx;
            for(SizeType i = 0; i < n; ++i) {
                const SizeType next = current / corners_extent_[i];
                tensor_index[i] = current - next * corners_extent_[i] + corners_begin_[i];
                current = next;
            }
        }

        UTOPIA_INLINE_FUNCTION UTOPIA_CONSTEXPR_VOID cell_size(const SizeType &elem_idx, Point &s) const
        {
            const SizeType n = dim();
            for(SizeType d = 0; d < n; ++d) {
                s[d] = cell_size(elem_idx, d);
            }
        }

        UTOPIA_INLINE_FUNCTION UTOPIA_CONSTEXPR_VOID point(const SizeType &local_node_idx, Point &p) const
        {
            const SizeType n = dim();

            SizeType current = local_node_idx;
            for(SizeType i = 0; i < n; ++i) {
                const SizeType next = current / ghost_corners_extent_[i];
                const SizeType tensor_index = current - next * ghost_corners_extent_[i];
                const SizeType c = tensor_index + ghost_corners_begin_[i];

                p[i] = c * (box_max_[i] - box_min_[i])/(dims_[i] - 1) + box_min_[i];

                //next
                current = next;
            }
        }

        ///////////////////////////////// non-constexpr methods ////////////////////////
        UTOPIA_INLINE_FUNCTION StructuredGrid() {}
        // UTOPIA_FUNCTION virtual ~StructuredGrid() {}

        UTOPIA_INLINE_FUNCTION bool is_node_on_boundary(const SizeType &idx, SideSet::BoundaryIdType b_id) const //_local_no_ghost
        {
            IntArray tensor_index;
            node_to_grid_coord(idx, tensor_index); //_local_no_ghost

            //FIXME use dim-dependent version
            return SideSets<StaticDim>::on_side(b_id, tensor_index, dims_);
        }

    private:
        IntArray dims_;
        IntArray corners_begin_;
        IntArray corners_extent_;

        IntArray ghost_corners_begin_;
        IntArray ghost_corners_extent_;

        // //elements
        SizeType n_components_;
        SizeType elements_x_cell_;

        // //dofs
        SizeType dof_range_begin_;
        SizeType dof_range_end_;

        // //geometry
        Point box_min_, box_max_;
    };

}

#endif //UTOPIA_STRUCTURED_GRID_HPP
