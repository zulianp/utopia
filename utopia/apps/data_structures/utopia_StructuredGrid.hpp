#ifndef UTOPIA_STRUCTURED_GRID_HPP
#define UTOPIA_STRUCTURED_GRID_HPP

#include "utopia_Base.hpp"
#include "utopia_Traits.hpp"
#include "utopia_Algorithms.hpp"
#include "utopia_SideSets.hpp"

#include <type_traits>
#include <utility>

namespace utopia {

#define UTOPIA_CONSTEXPR_VOID constexpr void

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
        using SizeType = typename Traits<IntArray>::SizeType;
        using Scalar   = typename Traits<Point>::Scalar;

        UTOPIA_INLINE_FUNCTION constexpr SizeType dim() const
        {
            return dims_.size();
        }

        UTOPIA_INLINE_FUNCTION constexpr SizeType n_components() const
        {
            return n_components_;
        }

        UTOPIA_INLINE_FUNCTION constexpr SizeType n_nodes() const
        {
            return prod(dims_, dims_.size());
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

        // bool is_local_node_on_boundary(const SizeType &idx, SideSet::BoundaryIdType b_id) const
        // {
        //     IntArray tensor_index = {0, 0, 0};
        //     impl_->local_node_grid_coord_no_ghost(idx, tensor_index);

        //     const auto &mirror = impl_->mirror;

        //     switch(b_id) {
        //         case SideSet::left():
        //         {
        //             return tensor_index[0] == 0;
        //         }

        //         case SideSet::right():
        //         {
        //             return tensor_index[0] == (mirror.dims[0] - 1);
        //         }

        //         case SideSet::bottom():
        //         {
        //             return tensor_index[1] == 0;
        //         }

        //         case SideSet::top():
        //         {
        //             return tensor_index[1] == (mirror.dims[1] - 1);
        //         }

        //         case SideSet::back():
        //         {
        //             return tensor_index[2] == 0;
        //         }

        //         case SideSet::front():
        //         {
        //             return tensor_index[2] == (mirror.dims[2] - 1);
        //         }

        //         default:
        //         {
        //             break;
        //         }
        //     }

        //     return false;
        // }


        // constexpr static typename SideSets::Sides sides()
        // {
        //     return SideSets::sides();
        // }

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


        ///////////////////////////////// non-constexpr methods with gcc ////////////////////////

        template<class Array>
        UTOPIA_CONSTEXPR_VOID node_to_grid_coord_local_no_ghost(
            const SizeType &idx,
            Array &tensor_index
        ) const
        {
            const SizeType n = dim();

            SizeType current = idx;
            for(SizeType i = 0; i < n; ++i) {
                const SizeType next = current / corners_extent_[i];
                tensor_index[i] = current - next * corners_extent_[i] + corners_begin_[i];
                current = next;
            }
        }


        ///////////////////////////////// non-constexpr methods ////////////////////////
        UTOPIA_INLINE_FUNCTION StructuredGrid() {}

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
