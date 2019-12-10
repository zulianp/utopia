#ifndef UTOPIA_ASSEMBLY_VIEW_HPP
#define UTOPIA_ASSEMBLY_VIEW_HPP

#include "utopia_FunctionSpaceView.hpp"

namespace utopia {

    template<class Elem, class Quadrature, class MemType = typename Elem::MemType, typename...>
    class PhysicalPoint {
    public:
        static const int Dim = Elem::Dim;

        UTOPIA_INLINE_FUNCTION PhysicalPoint(const Elem *elem = nullptr, const Quadrature *q = nullptr)
        : elem_(elem), q_(q)
        {}

        template<class Point>
        UTOPIA_INLINE_FUNCTION void get(const int qp_idx, Point &p) const
        {
            Point temp;
            q_->point(qp_idx, temp);
            elem_->point(temp, p);
        }

        UTOPIA_INLINE_FUNCTION std::size_t size() const
        {
            return q_->n_points();
        }

        UTOPIA_INLINE_FUNCTION void update(const Elem *elem, const Quadrature *q)
        {
            elem_ = elem;
            q_    = q;
        }

    private:
        const Elem *elem_;
        const Quadrature *q_;
    };

    template<class Mesh, int NComponents>
    class PhysicalPoint< FunctionSpace<Mesh, NComponents> > {
    public:

    };

    template<class Elem, class Quadrature, class MemType = typename Elem::MemType, typename...>
    class PhysicalGradient {
    public:
        static const int Dim = Elem::Dim;
        using Scalar = typename Elem::Scalar;
        using Point  = utopia::StaticVector<Scalar, Dim>;
        using Grad   = utopia::StaticVector<Scalar, Dim>;

        UTOPIA_INLINE_FUNCTION PhysicalGradient(const Elem *elem = nullptr, const Quadrature *q = nullptr)
        : elem_(elem), q_(q)
        {}

        template<class Grad>
        UTOPIA_INLINE_FUNCTION void get(const int fun_num, const int qp_idx, Grad &g) const
        {
            Point temp;
            q_->point(qp_idx, temp);
            elem_->grad(fun_num, temp, g);
        }

        UTOPIA_INLINE_FUNCTION Grad operator()(const int fun_num, const int qp_idx)
        {
            Grad g;
            get(fun_num, qp_idx, g);
            return g;
        }

        UTOPIA_INLINE_FUNCTION std::size_t size() const
        {
            return q_->n_points();
        }

        UTOPIA_INLINE_FUNCTION void update(const Elem *elem, const Quadrature *q)
        {
            elem_ = elem;
            q_    = q;
        }

    private:
        const Elem *elem_;
        const Quadrature *q_;
    };

    template<class Elem, class Quadrature, class MemType = typename Elem::MemType, typename...>
    class ShapeFunction {
    public:
        static const int Dim = Elem::Dim;
        using Scalar = typename Elem::Scalar;
        using Point  = utopia::StaticVector<Scalar, Dim>;

        UTOPIA_INLINE_FUNCTION ShapeFunction(const Elem *elem = nullptr, const Quadrature *q = nullptr)
        : elem_(elem), q_(q)
        {}

        UTOPIA_INLINE_FUNCTION Scalar operator()(const int fun_num, const int qp_idx) const
        {
            return get(fun_num, qp_idx);
        }

        UTOPIA_INLINE_FUNCTION Scalar get(const int fun_num, const int qp_idx) const
        {
            Point temp;
            q_->point(qp_idx, temp);
            return elem_->fun(fun_num, temp);
        }

        UTOPIA_INLINE_FUNCTION std::size_t size() const
        {
            return q_->n_points();
        }

        UTOPIA_INLINE_FUNCTION void update(const Elem *elem, const Quadrature *q)
        {
            elem_ = elem;
            q_    = q;
        }

    private:
        const Elem *elem_;
        const Quadrature *q_;
    };

    template<class Elem, class Quadrature, class MemType = typename Elem::MemType, typename...>
    class Differential {
    public:
        static const int Dim = Elem::Dim;
        using Scalar = typename Elem::Scalar;
        using Point  = utopia::StaticVector<Scalar, Dim>;

        UTOPIA_INLINE_FUNCTION Differential(const Elem *elem = nullptr, const Quadrature *q = nullptr)
        : elem_(elem), q_(q)
        {}

        UTOPIA_INLINE_FUNCTION Scalar operator()(const int qp_idx) const
        {
            return get(qp_idx);
        }

        UTOPIA_INLINE_FUNCTION Scalar get(const int qp_idx) const
        {
            if(elem_->is_affine()) {
                return q_->weight(qp_idx) * elem_->reference_measure();
            } else {
                UTOPIA_DEVICE_ASSERT(false);
                return -1.0;
            }
        }

        UTOPIA_INLINE_FUNCTION std::size_t size() const
        {
            return q_->n_points();
        }

        UTOPIA_INLINE_FUNCTION void update(const Elem *elem, const Quadrature *q)
        {
            elem_ = elem;
            q_    = q;
        }

    private:
        const Elem *elem_;
        const Quadrature *q_;
    };

    template<class FunctionSpace, class Vector, class Elem, class MemType = typename Elem::MemType, typename...>
    class NodalInterpolate {
    public:
        using DofIndex = typename FunctionSpace::DofIndex;
        using SizeType = typename utopia::Traits<Vector>::SizeType;

        //convert nodal condensed to element-wise map

        NodalInterpolate(const FunctionSpace &space)
        : space_(space), element_wise_values_(local_zeros(space.n_local_elements() * Elem::NNodes))
        {}

        NodalInterpolate(const FunctionSpace &space, const Vector &values)
        : space_(space), element_wise_values_(local_zeros(space.n_local_elements() * Elem::NNodes))
        {
            update(values);
        }

        //values must have ghost values
        void update(const Vector &values)
        {
            auto other_view = const_device_view(values);
            auto this_view  = device_view(element_wise_values_);

            space_.each_element(UTOPIA_LAMBDA(const SizeType &i, const Elem *elem) {
                DofIndex index;
                space_.dofs(i, index);

                const auto n = index.size();
                const auto offset = i * Elem::NNodes;

                for(SizeType k = 0; k < n; ++k) {
                    this_view.set(offset + k, other_view.get(index[k]));
                }
            });
        }

        // ElementNodalInterpolate get(const SizeType &idx) const
        // {

        // }

    private:
        const FunctionSpace &space_;
        Vector element_wise_values_;
    };


}

#endif //UTOPIA_ASSEMBLY_VIEW_HPP
