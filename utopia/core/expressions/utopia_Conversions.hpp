#ifndef UTOPIA_UTOPIA_CONVERSIONS_HPP
#define UTOPIA_UTOPIA_CONVERSIONS_HPP

#include "utopia_Range.hpp"
#include "utopia_Factory.hpp"
#include "utopia_Base.hpp"
#include "utopia_Tensor.hpp"

namespace utopia {

    template<class TFrom, class TTo>
    void backend_convert_sparse(const Tensor<TFrom, 2> &t_from, Tensor<TTo, 2> &t_to)
    {
        using SizeType = typename Traits<TFrom>::SizeType;

        const auto &from = t_from.derived();
        auto &to = t_to.derived();

        auto ls = local_size(from);
        auto n_row_local = ls.get(0);
        std::vector<int> nnzxrow(n_row_local, 0);
        auto r = row_range(from);

        each_read(from, [&nnzxrow,&r](const SizeType i, const SizeType &, const double &) {
            ++nnzxrow[i - r.begin()];
        });


        auto nnz = *std::max_element(nnzxrow.begin(), nnzxrow.end());

        //FIXME use nnzxrow instead
        to.sparse(layout(from), nnz, nnz);

        {
            Write<TTo> w_t(to);
            each_read(from, [&to](const SizeType i, const SizeType j, const double val) {
                to.set(i, j, val);
            });
        }


        assert(size(from) == size(to));
        assert(local_size(from) == local_size(to));
    }

    template<class TFrom, class TTo>
    void backend_convert(const Tensor<TFrom, 1> &t_from, Tensor<TTo, 1> &t_to)
    {
        using SizeType = typename Traits<TFrom>::SizeType;

        const auto &from = t_from.derived();
        auto &to = t_to.derived();

        // auto ls = local_size(from).get(0);
        to.zeros(layout(from));

        Write<TTo> w_t(to);
        each_read(from, [&to](const SizeType i, const double val) {
            to.set(i, val);
        });

        assert(size(from) == size(to));
        assert(local_size(from) == local_size(to));
    }

//    template<class Vector, class Matrix>
//    void vec2mat(const Vector &v, Matrix &m, const bool transpose = false)
//    {
//        if(transpose) {
//            m = values(1, v.size().get(0), 0);
//            const Read<Vector> read(v);
//            const Write<Matrix> write(m);
//            const Range r = range(v);
//
//            for(int i = r.begin(); i < r.end(); ++i) {
//                m.set(0, i, v.get(i));
//            }
//
//        } else {
//            m = values(v.size().get(0), 1, 0);
//            const Read<Vector> read(v);
//            const Write<Matrix> write(m);
//            const Range r = range(v);
//
//            for(int i = r.begin(); i < r.end(); ++i) {
//                m.set(i, 0, v.get(i));
//            }
//        }
//    }

    template<class T1, class T2>
    bool cross_backend_approxeq(const Tensor<T1, 1> &l, const Tensor<T2, 1> &r)
    {
        T1 r_copy;
        backend_convert(r, r_copy);
        return approxeq(l, r_copy, 1e-10);
    }

    template<class T1, class T2>
    bool cross_backend_approxeq(const Tensor<T1, 2> &l, const Tensor<T2, 2> &r)
    {
        T1 r_copy;
        backend_convert_sparse(r, r_copy);
        return approxeq(l, r_copy, 1e-10);
    }
}

#endif //UTOPIA_UTOPIA_CONVERSIONS_HPP


