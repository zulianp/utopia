#ifndef UTOPIA_EACH_HPP
#define UTOPIA_EACH_HPP

#include "utopia_Base.hpp"
#include "utopia_For.hpp"
#include "utopia_ForwardDeclarations.hpp"
#include "utopia_Range.hpp"
#include "utopia_Readable.hpp"
#include "utopia_RowView.hpp"
#include "utopia_Size.hpp"
#include "utopia_Traits.hpp"

// #define UTOPIA_DISABLE_UNROLLING

namespace utopia {

    template <class Tensor, int Order = Tensor::Order, int FILL_TYPE = Traits<Tensor>::FILL_TYPE>
    class Each {};

    template <class Tensor, int FILL_TYPE>
    class Each<Tensor, 1, FILL_TYPE> {
    public:
        using SizeType = typename Traits<Tensor>::SizeType;
        using Scalar = typename Traits<Tensor>::Scalar;

        template <class Fun>
        inline static void apply_read(const Tensor &v, Fun fun) {
            const Range r = range(v);

#ifdef UTOPIA_DISABLE_UNROLLING
            for (auto i = r.begin(); i != r.end(); ++i) {
                fun(i, v.get(i));
            }
#else

            Read<Tensor> read_lock(v);
            For<>::apply(r.begin(), r.end(), [&v, &fun](const std::size_t i) { fun(i, v.get(i)); });
#endif  // UTOPIA_DISABLE_UNROLLING
        }

        template <class Fun>
        inline static void apply_write(Tensor &v, Fun fun) {
            Range r = range(v);

            Write<Tensor> write_lock(v);

#ifdef UTOPIA_DISABLE_UNROLLING

            for (auto i = r.begin(); i != r.end(); ++i) {
                v.set(i, fun(i));
            }

#else

            For<>::apply(r.begin(), r.end(), [&v, &fun](const std::size_t i) { v.set(i, fun(i)); });

#endif  // UTOPIA_DISABLE_UNROLLING
        }

        template <class Fun>
        inline static void apply_transform(Tensor &in_out, Fun fun) {
            apply_transform(in_out, in_out, fun);
        }

        template <class Fun>
        inline static void apply_transform(const Tensor &in, Tensor &out, Fun fun) {
            Range r = range(in);

            if (&in == &out) {
                ReadAndWrite<Tensor> rw(out);

                For<>::apply(r.begin(), r.end(), [&out, &fun](const std::size_t i) { out.set(i, fun(i, out.get(i))); });

                return;
            }

            if (size(in) != size(out)) {
                out.zeros(layout(in));
            }

            Read<Tensor> read_lock(in);
            Write<Tensor> write_lock(out);

#ifdef UTOPIA_DISABLE_UNROLLING
            for (auto i = r.begin(); i != r.end(); ++i) {
                out.set(i, fun(i, in.get(i)));
            }
#else
            For<>::apply(r.begin(), r.end(), [&in, &out, &fun](const std::size_t i) { out.set(i, fun(i, in.get(i))); });
#endif  // UTOPIA_DISABLE_UNROLLING
        }
    };

    template <class Tensor>
    class Each<Tensor, 2, FillType::DENSE> {
    public:
        using SizeType = typename Traits<Tensor>::SizeType;
        using Scalar = typename Traits<Tensor>::Scalar;

        template <class Fun>
        inline static void apply_read(const Tensor &m, Fun fun) {
            Range r = row_range(m);

            Size s = size(m);
            Read<Tensor> read_lock(m);

            for (auto i = r.begin(); i != r.end(); ++i) {
                for (auto j = 0; j != s.get(1); ++j) {
                    fun(i, j, m.get(i, j));
                }
            }
        }

        template <class Fun>
        inline static void apply_write(Tensor &m, Fun fun) {
            Range r = row_range(m);

            Size s = size(m);
            Write<Tensor> write_lock(m);

            for (auto i = r.begin(); i != r.end(); ++i) {
                for (auto j = 0; j != s.get(1); ++j) {
                    m.set(i, j, fun(i, j));
                }
            }
        }

        template <class Fun>
        inline static void apply_transform(Tensor &mat, Fun fun) {
            const SizeType cols = mat.cols();
            auto r = row_range(mat);

            ReadAndWrite<Tensor> rw(mat, utopia::LOCAL);
            for (SizeType i = r.begin(); i < r.end(); ++i) {
                for (SizeType j = 0; j < cols; ++j) {
                    mat.set(i, j, fun(i, j, mat.get(i, j)));
                }
            }
        }

        template <class Fun>
        inline static void apply_transform(const Tensor &in, Tensor &out, Fun fun) {
            assert(raw_type(in) != raw_type(out));

            Write<Tensor> w(out, utopia::LOCAL);

            auto r = row_range(in);
            for (auto i = r.begin(); i < r.end(); ++i) {
                RowView<const Tensor> row(in, i);
                for (size_t c = 0; c < row.n_values(); ++c) {
                    const SizeType j = row.col(c);
                    out.set(i, j, fun(i, j, row.get(c)));
                }
            }
        }
    };

    template <class Tensor>
    class Each<Tensor, 2, FillType::SPARSE> {
    public:
        using SizeType = typename Traits<Tensor>::SizeType;
        using Scalar = typename Traits<Tensor>::Scalar;

        template <class Fun>
        inline static void apply_read(const Tensor &m, Fun fun) {
            Range r = row_range(m);

            for (auto i = r.begin(); i != r.end(); ++i) {
                RowView<const Tensor> row(m, i);
                for (size_t c = 0; c < row.n_values(); ++c) {
                    const SizeType j = row.col(c);
                    fun(i, j, row.get(c));
                }
            }
        }

        template <class Fun>
        inline static void apply(Tensor &mat, Fun fun) {
            apply_transform(mat, [&fun](SizeType i, SizeType j, Scalar v) -> Scalar {
                UTOPIA_UNUSED(i);
                UTOPIA_UNUSED(j);
                return fun(v);
            });
        }

        template <class Fun>
        inline static void apply_transform(Tensor &mat, Fun fun) {
            Tensor mat_copy(mat);  // copy-constructor ensures copy of communicator
            apply_transform(mat_copy, mat, fun);
        }

        template <class Fun>
        inline static void apply_transform(const Tensor &in, Tensor &out, Fun fun) {
            assert(raw_type(in) != raw_type(out));

            out.scale(0.0);

            Write<Tensor> w(out, utopia::LOCAL);

            auto r = row_range(in);
            for (auto i = r.begin(); i < r.end(); ++i) {
                RowView<const Tensor> row(in, i);
                for (size_t c = 0; c < row.n_values(); ++c) {
                    const SizeType j = row.col(c);
                    out.set(i, j, fun(i, j, row.get(c)));
                }
            }
        }
    };

    template <class Tensor>
    class Each<Tensor, 2, FillType::POLYMORPHIC> : public Each<Tensor, 2, FillType::SPARSE> {
    public:
        template <class Fun>
        inline static void apply_write(Tensor &m, Fun fun) {
            if (m.is_dense()) {
                Each<Tensor, 2, FillType::DENSE>::apply_write(m, fun);
            } else {
                // Each<Tensor, 2, FillType::SPARSE>::apply_write(m, fun);
                assert(false && "IMPLEMENT ME");
                // m_utopia_error("NOT IMPLEMENTED FOR SPARSE MATRICES");
            }
        }
    };

    /** 	@defgroup element_acess Element Acess
     * 		@ingroup read_write
     *  	@brief  Actual acess to the elements of tensor
     */

    /**
     * @ingroup element_acess
     * @brief      Creates read lock on the tensor, iterates over all elements and applies provided function on them. \n
     * 			   Example usage: Printing vector v.
     *
         \code{.cpp}
         each_read(v, [](const SizeType i, const double entry) { utopia::out() <<"v(" << i << ") = " << entry <<
     std::endl;
     }); \endcode
     *
     *
     * @param[in]  v       The tensor.
     * @param[in]  fun     The  function with desirable action.
     */
    template <class T, int Order, class Fun>
    inline void each_read(const Tensor<T, Order> &v, Fun fun) {
        Each<T>::apply_read(v.derived(), fun);
    }

    /**
     * @ingroup element_acess
     * @brief      Creates write lock on the tensor, iterates over all elements and applies provided function on them.
     \n
     * 			   Example usage: Writing prescribed value to all elements of vector, but the first and the
     last. \code{.cpp} Vector v = zeros(10); const double value = 6.0;

            each_write(rhs, [value](const SizeType i) -> double {
                //The returned value will be written in the vector
                if(i == 0 || i == 10 - 1) {
                    return 0;}
                return value;
            });
        \endcode
     * @warning    If tensor is a sparse matrix, it will iterate only following the sparsity pattern.
     * @param[in]  v       The tensor.
     * @param[in]  fun     The  function with desirable action.
     */
    template <class T, int Order, class Fun>
    inline void each_write(Tensor<T, Order> &v, Fun fun) {
        Each<T>::apply_write(v.derived(), fun);
    }

    /**
     * @ingroup element_acess
     * @brief      Creates read lock on the tensor a and write lock on the tensor b, then applies provided function.
     *
     *
     * 			   Example usage: Applying a filter to the content of a and writing it into b.

         \code{.cpp}
             {
             const double n = 10;
                VectorT a = zeros(n);
                VectorT b = zeros(n);
                VectorT c = values(n, 0.5);

                //writing i/n in the vector
                each_write(a, [](const SizeType i) -> double  { return i/double(n); }   );

                {
                    //if another vector is needed, just provide a lock and pass it to the lambda functor
                    Read<VectorT> r(c);

                    //applying a filter to the content of a and writing it into b. a cannot be equal to b (for the
     moment) each_transform(a, b, [&c](const SizeType i, const double entry) -> double  { return exp(-entry) * c.get(i);
     }    );
                }
            }
        \endcode
     * @warning    Tensor a cannot be equal to the tensor b (for the moment).
     * @param[in]  a       The tensor to be red from/ input.
     * @param[in]  b       The tensor to be write into/ output.
     * @param[in]  fun     The  function with desirable action.
     */
    template <class T, int Order, class Fun>
    inline void each_transform(const Tensor<T, Order> &a, Tensor<T, Order> &b, Fun fun) {
        Each<T>::apply_transform(a.derived(), b.derived(), fun);
    }

    template <class T, int Order, class Fun>
    inline void each_transform(Tensor<T, Order> &t, Fun fun) {
        Each<T>::apply_transform(t.derived(), fun);
    }

    template <class T, int Order, class Fun>
    inline void each_apply(Tensor<T, Order> &t, Fun fun) {
        Each<T>::apply(t.derived(), fun);
    }
}  // namespace utopia

#endif  // UTOPIA_EACH_HPP
