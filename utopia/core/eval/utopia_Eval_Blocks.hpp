#ifndef UTOPIA_EVAL_BLOCKS_HPP
#define UTOPIA_EVAL_BLOCKS_HPP

#include "utopia_Eval_Empty.hpp"
#include "utopia_Blocks.hpp"

namespace utopia {

    template<class Left, class Right, class Traits, int Backend>
    class Eval< Construct<Wrapper<Left, 1>, Blocks<Right> >, Traits, Backend> {
    public:
        using Tensor = utopia::Wrapper<Left, 1>;

        inline static bool apply(const Construct<Tensor, Blocks<Right> > &expr)
        {
            UTOPIA_TRACE_BEGIN(expr);
            auto &l = expr.left();
            const auto &b = expr.right().blocks();

            SizeType n = 0;

            for(auto b_ptr : b) {
                n += local_size(*b_ptr).get(0);
            }

            l = local_zeros(n);

            SizeType index = 0;

            {
                Write<Tensor> w_(l);
               
                for(auto b_ptr : b) {
                    auto rr = range(*b_ptr);
                    Read<Tensor> r_(*b_ptr);

                    for(auto i = rr.begin(); i < rr.end(); ++i) {
                        l.set(index++, b_ptr->get(i));
                    }
                }
            }

            UTOPIA_TRACE_END(expr);
            return true;
        }
    };

    template<class Left, class Right, class Traits, int Backend>
    class Eval< Construct<Wrapper<Left, 2>, Blocks<Right> >, Traits, Backend> {
    public:
        using Tensor = utopia::Wrapper<Left, 2>;
        using Scalar = UTOPIA_SCALAR(Tensor);

        inline static bool apply(const Construct<Tensor, Blocks<Right> > &expr)
        {
            UTOPIA_TRACE_BEGIN(expr);

            auto &l = expr.left();
            const auto &r = expr.right();

            SizeType rows = 0;
            SizeType cols = 0;

            std::vector<SizeType> row_offset(r.rows()+1, 0);
            std::vector<SizeType> col_offset(r.cols()+1, 0);

            for(SizeType i = 0; i < r.rows(); ++i) {
                for(SizeType j = 0; j < r.cols(); ++j) {
                    if(!r.block_is_null(i, j)) {
                        rows += local_size(r.block(i, j)).get(0);
                        row_offset[i+1] = rows;
                        break;
                    }
                }
            }

            for(SizeType j = 0; j < r.cols(); ++j) {
                for(SizeType i = 0; i < r.rows(); ++i) {
                    if(!r.block_is_null(i, j)) {
                        cols += local_size(r.block(i, j)).get(1);
                        col_offset[j+1] = cols;
                        break;
                    }
                }
            }

            SizeType max_nnz = 0;
            for(SizeType i = 0; i < r.rows(); ++i) {

                SizeType block_row_nnz = 0;

                for(SizeType j = 0; j < r.cols(); ++j) {
                    if(!r.block_is_null(i, j)) {
                        const auto &b = r.block(i, j);

                        std::vector<SizeType> nnz(local_size(b).get(0), 0);
                        auto rr = row_range(b);

                        each_read(b, [&](const SizeType r, const SizeType c, const Scalar val) {
                            nnz[r - rr.begin()];
                        });


                        block_row_nnz += *std::max_element(std::begin(nnz), std::end(nnz));
                    }
                }

                max_nnz = std::max(max_nnz, block_row_nnz);
            }

            l = local_sparse(rows, cols, max_nnz);

            {
                Write<Tensor> w_(l);

                for(SizeType i = 0; i < r.rows(); ++i) {
                    for(SizeType j = 0; j < r.cols(); ++j) {
                        if(!r.block_is_null(i, j)) {

                            each_read(r.block(i, j), [&](const SizeType r, const SizeType c, const Scalar val) {
                                l.set(row_offset[i] + r, col_offset[j] + c, val);
                            });
                        }
                    }
                }
            }

            UTOPIA_TRACE_END(expr);
            return true;
        }
    };

}

#endif //UTOPIA_EVAL_BLOCKS_HPP
