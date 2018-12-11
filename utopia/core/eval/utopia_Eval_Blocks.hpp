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
                assert((b_ptr));
                
                n += local_size(*b_ptr).get(0);
            }

            l = local_zeros(n);
            auto r = range(l);

            SizeType index = 0;

            {
                Write<Tensor> w_(l);
               
                for(auto b_ptr : b) {
                    auto rr = range(*b_ptr);
                    Read<Tensor> r_(*b_ptr);

                    for(auto i = rr.begin(); i < rr.end(); ++i) {
                        assert(index < n);
                        l.set(r.begin() + index++, b_ptr->get(i));
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

            utopia_test_assert(mpi_world_size() == 1 && "can only be used in serial");

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

            SizeType local_cols = 0;
            for(SizeType j = 0; j < r.cols(); ++j) {
                for(SizeType i = 0; i < r.rows(); ++i) {
                    if(!r.block_is_null(i, j)) {
                        local_cols += local_size(r.block(i, j)).get(1);
                        cols += size(r.block(i, j)).get(1);
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

                        each_read(b, [&](const SizeType r, const SizeType &, const Scalar &) {
                            ++nnz[r - rr.begin()];
                        });


                        if(!nnz.empty()) {
                            block_row_nnz += *std::max_element(std::begin(nnz), std::end(nnz));
                        }
                    }
                }

                max_nnz = std::max(max_nnz, block_row_nnz);
            }

            l = local_sparse(rows, local_cols, max_nnz);

            {
                Write<Tensor> w_(l);
                auto l_rr = row_range(l);

                for(SizeType i = 0; i < r.rows(); ++i) {
                    for(SizeType j = 0; j < r.cols(); ++j) {
                        if(!r.block_is_null(i, j)) {
                            const auto &b = r.block(i, j);
                            const auto b_rr = row_range(b);

                            const auto global_row_offset = l_rr.begin() - b_rr.begin() + row_offset[i];

                            each_read(b, [&](const SizeType r, const SizeType c, const Scalar val) {
                                l.set(
                                    global_row_offset + r,
                                    col_offset[j] + c, //BUG (the columns should be staggered to reflect the parallel decomposition)
                                    val
                                );
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
