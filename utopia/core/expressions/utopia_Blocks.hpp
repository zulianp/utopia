#ifndef UTOPIA_BLOCKS_HPP
#define UTOPIA_BLOCKS_HPP

#include "utopia_ForwardDeclarations.hpp"
#include "utopia_Size.hpp"

#include <vector>
#include <memory>
#include <initializer_list>

namespace utopia {
    template<class Tensor, int Order = Tensor::Order>
    class Blocks {};

    template<class Matrix>
    class Blocks<Matrix, 2> : public Expression<Blocks<Matrix, 2>> {
    public:
        using MatrixPtrT = std::shared_ptr<const Matrix>;
        using SizeType = typename utopia::Traits<Matrix>::SizeType;

        Blocks(const SizeType rows, const SizeType cols, const std::vector<MatrixPtrT> &blocks)
        : rows_(rows), cols_(cols), blocks_(blocks)
        {
            assert(rows_*cols_ == SizeType(blocks_.size()));
        }

        Blocks(const SizeType rows, const SizeType cols)
        : rows_(rows), cols_(cols)
        {
            blocks_.resize(rows_*cols_);
        }

        const std::vector<MatrixPtrT> &blocks() const
        {
            return blocks_;
        }

        inline SizeType rows() const
        {
            return rows_;
        }

        inline SizeType cols() const
        {
            return cols_;
        }

        const Matrix &block(const SizeType i, const SizeType j) const
        {
            assert(i < rows());
            assert(j < cols());
            assert((blocks_[i*cols_ + j]));

            return *blocks_[i*cols_ + j];
        }

        const MatrixPtrT &block_ptr(const SizeType i, const SizeType j) const
        {
            assert(i < rows());
            assert(j < cols());
            return blocks_[i*cols_ + j];
        }

        bool block_is_null(const SizeType i, const SizeType j) const
        {
            assert(i < rows());
            assert(j < cols());
            return !bool(blocks_[i*cols_ + j]);
        }

        void set_block(const SizeType i, const SizeType j, const MatrixPtrT &b)
        {
            assert(i < rows());
            assert(j < cols());
            blocks_[i*cols_ + j] = b;
        }

    private:
        SizeType rows_, cols_;
        std::vector<MatrixPtrT> blocks_;
    };

    template<class Derived>
    Blocks<Derived> block2x2(
        const Tensor<Derived, 2> &a00, const Tensor<Derived, 2> &a01,
        const Tensor<Derived, 2> &a10, const Tensor<Derived, 2> &a11
        )
    {
        using MatrixPtrT = typename Blocks<Derived>::MatrixPtrT;
        std::vector<MatrixPtrT> vec = {
            make_ref(a00.derived()),
            make_ref(a01.derived()),
            make_ref(a10.derived()),
            make_ref(a11.derived()),
        };

        return Blocks<Derived>(2, 2, vec);
    }

    template<class Derived>
    Blocks<Derived, 2> block3x3(
        const Tensor<Derived, 2> &a00, const Tensor<Derived, 2> &a01, const Tensor<Derived, 2> &a02,
        const Tensor<Derived, 2> &a10, const Tensor<Derived, 2> &a11, const Tensor<Derived, 2> &a12,
        const Tensor<Derived, 2> &a20, const Tensor<Derived, 2> &a21, const Tensor<Derived, 2> &a22
        )
    {
        using MatrixPtrT = typename Blocks<Derived, 2>::MatrixPtrT;
        std::vector<MatrixPtrT> vec = {
            make_ref(a00.derived()),
            make_ref(a01.derived()),
            make_ref(a02.derived()),
            make_ref(a10.derived()),
            make_ref(a11.derived()),
            make_ref(a12.derived()),
            make_ref(a20.derived()),
            make_ref(a21.derived()),
            make_ref(a22.derived())
        };

        return Blocks<Derived, 2>(3, 3, vec);
    }

    //////////////////////////////////////////////////////////////////


    template<class Vector>
    class Blocks<Vector, 1> : public Expression<Blocks<Vector, 1>> {
    public:
        using VectorPtrT = std::shared_ptr<const Vector>;
        using SizeType = typename utopia::Traits<Vector>::SizeType;

        Blocks(const std::vector<VectorPtrT> &blocks)
        : blocks_(blocks)
        {

        }

        Blocks(const SizeType size)
        : blocks_(size)
        {
        }

        const std::vector<VectorPtrT> &blocks() const
        {
            return blocks_;
        }

        inline SizeType size() const
        {
            return blocks_.size();
        }

        const Vector &block(const SizeType i) const
        {
            assert(i < size());
            return *blocks_[i];
        }

        const VectorPtrT &block_ptr(const SizeType i) const
        {
            assert(i < size());
            return blocks_[i];
        }

        bool block_is_null(const SizeType i) const
        {
            assert(i < size());
            return !bool(blocks_[i]);
        }

        void set_block(const SizeType i, const VectorPtrT &b)
        {
            assert(i < size());
            blocks_[i] = b;
        }

    private:
        std::vector<VectorPtrT> blocks_;
    };

    template<class Derived>
    Blocks<Derived> blocks(
        const Tensor<Derived, 1> &a0,
        const Tensor<Derived, 1> &a1
        )
    {
        using VectorPtrT = typename Blocks<Derived>::VectorPtrT;
        std::vector<VectorPtrT> vec = {
            make_ref(a0.derived()),
            make_ref(a1.derived()),
        };

        return Blocks<Derived>(vec);
    }

    template<class Vector>
    void undo_blocks(
        const Tensor<Vector, 1> &block_vec,
        Tensor<Vector, 1> &a0,
        Tensor<Vector, 1> &a1
        )
    {
        const auto &block_vec_d = block_vec.derived();

        auto &a0_d = a0.derived();
        auto &a1_d = a1.derived();

        auto r = range(block_vec_d);

        auto r0 = range(a0_d);
        auto r1 = range(a1_d);

        {
            Read<Vector> r_(block_vec_d);
            Write<Vector> w0(a0_d), w1(a1_d);

            SizeType index = r.begin();
            for(auto i = r0.begin(); i < r0.end(); ++i) {
                a0_d.set(i, block_vec_d.get(index++));
            }

            for(auto i = r1.begin(); i < r1.end(); ++i) {
                a1_d.set(i, block_vec_d.get(index++));
            }
        }
    }

    template<class Derived>
    Blocks<Derived, 1> blocks(
        const Tensor<Derived, 1> &a0,
        const Tensor<Derived, 1> &a1,
        const Tensor<Derived, 1> &a2
        )
    {
        using VectorPtrT = typename Blocks<Derived, 1>::VectorPtrT;
        std::vector<VectorPtrT> vec = {
            make_ref(a0.derived()),
            make_ref(a1.derived()),
            make_ref(a2.derived()),
        };

        return Blocks<Derived, 1>(vec);
    }

    template<class Derived>
    void undo_blocks(
        const Tensor<Derived, 1> &block_vec,
        Tensor<Derived, 1> &a0,
        Tensor<Derived, 1> &a1,
        Tensor<Derived, 1> &a2
        )
    {

        auto const &block_vec_d = block_vec.derived();
        auto &a0_d = a0.derived();
        auto &a1_d = a1.derived();
        auto &a2_d = a2.derived();

        auto r = range(block_vec);

        auto r0 = range(a0_d);
        auto r1 = range(a1_d);
        auto r2 = range(a2_d);

        {
            Read<Derived> r_(block_vec);
            Write<Derived> w0(a0_d), w1(a1_d), w2(a2_d);

            SizeType index = r.begin();
            for(auto i = r0.begin(); i < r0.end(); ++i) {
                a0_d.set(i, block_vec.get(index++));
            }

            for(auto i = r1.begin(); i < r1.end(); ++i) {
                a1_d.set(i, block_vec.get(index++));
            }

            for(auto i = r2.begin(); i < r2.end(); ++i) {
                a2_d.set(i, block_vec.get(index++));
            }
        }
    }

    /////////////////////////////////////////////////////////////////


    template<class Expr>
    class Traits< Blocks<Expr> > : public Traits<Expr> {};

    template<class Tensor>
    Size size(const Blocks<Tensor, 2> &expr)
    {
        Size s(2);

        for(SizeType i = 0; i < expr.rows(); ++i) {
            for(SizeType j = 0; j < expr.cols(); ++j) {
                if(!expr.block_is_null(i, j)) {
                    s.set(0, s.get(0) + size(expr.block(i, j)).get(0));
                    break;
                }
            }
        }

        for(SizeType j = 0; j < expr.cols(); ++j) {
            for(SizeType i = 0; i < expr.rows(); ++i) {
                if(!expr.block_is_null(i, j)) {
                    s.set(1, s.get(1) + size(expr.block(i, j)).get(1));
                    break;
                }
            }
        }

        assert(s.get(0) != 0);
        assert(s.get(1) != 0);
        return s;
    }

    template<class Tensor>
    Size size(const Blocks<Tensor, 1> &expr)
    {
        Size s(1);

        for(SizeType i = 0; i < expr.size(); ++i) {
            assert(!expr.block_is_null(i));
            s.set(0, s.get(0) + size(expr.block(i)).get(0));
        }

        assert(s.get(0) != 0);
        return s;
    }
}

#endif //UTOPIA_BLOCKS_HPP
