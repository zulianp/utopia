#ifndef UTOPIA_BLOCKS_HPP
#define UTOPIA_BLOCKS_HPP

#include "utopia_ForwardDeclarations.hpp"
#include "utopia_Size.hpp"

#include <vector>
#include <memory>
#include <initializer_list>

namespace utopia {
	template<class Tensor>
	class Blocks {};

	template<class Tensor>
	class Blocks<Wrapper<Tensor, 2>> : public Expression<Blocks<Wrapper<Tensor, 2>>> {
	public:
		using MatrixT = utopia::Wrapper<Tensor, 2>;
		using MatrixPtrT = std::shared_ptr<const MatrixT>;
		using SizeType = typename utopia::Traits<MatrixT>::SizeType;

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

		const MatrixT &block(const SizeType i, const SizeType j) const
		{
			assert(i < rows());
			assert(j < cols());
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

	template<class Tensor>
	Blocks<Wrapper<Tensor, 2>> block2x2(
		const Wrapper<Tensor, 2> &a00, const Wrapper<Tensor, 2> &a01,
		const Wrapper<Tensor, 2> &a10, const Wrapper<Tensor, 2> &a11
		)
	{
		using MatrixPtrT = typename Blocks<Wrapper<Tensor, 2>>::MatrixPtrT;
		std::vector<MatrixPtrT> vec = {
			make_ref(a00),
			make_ref(a01),
			make_ref(a10),
			make_ref(a11),
		};

		return Blocks<Wrapper<Tensor, 2>>(2, 2, vec);
	}

	template<class Tensor>
	Blocks<Wrapper<Tensor, 2>> block3x3(
		const Wrapper<Tensor, 2> &a00, const Wrapper<Tensor, 2> &a01, const Wrapper<Tensor, 2> &a02,
		const Wrapper<Tensor, 2> &a10, const Wrapper<Tensor, 2> &a11, const Wrapper<Tensor, 2> &a12,
		const Wrapper<Tensor, 2> &a20, const Wrapper<Tensor, 2> &a21, const Wrapper<Tensor, 2> &a22
		)
	{
		using MatrixPtrT = typename Blocks<Wrapper<Tensor, 2>>::MatrixPtrT;
		std::vector<MatrixPtrT> vec = {
			make_ref(a00),
			make_ref(a01),
			make_ref(a02),
			make_ref(a10),
			make_ref(a11),
			make_ref(a12),
			make_ref(a20),
			make_ref(a21),
			make_ref(a22)
		};

		return Blocks<Wrapper<Tensor, 2>>(3, 3, vec);
	}

	//////////////////////////////////////////////////////////////////


	template<class Tensor>
	class Blocks<Wrapper<Tensor, 1>> : public Expression<Blocks<Wrapper<Tensor, 1>>> {
	public:
		using VectorT = utopia::Wrapper<Tensor, 1>;
		using VectorPtrT = std::shared_ptr<const VectorT>;
		using SizeType = typename utopia::Traits<VectorT>::SizeType;

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

		const VectorT &block(const SizeType i) const
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

	template<class Tensor>
	Blocks<Wrapper<Tensor, 1>> block2(
		const Wrapper<Tensor, 1> &a0,
		const Wrapper<Tensor, 1> &a1
		)
	{
		using VectorPtrT = typename Blocks<Wrapper<Tensor, 1>>::VectorPtrT;
		std::vector<VectorPtrT> vec = {
			make_ref(a0),
			make_ref(a1),
		};

		return Blocks<Wrapper<Tensor, 1>>(vec);
	}

	template<class Tensor>
	Blocks<Wrapper<Tensor, 1>> block3(
		const Wrapper<Tensor, 1> &a0,
		const Wrapper<Tensor, 1> &a1,
		const Wrapper<Tensor, 1> &a2
		)
	{
		using VectorPtrT = typename Blocks<Wrapper<Tensor, 1>>::VectorPtrT;
		std::vector<VectorPtrT> vec = {
			make_ref(a0),
			make_ref(a1),
			make_ref(a2),
		};

		return Blocks<Wrapper<Tensor, 1>>(vec);
	}

	/////////////////////////////////////////////////////////////////


	template<class Expr>
	class Traits< Blocks<Expr> > : public Traits<Expr> {};

	template<class Tensor>
	Size size(const Blocks<Wrapper<Tensor, 2>> &expr)
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
	Size size(const Blocks<Wrapper<Tensor, 1>> &expr)
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
