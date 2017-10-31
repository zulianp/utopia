#ifndef UTOPIA_FE_LANG_HPP
#define UTOPIA_FE_LANG_HPP

#include "utopia.hpp"

#include <set>
#include <cassert>

namespace utopia {

	template<class Traits, int Backend = Traits::Backend>
	class FESpace {
	public:

	};

	template<class Traits, int Backend = Traits::Backend>
	class VectorFESpace {
	public:

	};


	// template<typename Scalar_>
	// class FEFunction : public Expression< FEFunction<Scalar_> > {
	// public: 
	// 	enum {
	// 		Order = 0
	// 	};

	// 	typedef Scalar_ Scalar;

	// 	std::string getClass() const { return "FEFunction"; }
	// };


	// template<typename Scalar_>
	// class Traits< FEFunction<Scalar_> > {
	// public:
	// 	typedef Scalar_ Scalar;

	// 	enum {
	// 		FILL_TYPE = utopia::FillType::DENSE
	// 	};
	// };









	template<class FEFun, int N>
	class VectorFE : public Expression< VectorFE<FEFun, N> > {
	public:
		typedef typename utopia::Traits<FEFun>::Scalar Scalar;

		enum {
			Order = 1
		};

		constexpr static int size()
		{
			return N;
		}

		VectorFE(std::initializer_list<FEFun *> funs) 
		{
			assert(N == funs.size());

			int i = 0;
			for(auto fe : funs) {
				funs_[i++] = fe;
			}
		}

		const FEFun &get(int index) const
		{
			assert(index < N);
			return *funs_[index];
		}

		FEFun &get(int index) 		{
			assert(index < N);
			return *funs_[index];
		}

	private:
		FEFun * funs_[N];
	};

	template<class FEFun, int N>
	class Traits< VectorFE<FEFun, N> > : public Traits<FEFun> {};

	template<class FEFun>
	inline VectorFE<FEFun, 2> prod(FEFun &fe_1, FEFun &fe_2)
	{
		return VectorFE<FEFun, 2>({&fe_1, &fe_2});
	}

	template<class FEFun>
	inline VectorFE<FEFun, 3> prod(FEFun &fe_1, FEFun &fe_2, FEFun &fe_3)
	{
		return VectorFE<FEFun, 3>({&fe_1, &fe_2, &fe_3});
	}

	template<class FEFun>
	inline VectorFE<FEFun, 4> prod(FEFun &fe_1, FEFun &fe_2, FEFun &fe_3, FEFun &fe_4)
	{
		return VectorFE<FEFun, 4>({&fe_1, &fe_2, &fe_3, &fe_4});
	}

	template<class ScalarSpace, int ROWS, int COLS>
	class TensorProductFESpace {
	public:
		template<class... Args>
		TensorProductFESpace(Args&... args)
		{
			for(int i = 0; i < ROWS*COLS; ++i) {
				spaces_[i] = std::make_shared<ScalarSpace>(args...);
			}
		}

		ScalarSpace &get(const int index)
		{
			assert(index < ROWS * COLS);
			assert(spaces_[index]);
			return *spaces_[index];
		}

		const ScalarSpace &get(const int index) const
		{
			assert(index < ROWS * COLS);
			assert(spaces_[index]);
			return *spaces_[index];
		}

	private:
		std::shared_ptr<ScalarSpace> spaces_[ROWS * COLS];
	};

	template<class FEFun, int ROWS, int COLS>
	class MatrixFE : public Expression< MatrixFE<FEFun, ROWS, COLS> > {
	public:
		typedef typename utopia::Traits<FEFun>::Scalar Scalar;

		enum {
			Order = 2
		};

		enum {
			FILL_TYPE = utopia::FillType::DENSE
		};

		template<class ScalarSpace>
		MatrixFE(const TensorProductFESpace<ScalarSpace, ROWS, COLS> &space)
		{
			for(int i = 0; i < ROWS*COLS; ++i) {
				funs_[i] = std::make_shared<FEFun>(space.get(i));
			}
		}

		MatrixFE(std::initializer_list< std::shared_ptr<FEFun> > fes) 
		{
			assert(COLS * COLS == fes.size());

			int i = 0;
			for(auto fe : fes) {
				funs_[i++] = fe;
			}
		}

		const FEFun &get(const int r, const int c) const
		{
			assert(r < ROWS);
			assert(c < COLS);
			assert(funs_[r * COLS + c]);
			return *funs_[r * COLS + c];
		}

		const FEFun &get(const int index) const
		{
			assert(index < ROWS * COLS);
			assert(funs_[index]);
			return *funs_[index];
		}

		FEFun &get(const int index)
		{
			assert(index < ROWS * COLS);
			assert(funs_[index]);
			return *funs_[index];
		}


		static constexpr int size() 
		{
			return ROWS * COLS;
		}

	private:
		std::shared_ptr<FEFun> funs_[ROWS * COLS];
	};

	template<class FEFun, int ROWS, int COLS>
	class Traits< MatrixFE<FEFun, ROWS, COLS> > : public Traits<FEFun> {
	public:
		enum {
			FILL_TYPE = utopia::FillType::DENSE
		};
	};

	template<class FEFun>
	inline MatrixFE<FEFun, 2, 2> tensor_product_fe_2x2(FEFun &fe_1, FEFun &fe_2, FEFun &fe_3, FEFun &fe_4)
	{
		return MatrixFE<FEFun, 2, 2>({ make_ref(fe_1), make_ref(fe_2), 
									   make_ref(fe_3), make_ref(fe_4) });
	}

	template<class FEFun>
	inline MatrixFE<FEFun, 3, 3> tensor_product_fe_3x3(FEFun &fe_1, FEFun &fe_2, FEFun &fe_3, 
													   FEFun &fe_4, FEFun &fe_5, FEFun &fe_6,
													   FEFun &fe_7, FEFun &fe_8, FEFun &fe_9
													   )
	{
		return MatrixFE<FEFun, 3, 3>({ make_ref(fe_1), make_ref(fe_2), make_ref(fe_3), 
									   make_ref(fe_4), make_ref(fe_5), make_ref(fe_6), 
									   make_ref(fe_7), make_ref(fe_8), make_ref(fe_9) });
	}

	template<class Traits>
	inline auto tensor_product_fe_2x2(FESpace<Traits> &fe_1, 
									  FESpace<Traits> &fe_2, 
									  FESpace<Traits> &fe_3, 
									  FESpace<Traits> &fe_4) -> MatrixFE<decltype(fe_function(fe_1)), 2, 2>
	{
		using std::make_shared;
		typedef decltype(fe_function(fe_1)) FEFun;
		return MatrixFE<FEFun, 2, 2>({ make_shared<FEFun>(fe_1), make_shared<FEFun>(fe_2), 
									   make_shared<FEFun>(fe_3), make_shared<FEFun>(fe_4) });
	}

	template<class FESpace, int ROWS, int COLS>
	MatrixFE<typename FESpace::Function, ROWS, COLS> fe_function(const TensorProductFESpace<FESpace, ROWS, COLS> &space)
	{
		return MatrixFE<typename FESpace::Function, ROWS, COLS>(space);
	}

	template<class Expr>
	class DirichletBoundaryCondition {
	public:
		DirichletBoundaryCondition(
			const Expr &expr, 
			std::initializer_list<int> boundary_tags)
		: expr_(expr)
		{
			boundary_tags_.insert(boundary_tags.begin(), boundary_tags.end());
		}

		inline Expr &expr()
		{
			return expr_;
		}

		inline const Expr &expr() const
		{
			return expr_;
		}

		const std::set<int> &boundary_tags() const
		{
			return boundary_tags_;
		}

	private:
		Expr expr_;
		std::set<int> boundary_tags_;
	};

	template<class Left, class Right>
	DirichletBoundaryCondition<Equality<Left, Right> > boundary_conditions(const Equality<Left, Right> &expr, std::initializer_list<int> boundary_tags)
	{
		return DirichletBoundaryCondition<Equality<Left, Right> > (expr, boundary_tags);
	}

}


#endif //UTOPIA_FE_LANG_HPP
