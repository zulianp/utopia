#ifndef UTOPIA_FE_HPP
#define UTOPIA_FE_HPP 


#include <set>
#include <assert.h>
// #include "utopia_Unary.hpp"
// #include "utopia_StoreAs.hpp"
#include "utopia.hpp"

namespace utopia {

	template<class Traits, int Backend = Traits::Backend>
	class FESpace {
	public:

	};

	template<class Traits, int Backend = Traits::Backend>
	class VectorFESpace {
	public:

	};


	template<typename Scalar_>
	class FEFunction : public Expression< FEFunction<Scalar_> > {
	public: 
		enum {
			Order = 0
		};

		typedef Scalar_ Scalar;

		std::string getClass() const { return "FEFunction"; }
	};


	template<typename Scalar_>
	class Traits< FEFunction<Scalar_> > {
	public:
		typedef Scalar_ Scalar;

		enum {
			FILL_TYPE = utopia::FillType::DENSE
		};
	};

	template<class Derived>
	class DifferentialOperator : public Expression<Derived> {
	public:
		/*virtual*/ ~DifferentialOperator() {}
	};


	template<class Expr>
	class Gradient : public DifferentialOperator< Gradient<Expr> > {
	public:
		enum {
			Order = Expr::Order + 1
		};

		typedef typename Expr::Scalar Scalar;

		std::string getClass() const { return "Gradient<" + expr_.getClass() + ">"; }

		inline const Expr &expr() const
		{
			return expr_;
		}

		Gradient(const Expr &expr)
		: expr_(expr)
		{}

	private:
		UTOPIA_STORE_CONST(Expr) expr_;
	};

	template<class Derived>
	inline Gradient<Derived> grad(const Expression<Derived> &expr) {
		return Gradient<Derived>(expr);
	}

	template<class Expr>
	class Traits< Gradient<Expr> > : public Traits<Expr> {
	public:
		enum {
			FILL_TYPE = utopia::FillType::DENSE
		};
	};


	template<class Expr>
	class Divergence : public DifferentialOperator< Divergence<Expr> > {
	public:
		enum {
			Order = Expr::Order - 1
		};

		typedef typename Expr::Scalar Scalar;

		std::string getClass() const { return "Divergence<" + expr_.getClass() + ">"; }

		inline const Expr &expr() const
		{
			return expr_;
		}

		Divergence(const Expr &expr)
		: expr_(expr)
		{}

	private:
		UTOPIA_STORE_CONST(Expr) expr_;
	};

	template<class Derived>
	inline Divergence<Derived> div(const Expression<Derived> &expr) {
		return Divergence<Derived>(expr);
	}

	template<class Expr>
	class Traits< Divergence<Expr> > : public Traits<Expr> {
	public:
		enum {
			FILL_TYPE = utopia::FillType::DENSE
		};
	};



	template<class Expr>
	class Curl : public DifferentialOperator< Curl<Expr> > {
	public:
		enum {
			Order = Expr::Order - 1
		};

		typedef typename Expr::Scalar Scalar;

		std::string getClass() const { return "Curl<" + expr_.getClass() + ">"; }

		inline const Expr &expr() const
		{
			return expr_;
		}

		Curl(const Expr &expr)
		: expr_(expr)
		{}

	private:
		UTOPIA_STORE_CONST(Expr) expr_;
	};

	template<class Derived>
	inline Curl<Derived> curl(const Expression<Derived> &expr) {
		return Curl<Derived>(expr);
	}

	template<class Expr>
	class Traits< Curl<Expr> > : public Traits<Expr> {
	public:
		enum {
			FILL_TYPE = utopia::FillType::DENSE
		};
	};


	template<class Expr_>
	class Integral : public Expression< Integral<Expr_> > {
	public:
		typedef Expr_ Expr;

		enum {
			Order = Expr::Order
		};

		typedef typename Expr::Scalar Scalar;

		std::string getClass() const { return "Integral<" + expr_.getClass() + ">"; }

		Integral(const Expr &expr)
		: expr_(expr)
		{}

		inline const Expr &expr() const
		{
			return expr_;
		}

	private:
		UTOPIA_STORE_CONST(Expr) expr_;
	};


	template<class Derived>
	inline Integral<Derived> integral(const Expression<Derived> &expr) {
		return Integral<Derived>(expr);
	}

	template<class Expr>
	class Traits< Integral<Expr> > : public Traits<Expr> {
	public:
		enum {
			FILL_TYPE = utopia::FillType::DENSE
		};
	};



	template<class Expr_, int OrderOfDifferentiation_>
	class TimeDerivative : public Expression< TimeDerivative<Expr_, OrderOfDifferentiation_> > {
	public:

		typedef Expr_ Expr;
		
		enum {
			Order = Expr::Order
		};

		enum {
			OrderOfDifferentiation = OrderOfDifferentiation_
		};

		typedef typename Expr::Scalar Scalar;

		TimeDerivative(const Expr &expr)
		: expr_(expr)
		{}

		std::string getClass() const { return "TimeDerivative<" + expr_.getClass() + ">"; }

		
		inline const Expr &expr() const
		{
			return expr_;
		}

	private:
		UTOPIA_STORE_CONST(Expr) expr_;
	};

	template<class Derived>
	TimeDerivative<Derived, 1> dt(const Expression<Derived> &expr)
	{
		return expr.derived();
	}

	template<class Expr, int OrderOfDifferentiation>
	class Traits< TimeDerivative<Expr, OrderOfDifferentiation> > : public Traits<Expr> {
	public:
		enum {
			FILL_TYPE = utopia::FillType::DENSE
		};
	};

	template<typename T>
	class DeltaT : public Number<T> {
	public:
		DeltaT(const T dt = 0.1)
		: Number<T>(dt){}
	};

	template<typename T>
	class Traits< DeltaT<T> > {
	public:
		enum {
			FILL_TYPE = utopia::FillType::SCALAR
		};
	};

	template<typename T, int Order_>
	class ConstantCoefficient : public Expression< ConstantCoefficient<T, Order_> > {
	public:

		enum {
			Order = Order_
		};

		typedef typename Traits<T>::Scalar Scalar;	

		ConstantCoefficient(const T &value)
		: value_(value)
		{}


		inline const T &expr() const
		{
			return value_;
		}

		inline operator const T&() const
		{
			return value_;
		}

		inline std::string getClass() const override
		{
			return "ConstantCoefficient" + std::to_string(Order);
		}


		//always returns value_
		inline const T &operator[](const int) const
		{
			return value_;
		}

		template<class Point, class Result>
		void eval(const Point &, Result &result) const {
			result = value_;
		}

	private:
		T value_;
	};

	template<typename T>
	class ConstantCoefficient<T, 0> : public Expression< ConstantCoefficient<T, 0> > {
	public:

		enum {
			Order = 0
		};

		typedef T Scalar;	

		ConstantCoefficient(const T &value)
		: value_(value)
		{}


		inline const T &expr() const
		{
			return value_;
		}

		inline operator const T&() const
		{
			return value_;
		}

		//always returns value_
		inline const T &operator[](const int) const
		{
			return value_;
		}

		template<class Point, class Result>
		void eval(const Point &, Result &result) const {
			result = value_;
		}

		inline std::string getClass() const override
		{
			return "ConstantCoefficient0";
		}
	private:
		T value_;
	};

	template<class T>
	inline ConstantCoefficient<T, 0> coeff(const T &expr) {
		return ConstantCoefficient<T, 0>(expr);
	}


	template<class T>
	inline ConstantCoefficient<T, 1> vec_coeff(const T &expr) {
		return ConstantCoefficient<T, 1>(expr);
	}


	template<class T>
	inline ConstantCoefficient<T, 2> mat_coeff(const T &expr) {
		return ConstantCoefficient<T, 2>(expr);
	}




	template<class T>
	class Traits< ConstantCoefficient<T, 0> > {
	public:
		enum {
			FILL_TYPE = utopia::FillType::DENSE
		};

		enum {
			Order = 0
		};

		typedef T Scalar;	
	};


	template<class T>
	class Traits< ConstantCoefficient<T, 1> > {
	public:
		enum {
			FILL_TYPE = utopia::FillType::DENSE
		};

		enum {
			Order = 1
		};

		typedef typename Traits<T>::Scalar Scalar;	
	};


	template<class T>
	class Traits< ConstantCoefficient<T, 2> > {
	public:
		enum {
			FILL_TYPE = utopia::FillType::DENSE
		};

		enum {
			Order = 2
		};

		typedef typename Traits<T>::Scalar Scalar;	
	};


	template<class Fun, typename T, int Order_>
	class FunctionCoefficient : public Expression< FunctionCoefficient<Fun, T, Order_> > {
	public:
			enum {
				Order = Order_
			};

			typedef typename Traits<T>::Scalar Scalar;	

			FunctionCoefficient(const Fun &fun)
			: fun_(fun)
			{}


			template<class Point, class Result>
			void eval(const Point &p, Result &result) const {
				fun_(p, result);
			}

			inline std::string getClass() const override
			{
				return "FunctionCoefficient" + std::to_string(Order);
			}

			inline const Fun &fun() const
			{
				return fun_;
			}

			inline Fun &fun()
			{
				return fun_;
			}

		private:
			Fun fun_;
	};

	template<class Fun, typename T>
	class FunctionCoefficient<Fun, T, 0>  : public Expression< FunctionCoefficient<Fun, T, 0> > {
	public:
			enum {
				Order = 0
			};

			typedef T Scalar;	

			FunctionCoefficient(const Fun &fun)
			: fun_(fun)
			{}


			template<class Point, class Result>
			void eval(const Point &p, Result &result) const {
				result = fun_(p);
			}


			inline std::string getClass() const override
			{
				return "FunctionCoefficient0";
			}

			inline const Fun &fun() const
			{
				return fun_;
			}

			inline Fun &fun()
			{
				return fun_;
			}

		private:
			Fun fun_;
	};

	template<class Fun, class T, int Order_>
	class Traits< FunctionCoefficient<Fun, T, Order_> > {
	public:
		enum {
			FILL_TYPE = utopia::FillType::DENSE
		};

		enum {
			Order = Order_
		};

		typedef typename Traits<T>::Scalar Scalar;	
	};

	template<class Fun, class T>
	class Traits< FunctionCoefficient<Fun, T, 0> > {
	public:
		enum {
			FILL_TYPE = utopia::FillType::DENSE
		};

		enum {
			Order = 0
		};

		typedef T Scalar;	
	};

	template<class Point, class T>
	inline FunctionCoefficient<T (*)(const Point&), T, 0> coeff(T (*fun)(const Point&)) {
		return fun;
	}

	template<class Point, class T>
	inline FunctionCoefficient<std::function<T(const Point&)> , T, 0> coeff(std::function<T(const Point&)> fun) {
		return fun;
	}

	template<class Point, class T>
	inline FunctionCoefficient<void (*)(const Point&, T &), T, 1> vec_coeff(void (*fun)(const Point&, T &)) {
		return fun;
	}

	template<class Point, class T>
	inline FunctionCoefficient<std::function<void(const Point&, T &)> , T, 1> vec_coeff(std::function<void(const Point&, T&)> fun) {
		return fun;
	}

	template<class Point, class T>
	inline FunctionCoefficient<void (*)(const Point&, T &), T, 2> mat_coeff(void (*fun)(const Point&, T &)) {
		return fun;
	}

	template<class Coefficient, class FESpace>
	class Interpolate : public Expression< Interpolate<Coefficient, FESpace> > {
	public:

		typedef typename Traits<Coefficient>::Scalar Scalar;

		enum {
			Order = Traits<Coefficient>::Order
		};

		Interpolate(const Coefficient &coeff, FESpace &fe)
		: coeff_(coeff), fe_(fe)
		{}


	private:
		Coefficient coeff_;
		UTOPIA_STORE(FESpace) fe_;	
	};

	template<class Coefficient, class FESpace>
	class Traits< Interpolate<Coefficient, FESpace> > : public Traits<Coefficient> {
	public:
	
	};

		

	template<class Coefficient, class FESpace>
	inline Interpolate<Coefficient, FESpace> interpolate(const Coefficient &coeff, FESpace &fe)
	{
		return Interpolate<Coefficient, FESpace>(coeff, fe);
	}


	template<class Coefficient, class FESpace, class Tensor>
	inline Interpolate<Coefficient, FESpace> interpolate(const Coefficient &coeff, FESpace &fe, Tensor &&tensor)
	{
		return Interpolate<Coefficient, FESpace>(coeff, fe, tensor);
	}

	template<class FEFun, int N>
	class VectorFE : public Expression< VectorFE<FEFun, N> > {
	public:
		typedef typename utopia::Traits<FEFun>::Scalar Scalar;

		enum {
			Order = 1
		};

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

#include "utopia_LibMeshBackend.hpp"


#endif	//UTOPIA_FE_HPP
