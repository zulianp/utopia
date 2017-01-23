#ifndef UTOPIA_CL_PROGRAM_HPP
#define UTOPIA_CL_PROGRAM_HPP

#include "utopia_CLKernel.hpp"
// #include "utopia_Expressions.hpp"


namespace utopia {
	namespace opencl {
		
		template<class Left, class Right, int Order = (Left::Order == 2 || Right::Order == 2)? 2 : 1>
		class MultiplyTransform {
		public:
			static const Multiply<Left, Right> & transform(const Multiply<Left, Right> &expr)
			{
				return expr;
			}
		};
		
		template<class Left, class Right>
		class MultiplyTransform<Transposed<Left> , Right, 1> {
		public:
			typedef utopia::Reduce<utopia::Binary<Left, Right, utopia::EMultiplies>, utopia::Plus> Type;
			
			static const Type &transform(const Multiply<Left, Right> &expr)
			{
				return Reduce<Binary<Left, Right, EMultiplies>, Plus>(Binary<Left, Right, EMultiplies>(expr.left().expr(), expr.right()));
			}
		};
		
		template<class Expr>
		class IsTransposed {
		public:
			static const bool value = false;
		};
		
		template<class Expr>
		class IsTransposed< utopia::Transposed<Expr> > {
		public:
			static const bool value = !IsTransposed<Expr>::value;
		};
		
		template<class Expr, class Operation>
		class IsTransposed< utopia::Unary<Expr, Operation> > {
		public:
			static const bool value = IsTransposed<Expr>::value;
		};
		
		template<class Expr>
		constexpr bool is_transposed()
		{
			return IsTransposed<Expr>::value;
		}
		
		
		template<class _Expr>
		class MultiKernelTransform {
		public:
			typedef _Expr Type;
			
			static const _Expr &transform(const _Expr &expr)
			{
				return expr;
			}
			
			static _Expr &transform(_Expr &expr)
			{
				return expr;
			}
		};
		
		template<class Expr>
		inline typename MultiKernelTransform<Expr>::Type mk_transform(const Expr &expr)
		{
			return MultiKernelTransform<Expr>::transform(expr);
		}
		
		template<class Tensor, int Order>
		inline const Wrapper<Tensor, Order> & mk_transform(const Wrapper<Tensor, Order> &expr)
		{
			return expr;
		}

		template<class Expr> 
		class EvaluateType {
		public:
			typedef utopia::Evaluate<Expr> Type;
		};

		template<class Tensor, int Order> 
		class EvaluateType< Wrapper<Tensor, Order> > {
		public:
			typedef utopia::Wrapper<Tensor, Order> Type;
		};
		
		
		template<class Left, class Right, int LeftOrder = Left::Order, int RightOrder = Right::Order>
		class MultiKernelTransformMultiplyAux {
		public:
			typedef typename MultiKernelTransform<Left>::Type  TLeft;
			typedef typename MultiKernelTransform<Right>::Type TRight;

			typedef typename EvaluateType<TLeft>::Type  ELeft;
			typedef typename EvaluateType<TRight>::Type ERight;

			
			//transform it to the evaluate expression
			// typedef utopia::Evaluate< utopia::Multiply<ELeft, ERight>  > Type;
			typedef utopia::Multiply<ELeft, ERight> Type;
			
			static Type transform(const Multiply<Left, Right > &expr)
			{
				return 
				// make_evaluate(
					make_evaluate( mk_transform(expr.left() ) ) * 
					make_evaluate( mk_transform(expr.right()) ) 
				// )
				;
			}
		};
		
		template<class Left, class Right>
		class MultiKernelTransformMultiplyAux< Transposed<Left>, Right, 1, 1> {
		public:
			typedef utopia::opencl::MultiKernelTransform< typename MultiplyTransform<Transposed<Left>, Right>::Type > Transform;
			typedef typename Transform::Type Type;
			
			static Type transform(const Multiply<Transposed<Left>, Right > &expr)
			{
				return Type(dot( mk_transform(expr.left().expr()),
								mk_transform(expr.right()) ));
			}
		};
		
		template<class Left, class Right>
		class MultiKernelTransform< Multiply<Left, Right > > {
		public:
			
			typedef typename utopia::opencl::MultiKernelTransformMultiplyAux<Left, Right>::Type Type;
			
			static Type transform(const Multiply<Left, Right > &expr)
			{
				return MultiKernelTransformMultiplyAux<Left, Right>::transform(expr);
			}
			
			// typedef typename MultiKernelTransform<Left>::Type  TLeft;
			// typedef typename MultiKernelTransform<Right>::Type TRight;
			
			// //transform it to the evaluate expression
			// typedef utopia::Multiply<
			// 				 Evaluate<TLeft>,
			// 				 Evaluate<TRight>
			// 				 > Type;
			
			// static Type transform(const Multiply<Left, Right > &expr)
			// {
			// 	return Type(make_evaluate(mk_transform(expr.left())), make_evaluate(mk_transform(expr.right())));
			// }
		};
		
		template<class Left, class Right, class Operation>
		class MultiKernelTransform< Binary<Left, Right, Operation> > {
		public:
			typedef utopia::Binary<Left, Right, Operation> Expr;
			typedef typename MultiKernelTransform<Left>::Type TLeft;
			typedef typename MultiKernelTransform<Right>::Type TRight;
			
			//transform it to the evaluate expression
			typedef utopia::Binary<
			TLeft,
			TRight,
			Operation
			> Type;
			
			static Type transform(const Expr &expr)
			{
				return Type(mk_transform(expr.left()), mk_transform(expr.right()), expr.operation());
			}
		};

		template<class Left, class Right>
		class MultiKernelTransform< Construct<Left, Right> > {
		public:
			typedef utopia::Construct<Left, Right> Expr;
			typedef typename MultiKernelTransform<Right>::Type TRight;
			
			//transform it to the evaluate expression
			typedef utopia::Construct<
			Left,
			TRight
			> Type;
			
			static Type transform(const Expr &expr)
			{
				return Type( expr.left(), mk_transform(expr.right()) );
			}
		};
		
		template<class InnerExpr, class Operation>
		class MultiKernelTransform<Unary<InnerExpr, Operation> > {
		public:
			//transform it to the evaluate expression
			typedef utopia::Unary<typename MultiKernelTransform<InnerExpr>::Type, Operation> Type;
			
			static Type transform(const Unary<InnerExpr, Operation>  &expr)
			{
				return Type(mk_transform(expr.expr()));
			}
		};
		
		template<class InnerExpr, class Operation>
		class MultiKernelTransform< Reduce<InnerExpr, Operation> > {
		public:
			typedef typename MultiKernelTransform<InnerExpr>::Type InnerExprT;
			typedef utopia::Evaluate< Reduce<InnerExprT, Operation> > Type;
			
			static Type transform(const Reduce<InnerExpr, Operation> &expr)
			{
				return Type(mk_transform(expr.expr()));
			}
		};
		
		template<class InnerExpr>
		class MultiKernelTransform< Norm<InnerExpr, 2> > {
		public:
			typedef typename utopia::Unfold< Norm<InnerExpr, 2> >::Type  Unfolded;
			
			typedef utopia::opencl::MultiKernelTransform<Unfolded> UnfoldedTransform;
			typedef Evaluate<typename UnfoldedTransform::Type> Type;
			
			static Type transform(const Norm<InnerExpr, 2> &expr)
			{
				return UnfoldedTransform::transform( shallow_unfold(expr) );
			}
		};
		
		template<class InnerExpr>
		class MultiKernelTransform< Transposed<InnerExpr> > {
		public:
			typedef utopia::Transposed<typename MultiKernelTransform<InnerExpr>::Type> Type;
			
			static Type transform(const Transposed<InnerExpr>  &expr)
			{
				return Type(mk_transform(expr.expr()));
			}
		};


		///FIXME 
		template<class Left, class Right>
		class MultiKernelTransform< Assign<Left, Right> > {
		public:
			typedef typename MultiKernelTransform<Left>::Type  LeftT;
			typedef typename MultiKernelTransform<Right>::Type RightT;

			typedef utopia::Construct<LeftT, RightT> Type;

			static Type transform(const Assign<Left, Right> &expr)
			{
				return Type(expr.left(), 
							mk_transform(expr.right()) );
			}

		};
		
		template<class Expr>
		class Program {
		public:
			std::vector<IKernel *> kernels_;
			
			class EvalAction {
			public:
				EvalAction(std::vector<IKernel *>  &kernels) : kernels_(kernels), next_kernel_(0) {}

				template<class Tensor>
				static void resize(const Size &size, Wrapper<Tensor, 1> &w)
				{
					raw_type(w).resize(size.get(0));
				}

				template<class Tensor>
				static void resize(const Size &size, Wrapper<Tensor, 2> &w)
				{
					raw_type(w).resize(size.get(0), size.get(1));
				}

				// template<class Tensor>
				// static void resize(const Size &size, Wrapper<Tensor, 0> &w)
				// {
				// 	assert(false && "should never be called"); 
				// }

				template<class EvaluatedExpr>
				void post_order_visit(const Evaluate<EvaluatedExpr> &expr)
				{
					auto w = wrap<EvaluatedExpr::Order>(expr.backend_tensor_ptr());
					resize(size(expr.expr()), w);
					post_order_visit( make_evaluate( construct(w, expr.expr()) ));
				}

				template<class Left, class Right>
				void post_order_visit(const Evaluate< Construct<Left, Right> > &expr)
				{
					typedef utopia::Construct<Left, Right>  ExprT;
					auto k = static_cast< Kernel<ExprT> * >( kernels_[next_kernel_] );
					k->call(expr.expr());
					++next_kernel_;	
				}

				template<class InnerExpr, class Operation>
				void post_order_visit(const Evaluate< Reduce<InnerExpr, Operation>, 0> &expr)
				{
					typedef utopia::Reduce<InnerExpr, Operation> ExprT;
					auto k = static_cast< Kernel<ExprT> * >( kernels_[next_kernel_] );	
					// k->call(expr.expr());
					k->call(expr);
					++next_kernel_;	
				}

				// template<typename T, class InnerExpr, class Operation>
				// void post_order_visit(const Evaluate< Construct<Number<T>, Evaluate< Reduce<InnerExpr, Operation> > > > &expr)
				// {
				// 	typedef utopia::Reduce<InnerExpr, Operation> ExprT;
				// 	auto k = static_cast< Kernel<ExprT> * >( kernels_[next_kernel_] );	
				// 	k->call(expr);
				// 	++next_kernel_;	
				// }

				template<class Ignore>
				inline void post_order_visit(const Expression<Ignore> &) {}

				template<class Ignore>
				inline void pre_order_visit(const Ignore &) {}

				template<class Ignore>
				inline void in_order_visit(const Ignore &) {}

			private:
				std::vector<IKernel *>  &kernels_;
				int next_kernel_;
			};


			template<class Derived>
			inline bool execute(const Expression<Derived> &all_expr)
			{
				const auto &expr = all_expr.derived();

				EvalAction ea(kernels_);
				TreeNavigator< EvalAction & > nav(ea);
				auto eval = make_evaluate(mk_transform(expr));
				nav.visit(eval);
				return false;
			}


			template<typename T, class Right>
			inline bool execute(const Construct< Number<T>, Right> &expr)
			{
				EvalAction ea(kernels_);
				TreeNavigator< EvalAction & > nav(ea);
				auto eval = mk_transform(expr.right());
				nav.visit(eval);
				expr.left() = eval.get_value();
				return false;
			}


			// template<typename T, class Right>
			// void initialize(const Construct< Number<T>, Right> &expr)
			// {
			// 	initialize(expr.right());
			// }
			
			// template<class Derived>
			// void initialize(const Expression<Derived> &all_expr)

			void initialize(const Expr &expr)
			{
				if(initialized_) return;

				// const auto &expr = all_expr.derived();

				kernels_.clear();
				
				Env env;
				Options options;
				TreeNavigator< Program & > nav(*this);
				
				auto eval = make_evaluate(mk_transform(expr));
				nav.visit(eval);
				

				cl::Program::Sources sources;
				for(auto k_ptr : kernels_) {
					// write(k_ptr->get_name() + ".cl", 
					// 	  k_ptr->get_code_string());

					if(!k_ptr->is_callable()) {
						sources.push_back({ k_ptr->get_code_string().c_str(), 
									    	k_ptr->get_code_string().size() });
					}
				}

				// std::cout << "-------------------------------------------" << std::endl;
				// std::cout << "compiling program" << std::endl;
			
				std::string flags = " -I" + utopia::Utopia::Instance().get("opencl_templates_path") + "/../kernels/ ";
				// std::cout << flags << std::endl;

				compile_opencl_programs(sources, program_, flags);
				initialized_ = true;

				// std::cout << "-------------------------------------------" << std::endl;

				//init kernels
				for(auto k_ptr : kernels_) {
					if(!k_ptr->is_callable()) {
						k_ptr->make_callable(program_);
					} else {
						// std::cout << "already callable" << std::endl;
					}
				}
			}
			
			template<class Tensor, int Order>
			void post_order_visit(const Wrapper<Tensor, Order> &expr)
			{
				// std::cout << "Data" << std::endl;
			}
			
			template<class EvaluatedExpr>
			void post_order_visit(const Evaluate<EvaluatedExpr> &expr)
			{
				auto w = wrap<EvaluatedExpr::Order>(expr.backend_tensor_ptr());
				kernels_.push_back( &get_kernel( construct(w, expr.expr()) ) );
			}

			template<class Left, class Right>
			void post_order_visit(const Evaluate< Construct<Left, Right> > &expr)
			{
				kernels_.push_back( &get_kernel( expr.expr() ) );
			}

			template<class InnerExpr, class Operation>
			void post_order_visit(const Evaluate< Reduce<InnerExpr, Operation> > &expr)
			{
				kernels_.push_back( &get_kernel( expr.expr() ) );
			}
			
			template<class Derived>
			void post_order_visit(const Expression<Derived> &expr)
			{
				// std::cout << "expression" << std::endl;
			}

			template<class Ignore>
			inline void pre_order_visit(const Ignore &){}

			template<class Ignore>
			inline void in_order_visit(const Ignore &) {}
			
			static Program &instance()
			{
				static Program instance;
				return instance;
			}
			
		private:
			cl::Program program_;
			bool initialized_;
			Program() : initialized_(false) {}
		};
		
		template<class Expr>
		Program<Expr> &get_program(const Expr &expr)
		{
			Program<Expr>::instance().initialize(expr);
			return Program<Expr>::instance();
		}
	}
}

#endif //UTOPIA_CL_PROGRAM_HPP

