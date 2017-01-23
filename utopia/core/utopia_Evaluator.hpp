//
// Created by Patrick Zulian on 15/05/15.
//

#ifndef utopia_utopia_EVALUATOR_HPP
#define utopia_utopia_EVALUATOR_HPP

#include <vector>

#include "utopia_ForwardDeclarations.hpp"

#include "utopia_Operators.hpp"
#include "utopia_Traits.hpp"
#include "utopia_Reduce.hpp"
#include "utopia_Backend.hpp"
#include "utopia_Assign.hpp"
#include "utopia_Structured.hpp"
#include "utopia_Factory.hpp"
#include "utopia_Ranged.hpp"
#include "utopia_Multiply.hpp"
#include "utopia_Transposed.hpp"
#include "utopia_Boolean.hpp"
#include "utopia_OuterProduct.hpp"
#include "utopia_Norm.hpp"
#include "utopia_FillTypeQuery.hpp"
#include "utopia_MPI.hpp"


#include "utopia_Eval.hpp"

namespace utopia {
    template<class Result, int BAKEND_FLAG = Traits<Result>::Backend>
    class Evaluator {
    public:
        typedef utopia::Traits<Result> Traits;
        typedef typename Traits::Scalar  Scalar;

        typedef utopia::Backend<Scalar, Traits::Backend> BackendT;

        template<class Expr, int Order>
        static bool eval(const Wrapper<Expr, Order> &expr, Size &size)
        {
            return BackendT::Instance().size(eval(expr), size);
        }

// #define USE_NEW_EVAL
// #ifdef USE_NEW_EVAL        

        template<class Derived>
        static auto eval(const Expression<Derived> &expr) -> decltype(Eval<Derived, Traits>::apply(expr.derived()))
        {
            return Eval<Derived, Traits>::apply(expr.derived());
        }

        template<class Left, class Right>
        static void eval(const Construct<Left, Right> &expr)
        {
            Eval<Construct<Left, Right>, Traits>::apply(expr);
        }

        template<class Left, class Right>
        static void eval(const Assign<Left, Right> &expr)
        {
            Eval<Assign<Left, Right>, Traits>::apply(expr);
        }

        template<class Left, class Right, class Operation>
        static void eval(const InPlace<Left, Right, Operation> &expr)
        {
            Eval<InPlace<Left, Right, Operation>, Traits>::apply(expr);
        }

// #else        

//         enum {
//             FILL_TYPE = Traits::FILL_TYPE
//         };

//         int n_recursions;

//         Evaluator() 
//         : n_recursions(0)
//         {}

//         class ScopedRecursionCounter {
//         public:
//             ScopedRecursionCounter(int &val)
//             : val_(val)
//             {
//                 ++val_;
//             }

//             ~ScopedRecursionCounter()
//             {
//                 --val_;
//             }

//             bool good() 
//             {
//                 return val_ == 1;
//             }

//             int &val_;
//         };

//         template<class Derived>
//         Result eval(const Expression<Derived> &expr)
//         {   
//             ScopedRecursionCounter rc(n_recursions);
           
//             if(!rc.good()) {
//                 if(mpi_world_rank() == 0) {
//                     std::cout << "\n\n\n" << std::endl;
//                     std::cout << "Expression not matched. Send the following expression to patrick.zulian@gmail.com :" << std::endl;
//                     std::cout << "-------------------------------------" << std::endl;
//                     std::cout << tree_format(expr.getClass()) << std::endl;
//                     std::cout << "-------------------------------------" << std::endl;
//                     std::cout << "If it is a composite expression, then as a temporary solution you can evaluate sub expression separately." << std::endl;
//                     std::cout << "\n\n\n" << std::endl;
//                 }

//                 assert(false  && "expression not matched.");
//                 return Result();
//             }

//             return eval(expr.derived());
//         }

   
        
//         template<class Left, class Right>
//         Result eval(const OuterProduct<Left, Right> &expr)
//         {
//             Result result;
//             BackendT::Instance().outer(eval(expr.left()), eval(expr.right()), result);
//             return result;
//         }

//         template<class Left, class Right>
//         void eval(const Assign<Left, Right> &assign)
//         {
//             BackendT::Instance().assign( assign.left().implementation(), eval(assign.right()) );
//         }

//         template<class Left, class Right>
//         void eval(const Assign<View<Left>, Right> &assign)
//         {
//             const auto &l = assign.left();
//             BackendT::Instance().assignToRange( l.expr().implementation(), eval(assign.right()), l.rowRange(), l.colRange() );
//         }

//         template<class Left, class Right>
//         void eval(const Assign<View< Wrapper<Left, 1> >, Right> &assign)
//         {
//             const auto &l = assign.left();
//             BackendT::Instance().assignToRange( l.expr().implementation(), eval(assign.right()), l.rowRange(), l.colRange());
//         }

//         template<class Left, class Right>
//         void eval(const Construct<Left, Right> &assign)
//         {
//             BackendT::Instance().assign( assign.left().implementation(), eval(assign.right()) );
//         }

//         template<class Left, class Right>
//         void eval(const Construct< Number<Left>, Right> &assign)
//         {
//            assign.left() = eval(assign.right());
//         }

//         template<class Left, class Right>
//         void eval(const Construct<Left, Transposed <Wrapper<Right, 2> > > &assign)
//         {
//             BackendT::Instance().assignTransposed( assign.left().implementation(), eval(assign.right().expr()) );
//         }

//         template<class Left, class Right>
//         void eval(const Assign<Left, Transposed <Wrapper<Right, 2> > > &assign)
//         {
//             BackendT::Instance().assignTransposed( assign.left().implementation(), eval(assign.right().expr()) );
//         }
		
// 		template<class Left, class Right>
// 		void eval(const Assign<Left, View<Right> > &constr)
// 		{

// 			BackendT::Instance().assignFromRange( constr.left().implementation(), eval(constr.right().expr()), constr.right().rowRange(), constr.right().colRange() );
// 		}


//         ////////////////////////////////////// DIAG ////////////////////////////////////////

//         // template<typename ScalarT>
//         // Diag<Vector> eval(const Binary< Number<ScalarT>, Factory<Identity, 2>, Multiplies> &expr)
//         // {
//         //     return diag( eval(values(expr.right().size().get(0), expr.left()) ));
//         // }

//         template<class Derived>
//         EXPR_TYPE(Traits, Diag<Derived>) eval(const Diag<Derived> & /*expr*/)
//         {
//             assert(false); //TODO
//             return EXPR_TYPE(Traits, Diag<Derived>)();
//         }

//         template<class Left, class Right>
//         EXPR_TYPE(Traits, Left) eval(const Multiply< Left, Diag<Right> > &expr)
//         {
//             static_assert(Right::Order == 1, "Right has to be a vector");
//             EXPR_TYPE(Traits, Left) result; 
//             BackendT::Instance().diag_scale_right(eval(expr.left()), eval(expr.right().expr()), result);  
//             return result; 
//         }

//         template<class Left, class Right> 
//         EXPR_TYPE(Traits, Right) eval(const Multiply< Diag<Right>, Right> &expr)
//         {
//             EXPR_TYPE(Traits, Right) result;
//             BackendT::Instance().diag_scale_left(eval(expr.left().expr()), eval(expr.right()), result);
//             return result;   
//         }

//         //// assign

//         template<class Left, class Right>
//         void eval(const Assign<Left, Diag<Right> > &assign)
//         {
//             BackendT::Instance().diag( assign.left().implementation(), eval(assign.right().expr()) );
//         }

//         template<class Left, class Right>
//         void eval(const Assign<Wrapper<Left, 2>, Diag< Diag<Right> > > &assign)
//         {
//             BackendT::Instance().diag(assign.left().implementation(), eval(assign.right().expr().expr()) );
//         }

//         template<class Left, class Right>
//         void eval(const Assign<Wrapper<Left, 1>, Diag< Diag<Right> > > &assign)
//         {
//             BackendT::Instance().assign(eval(assign.left()), eval(assign.right().expr().expr()));
//         }

//         //// construct 

//         template<class Left, class Right>
//         void eval(const Construct<Left, Diag<Right> > &construct)
//         {
//             BackendT::Instance().diag( construct.left().implementation(), eval(construct.right().expr()) );
//         }

//         template<class Left, class Right>
//         void eval(const Construct<Wrapper<Left, 2>, Diag< Diag<Right> > > &construct)
//         {
//             BackendT::Instance().diag( construct.left().implementation(), eval(construct.right().expr().expr()) );
//         }

//         template<class Left, class Right>
//         void eval(const Construct<Wrapper<Left, 1>, Diag< Diag<Right> > > &assign)
//         {
//             BackendT::Instance().assign(eval(assign.left()), eval(assign.right().expr().expr()));
//         }


//         ////////////////////////////////////// DIAG ////////////////////////////////////////
//         template<class Left, class Right>
//         void eval(const Construct<Left, View<Right> > &constr)
//         {
//             const auto &r = constr.right();
//             BackendT::Instance().assignFromRange( constr.left().implementation(), eval(r.expr()), r.rowRange(), r.colRange() );
//         }

//         template<class Left, class Right, int Order>
//         void eval(const Assign<View<Left>, Factory<Right, Order> > &assign)
//         {
//             const auto &l = assign.left();
//             BackendT::Instance().assignToRange( l.expr().implementation(), assign.right().type(),  l.rowRange(), l.colRange() );
//         }

//         template<class Type, int Order>
//         typename TypeAndFill<Traits, Factory<Type, Order> >::Type eval(const Factory<Type, Order> & expr)
//         {
//             typename TypeAndFill<Traits, Factory<Type, Order> >::Type ret;
//             BackendT::Instance().build(ret, expr.size(), expr.type());
//             return ret;
//         }

//         template<class Expr, int Order = Differentiable<Expr>::Order>
//         EXPR_TYPE(Traits, Expr) eval(const Differentiable<Expr> &expr)
//         {
//             return eval(expr.expr());
//         }

//         template<class Tensor, int Order>
//         const Tensor & eval(const Differentiable< Wrapper<Tensor, Order> > &expr)
//         {
//             return eval(expr.expr());
//         }

//         template<class Tensor, int Order>
//         const Tensor & eval(const Differentiable<const Wrapper<Tensor, Order> > &expr)
//         {
//             return eval(expr.expr());
//         }

//         template<class Left, class Right, int Order>
//         void eval(const Construct<Left, Factory<Right, Order> > &assign)
//         {
//             BackendT::Instance().build( assign.left().implementation(), assign.right().size(), assign.right().type() );
//         }

//         template<class Left, class Right, int Order>
//         void eval(const Assign<Left, Factory<Right, Order> > &assign)
//         {
//             BackendT::Instance().build( assign.left().implementation(), assign.right().size(), assign.right().type() );
//         }

// //        template<class Left, class Right, class Operation>
// //        Result eval(const Binary<Left, Right, Operation> &expr)
// //        {
// //            return BackendT::Instance().apply( eval(expr.left()), eval(expr.right()), expr.operation() );
// //        }


//     ////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//     ////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//     //////////////////////////////////////////////////AXPY////////////////////////////////////////////////////////


//         template<class Left, class Right, typename ScalarT>
//          // typename TensorQuery<Traits, Right::Order>::Type 
//         EXPR_TYPE(Traits, Right)
//          eval(const Binary< Binary<Number<ScalarT>, Left, Multiplies>, Right, Plus> &expr)
//         {
//             // typename TensorQuery<Traits, Right::Order>::Type result;
//             EXPR_TYPE(Traits, Right) result;

//             bool ok = zaxpy(expr.left().left(), eval(expr.left().right()), eval(expr.right()), result);
//             assert(ok);
//             return result;
//         }


//         template<class Left, typename ScalarT>
//         void eval(const Assign<Left, Binary< Number<ScalarT>, 
//                                              Factory<Identity, 2>, Multiplies> > &assign)
//         {
//             BackendT::Instance().build( assign.left().implementation(), 
//                                          size(assign.right().right()),
//                                          assign.right().right().type() );

//             BackendT::Instance().scal(assign.right().left(), assign.left().implementation(), assign.left().implementation());
            
//         }


//         template<class Left>
//         typename TypeAndFill<Traits, Left>::Type eval(const Binary<Left, Factory<Identity, 2>, Plus> &expr)
//         {
//             static_assert(Left::Order == 2, "can only be instantiated for 2nd order tensors");
//             typename TypeAndFill<Traits, Left>::Type result = eval(expr.left());
//             bool ok = BackendT::Instance().mat_diag_shift(result, 1.0);
//             assert(ok);
//             return result;
//         }

//         template<class Left, typename ScalarT>
//         // typename TensorQuery<Traits, 2>::Type 
//         EXPR_TYPE(Traits, Left)
//         eval(const Binary<
//                                 Left, 
//                                 Binary<  
//                                         Number<ScalarT>,
//                                         Factory<Identity, 2>,
//                                         Multiplies
//                                     >, 
//                                 Plus
//                                 > &expr)// -> decltype( eval(expr.left()) )
//         {
//             // typename TensorQuery<Traits, 2>::Type result = eval(expr.left());
//             EXPR_TYPE(Traits, Left) result = eval(expr.left());
//             bool ok = BackendT::Instance().mat_diag_shift(result, expr.right().left());
//             assert(ok);
//             return result;
//         }

//         template<class Left, class Right, typename ScalarT>
//         EXPR_TYPE(Traits, Right)
//         // Result 
//         eval(const Binary< Binary<Left, Number<ScalarT>, Multiplies>, Right, Plus > &expr)
//         {
//             // Result result;

//             EXPR_TYPE(Traits, Right) result;
//             bool ok = zaxpy(expr.left().right(),  eval(expr.left().left()),  eval(expr.right()), result);
//             assert(ok);
//             return result;
//         }

//         template<class Left, class Right, typename ScalarT>
//         // Result 
//         EXPR_TYPE(Traits, Right) 
//         eval(const Binary< Binary<Left, Number<ScalarT>, Multiplies>, Right, Minus > &expr)
//         {
//             // Result result;
//             EXPR_TYPE(Traits, Right) result;
//             bool ok =  zaxpy(-expr.left().right(),
//                         eval(expr.left().left()),
//                         eval(expr.right()));
//             assert(ok);
//             return result;
//         }

//         template<class ResultT>
//         bool zaxpy(const Scalar scaleFactor, const ResultT &left, const ResultT &right, ResultT &result)
//         {
//             return BackendT::Instance().zaxpy(scaleFactor, left, right, result);
//         }


//         ////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//         ////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//         ////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//         ////////////////////////////////////////////////////////////////////////////////////////////////////////////////

//         template<class Left, class Right, class Operation>
//         typename TypeAndFill<Traits, Binary<Left, Right, Operation> >::Type eval(const Binary<Left, Right, Operation> &expr)
//         {
//             typename TypeAndFill<Traits, Binary<Left, Right, Operation> >::Type result;
//             bool ok =  BackendT::Instance().apply( eval(expr.left()), eval(expr.right()), expr.operation(), result );
//             assert(ok);
//             return result;
//         }

//         template<class Left, class Right>
//         typename TypeAndFill<Traits, Multiply<Left, Right> >::Type eval(const Multiply<Left, Right> &expr)
//         {
//             typename TypeAndFill<Traits, Multiply<Left, Right> >::Type result;
//             bool ok = BackendT::Instance().apply( eval(expr.left()), eval(expr.right()), Multiplies(), result );
//             assert(ok);
//             return result;
//         }

//         template<class Left, class Right>
// 		typename TypeAndFill<Traits, Multiply<Transposed<Left>, Right> >::Type eval(const Multiply<Transposed<Left>, Right> &expr)
//         {
//             typename TypeAndFill<Traits, Multiply<Transposed<Left>, Right> >::Type result;
//             bool ok = BackendT::Instance().gemm(1.0, eval(expr.left().expr()), eval(expr.right()), true, false, 0.0, result);
//             assert(ok);
//             return result;
//         }

//         template<class Left, class Right>
// 		typename TypeAndFill<Traits, Multiply<Left, Transposed<Right> > >::Type eval(const Multiply<Left, Transposed<Right> > &expr)         {
//             typename TypeAndFill<Traits, Multiply<Left, Transposed<Right> > >::Type result;
//             bool ok = BackendT::Instance().gemm(1.0, eval(expr.left()), eval(expr.right().expr()), false, true, 0.0, result);
//             assert(ok);
//             return result;
//         }

//         template<class Left, class Right>
//         typename TypeAndFill<Traits, Multiply<Transposed<Left>, Transposed<Right>> >::Type  eval(const Multiply< Transposed<Left>, Transposed<Right> > &expr)
//         {
//             typename TypeAndFill<Traits, Multiply<Transposed<Left>, Transposed<Right> > >::Type result;
//             bool ok = BackendT::Instance().gemm(1.0, eval(expr.left().expr()), eval(expr.right().expr()), true, true, 0.0,
//                                                 result);
//             assert(ok);
//             return result;
//         }

//         template<class Left, class Right>
//         typename TypeAndFill<Traits, Transposed< Multiply<Left, Right> > >::Type  eval(const Transposed< Multiply<Left, Right> > &expr)
//         {
//             typename TypeAndFill<Traits, Transposed< Multiply<Left, Right> >  >::Type result;
//             bool ok = BackendT::Instance().gemm(1.0, eval(expr.expr().right()), eval(expr.expr().left()), true, true, 0.0,
//                                                 result);
//             assert(ok);
//             return result;
//         }

//         template<class Tensor>
//         typename TypeAndFill<Traits, Tensor>::Type
//         eval(const Transposed<Tensor> &t)
//         {
//             typename TypeAndFill<Traits, Tensor>::Type result;
//             bool ok = BackendT::Instance().transpose(eval(t.expr()), result);
//             assert(ok);
//             return result;
//         }

//         /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//         /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//         //////////////////////////////////////////////////OPERATIONS WITH SCALARS ///////////////////////////////////////////////////////

//         template<class Left, class ScalarT>
//         typename TypeAndFill<Traits, Binary<Left, Number<ScalarT>, Multiplies> >::Type eval(const Binary<Left, Number<ScalarT>, Multiplies> &expr)
//         {
//             typename TypeAndFill<Traits, Binary<Left, Number<ScalarT>, Multiplies> >::Type result;
//             bool ok = BackendT::Instance().apply(expr.right(), eval(expr.left()), expr.operation(), result);
//             assert(ok);
//             return result;
//         }

//         template<class ScalarT, class Right, class Operation>
//         typename TypeAndFill<Traits, Binary<Number<ScalarT>, Right, Operation> >::Type eval(const Binary<Number<ScalarT>, Right, Operation> &expr)
//         {
//             typename TypeAndFill<Traits, Binary<Number<ScalarT>, Right, Operation> >::Type result;
//             bool ok = BackendT::Instance().apply(expr.left(), eval(expr.right()), expr.operation(), result);
//             assert(ok);
//             return result;
//         }


//         template<typename L>
//         inline L eval(const Number<L> &lit) const
//         {
//             return lit;
//         }

//         template<class Implementation, int Order>
//         const Implementation &eval(const Wrapper<Implementation, Order> &expr)
//         {
//             return expr.implementation();
//         }

//         template<class Implementation, int Order>
//         Implementation &eval(Wrapper<Implementation, Order> &expr)
//         {
//             return expr.implementation();
//         }

//         template<class Implementation, int Order>
//         const Implementation &eval(const Wrapper<const Implementation &, Order> &expr)
//         {
//             return expr.implementation();
//         }

//         template<class Implementation, int Order>
//         const Implementation &eval(const Wrapper<Implementation &, Order> &expr)
//         {
//             return expr.implementation();
//         }

//         template <class Derived, class Operation>
//         Scalar eval(const Reduce<Derived, Operation> &expr)
//         {
//             return BackendT::Instance().reduce(eval(expr.expr()), expr.operation());
//         }

//         template<class Left, class Right>
//         Scalar eval(const Reduce<Binary<Left, Right, EMultiplies>, Plus> &expr)
//         {

//             //std::cerr << expr.getClass() << std::endl;
//             return BackendT::Instance().dot( eval(expr.expr().left()),
//                                              eval(expr.expr().right()));

//         }

//         template<class Expr>
//         Scalar eval(const Norm<Expr, 2> &expr)
//         {
//             return BackendT::Instance().norm2(eval(expr.expr()));
//         }

//         template<class Expr>
//         Scalar eval(const Norm<Expr, 1> &expr) {
//             return BackendT::Instance().norm1(eval(expr.expr()));

// //            //TODO
// //            assert(false);
// //            std::cerr << "Norm<Expr, 1> TODO in backend" << std::endl;
// //            return 0;
//         }

//         template<class Expr>
//         Scalar eval(const Norm<Expr, INFINITY_NORM_TAG> &expr)
//         {
//             return BackendT::Instance().norm_infty(eval(expr.expr()));
//         }


//         template<class Left, class Right, class Operation>
//         void eval(const InPlace<Left, Right, Operation> &expr)
//         {
//             //FIXME connect to backend without tree transformation
//             typedef utopia::Binary<Left, Right, Operation> ExprT;
//             eval(Assign<Left, ExprT>(expr.left(), ExprT(expr.left(), expr.right(), expr.operation())));
//         }




//         ////////////////////////////////////////////////// ////////////////////////////////////////////////////////////
//         ////////////////////////////////////////////////// ////////////////////////////////////////////////////////////
//         //////////////////////////////////////////////////BOOLEAN//////////////////////////////////////////////////////

//         template<class Expr>
//         bool eval(const Boolean<Expr> &expr) {
//             return eval(expr.expr());
//         }

//         template<class Left, class Right>
//         bool eval(const Reduce<Binary<Left, Right, ApproxEqual>, And> &expr) {
//             return BackendT::Instance().compare(eval(expr.expr().left()), eval(expr.expr().right()),
//                                                 expr.expr().operation());
//         }


//         template<class Expr>
//         typename TensorQuery<Traits, 1>::Type eval(const Unary< Wrapper<Expr, 1>, Reciprocal<Scalar> > &expr)
//         {
//             typename TensorQuery<Traits, 1>::Type result;
//             bool ok =  BackendT::Instance().apply( eval(expr.expr()), expr.operation(), result );
//             assert(ok);
//             return result;
//         }


//         template<class Expr, class Operation>
//         typename TypeAndFill<Traits, Unary<Expr, Operation> >::Type eval(const Unary<Expr, Operation> &expr)
//         {
//             typename TypeAndFill<Traits, Unary<Expr, Operation> >::Type result;
//             bool ok =  BackendT::Instance().apply( eval(expr.expr()), expr.operation(), result );
//             assert(ok);
//             return result;
//         }
        

        
//         template<class leftExpr, class RightExpr>
//         bool eval(const Construct<leftExpr, LocalDiagBlock<RightExpr> > & expr)
//         {
//             bool ok = BackendT::Instance().build_local_diag_block(eval(expr.left()), eval(expr.right().expr()));
//             assert(ok);
//             return ok;
//         }

//         template<class Left, class Right>
//         Result eval(const MatrixPtAPProduct<Left, Right> &expr)
//         {
//             Result result;
//             bool ok = BackendT::Instance().triple_product_PtAP(eval(expr.left()), eval(expr.right()), result);
//             return result;
//         }


// #ifdef WITH_PETSC
//         // local change of sizes 
//         template<class Left, class Right>
//         Result eval(const LocalRedistribute<Left, Right> &expr)
//         {
//             Result result;
//             bool ok = BackendT::Instance().build_local_redistribute(eval(expr.left()), eval(expr.right()), result);
//             return result;
//         }
// #endif //WITH_PETSC

//         template<class Expr>
//         Scalar eval(const Reduce< Diag<Expr>, Plus> &expr)
//         {
//             return BackendT::Instance().trace( eval(expr.expr().expr()) );
//         }

//         template<class Expr>
//         Scalar eval(const Trace<Expr> &expr)
//         {
//             return BackendT::Instance().trace( eval(expr.expr()) );
//         }
// #endif //USE_NEW_EVAL        

    };
}

#endif //utopia_utopia_EVALUATOR_HPP
