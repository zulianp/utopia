#ifndef UTOPIA_PETSC_EVAL_DOT_OP_DOT_HPP
#define UTOPIA_PETSC_EVAL_DOT_OP_DOT_HPP

#include "utopia_Eval_Empty.hpp"
#include "utopia_ForwardDeclarations.hpp"

//Divides
// |	Reduce
// |	|	EMultiplies
// |	|	|	Vector
// |	|	|	Vector
// |	|	Plus
// |	Reduce
// |	|	EMultiplies
// |	|	|	Multiply
// |	|	|	|	SparseMatrix
// |	|	|	|	Vector
// |	|	|	Vector
// |	|	Plus

namespace utopia {
    template<class X, class FunOfX, class Op, class Traits>
    class Eval< Binary<
                    Dot<Tensor<X, 1>, Tensor<X, 1> >,
                    Dot<FunOfX, Tensor<X, 1> >,
                    Op>, Traits, PETSC> {
    public:
        typedef Binary< Dot<Tensor<X, 1>, Tensor<X, 1> >,
                        Dot<FunOfX, Tensor<X, 1> >,
                        Op> Expr;

        using Scalar = typename Traits::Scalar;

        template<typename T>
        inline static void apply(const Expr &expr, Number<T> &num)
        {
            num = apply(expr);
        }


        inline static Scalar apply(const Expr &expr) {
            UTOPIA_TRACE_BEGIN(expr);

            const auto &x1 = expr.left().expr().left().derived();
            const auto &x2 = expr.left().expr().right().derived();

            auto &&x3 = Eval<FunOfX, Traits>::apply(expr.right().expr().left());
            const auto &x4 = expr.right().expr().right().derived();


            PetscScalar left_num = 0., right_num = 0.;

            PetscErrorCode ierr = 0;
            if(x1.raw_type() == x4.raw_type()) {
                Vec vecs[2] = { x2.raw_type(), x3.raw_type() };
                PetscScalar vals[2] = { 0., 0.};
                ierr = VecMDot(x1.raw_type(), 2, vecs, vals); assert(ierr == 0);
                left_num = vals[0];
                right_num = vals[1];
            } else if(x2.raw_type() == x4.raw_type()) {
                Vec vecs[2] = { x1.raw_type(), x3.raw_type() };
                PetscScalar vals[2] = { 0., 0.};
                ierr = VecMDot(x2.raw_type(), 2, vecs, vals); assert(ierr == 0);
                left_num = vals[0];
                right_num = vals[1];
            } else {
                ierr = VecDotBegin(x1.raw_type(), x2.raw_type(), &left_num); assert(ierr == 0);
                ierr = VecDotBegin(x3.raw_type(), x4.raw_type(), &right_num); assert(ierr == 0);

                ierr = VecDotEnd(x1.raw_type(), x2.raw_type(), &left_num); assert(ierr == 0);
                ierr = VecDotEnd(x3.raw_type(), x4.raw_type(), &right_num); assert(ierr == 0);
            }

            Scalar r = 0;
            if(std::is_same<Op, Divides>::value) {
                r = left_num/right_num;
                if(right_num == 0. || ierr != 0) {
                    r = 0.;
                }
            } else {
                r = expr.operation().apply(left_num, right_num);
            }

            UTOPIA_TRACE_END(expr);
            return r;
        }
    };


    // used for computing predicted reduction in TR methods
    // result = -num1 * dot(x1, x2) - num2 * dot(B * x1, x2)
    template<class Num, class X, class FunOfX, class Op, class Traits>
    class Eval< Binary<
                     Binary<Number<Num>,  Dot<Tensor<X, 1>, Tensor<X, 1> >, Multiplies>,
                    Binary<Number<Num>,  Dot<FunOfX, Tensor<X, 1> >, Multiplies>,
                    Op>, Traits, PETSC>

    {
        public:
            typedef Binary<Binary<Number<Num>,  Dot<Tensor<X, 1>, Tensor<X, 1> >, Multiplies>,Binary<Number<Num>,  Dot<FunOfX, Tensor<X, 1> >, Multiplies>,Op> Expr;

            using Scalar = typename Traits::Scalar;

            template<typename T>
            inline static void apply(const Expr &expr, Number<T> &num)
            {
                num = apply(expr);
            }

            inline static Scalar apply(const Expr &expr) {
                UTOPIA_TRACE_BEGIN(expr);

                const auto & multiplierLeft  = expr.left().left();
                const auto & multiplierRight  = expr.right().left();

                const auto &x1 = expr.left().right().expr().left().derived();
                const auto &x2 = expr.left().right().expr().right().derived();

                auto &&x3 = Eval<FunOfX, Traits>::apply(expr.right().right().expr().left());
                const auto &x4 = expr.right().right().expr().right().derived();


                PetscScalar left_num = 0., right_num = 0.;

                PetscErrorCode ierr = 0;
                if(x1.raw_type() == x4.raw_type()) {
                    Vec vecs[2] = { x2.raw_type(), x3.raw_type() };
                    PetscScalar vals[2] = { 0., 0.};
                    ierr = VecMDot(x1.raw_type(), 2, vecs, vals); assert(ierr == 0);
                    left_num = vals[0];
                    right_num = vals[1];
                } else if(x2.raw_type() == x4.raw_type()) {
                    Vec vecs[2] = { x1.raw_type(), x3.raw_type() };
                    PetscScalar vals[2] = { 0., 0.};
                    ierr = VecMDot(x2.raw_type(), 2, vecs, vals); assert(ierr == 0);
                    left_num = vals[0];
                    right_num = vals[1];
                } else {
                    ierr = VecDotBegin(x1.raw_type(), x2.raw_type(), &left_num); assert(ierr == 0);
                    ierr = VecDotBegin(x3.raw_type(), x4.raw_type(), &right_num); assert(ierr == 0);

                    ierr = VecDotEnd(x1.raw_type(), x2.raw_type(), &left_num); assert(ierr == 0);
                    ierr = VecDotEnd(x3.raw_type(), x4.raw_type(), &right_num); assert(ierr == 0);
                }

                left_num *= multiplierLeft;
                right_num *= multiplierRight;

                Scalar r = 0;
                if(std::is_same<Op, Divides>::value) {
                    r = left_num/right_num;
                    if(right_num == 0. || ierr != 0) {
                        r = 0.;
                    }
                } else {
                    r = expr.operation().apply(left_num, right_num);
                }

                UTOPIA_TRACE_END(expr);
                return r;
        }
    };


}


#endif //UTOPIA_PETSC_EVAL_DOT_OP_DOT_HPP
