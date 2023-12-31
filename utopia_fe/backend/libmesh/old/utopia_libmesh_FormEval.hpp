#ifndef UTOPIA_LIBMESH_LINEAR_FORM_EVAL_HPP
#define UTOPIA_LIBMESH_LINEAR_FORM_EVAL_HPP

#include "utopia_Base.hpp"
#include "utopia_libmesh_FEForwardDeclarations.hpp"

#include "utopia_FEEval.hpp"
#include "utopia_FormEval.hpp"

#include "utopia_libmesh_AssemblyContext.hpp"
#include "utopia_libmesh_FunctionSpace.hpp"

namespace utopia {
    template <class Expr>
    class OffsetFinder {
    public:
        template <class Any>
        inline constexpr static int visit(const Any &) {
            return TRAVERSE_CONTINUE;
        }

        template <class T>
        inline int visit(const TestFunction<T> &expr) {
            // TODO
            return TRAVERSE_CONTINUE;
        }

        template <class T>
        inline int visit(const TestFunction<ProductFunctionSpace<T>> &expr) {
            // TODO
            return TRAVERSE_CONTINUE;
        }

        template <class T>
        inline int visit(const TrialFunction<T> &expr) {
            // TODO
            return TRAVERSE_CONTINUE;
        }

        template <class T>
        inline int visit(const TrialFunction<ProductFunctionSpace<T>> &expr) {
            // TODO
            return TRAVERSE_CONTINUE;
        }

        template <class ExprTree>
        inline void apply(const ExprTree &expr) {
            traverse(expr, *this);
        }

        OffsetFinder(AssemblyContext<LIBMESH_TAG> &ctx) : ctx(ctx), row_offset(0), col_offset(0) {}

        AssemblyContext<LIBMESH_TAG> &ctx;

        int row_offset;
        int col_offset;
    };

    template <class Form>
    class FormEval<Form, LIBMESH_TAG> {
    public:
        typedef utopia::Traits<LibMeshFunctionSpace> Traits;
        FormEval() {}

        template <class FS, class Tensor>
        static void apply(LinearIntegrator<FS> &expr, Tensor &t, AssemblyContext<LIBMESH_TAG> &ctx) {
            if (expr.is_surface()) {
                if (ctx.n_sides() != 0) {
                    ctx.surface_integral_begin();

                    for (std::size_t i = 0; i < ctx.n_sides(); ++i) {
                        ctx.set_side(i);

                        if (expr.assemble(ctx, t)) {
                            ctx.set_has_assembled(true);
                        }
                    }

                    ctx.surface_integral_end();
                }

            } else {
                if (expr.assemble(ctx, t)) {
                    ctx.set_has_assembled(true);
                }
            }
        }

        template <class FS, class Tensor>
        static void apply(BilinearIntegrator<FS> &expr, Tensor &t, AssemblyContext<LIBMESH_TAG> &ctx) {
            if (expr.is_surface()) {
                if (ctx.n_sides() != 0) {
                    ctx.surface_integral_begin();

                    for (std::size_t i = 0; i < ctx.n_sides(); ++i) {
                        ctx.set_side(i);

                        if (expr.assemble(ctx, t)) {
                            ctx.set_has_assembled(true);
                        }
                    }

                    ctx.surface_integral_end();
                }

            } else {
                if (expr.assemble(ctx, t)) {
                    ctx.set_has_assembled(true);
                }
            }
        }

        template <class FS, class Matrix, class Vector>
        static void apply(EquationIntegrator<FS> &expr, Matrix &mat, Vector &vec, AssemblyContext<LIBMESH_TAG> &ctx) {
            if (expr.is_surface()) {
                if (ctx.n_sides() != 0) {
                    ctx.surface_integral_begin();

                    for (std::size_t i = 0; i < ctx.n_sides(); ++i) {
                        ctx.set_side(i);

                        if (expr.assemble(ctx, mat, vec)) {
                            ctx.set_has_assembled(true);
                        }
                    }

                    ctx.surface_integral_end();
                }

            } else {
                if (expr.assemble(ctx, mat, vec)) {
                    ctx.set_has_assembled(true);
                }
            }
        }

        template <class Expr, class Tensor>
        static void apply(const Integral<Expr> &expr, Tensor &t, AssemblyContext<LIBMESH_TAG> &ctx) {
            if (expr.is_surface()) {
                if (ctx.n_sides() != 0) {
                    ctx.surface_integral_begin();

                    for (std::size_t i = 0; i < ctx.n_sides(); ++i) {
                        ctx.set_side(i);

                        if (expr.has_block_id() && ctx.block_id() != expr.block_id()) continue;

                        // t += FEEval<Integral<Expr>, Traits, LIBMESH_TAG, QUAD_DATA_NO>::apply(expr, ctx);
                        add(FEEval<Integral<Expr>, Traits, LIBMESH_TAG, QUAD_DATA_NO>::apply(expr, ctx), t);
                        ctx.set_has_assembled(true);
                    }

                    ctx.surface_integral_end();
                }

            } else if (expr.is_inner_volume_only()) {
                if (!ctx.on_boundary()) {
                    if (expr.has_block_id() && ctx.block_id() != expr.block_id()) {
                        return;
                    }

                    auto &&r = FEEval<Integral<Expr>, Traits, LIBMESH_TAG, QUAD_DATA_NO>::apply(expr, ctx);
                    assign(r, t);

                    ctx.set_has_assembled(true);
                }

            } else {
                assert(expr.is_volume());
                if (expr.has_block_id() && ctx.block_id() != expr.block_id()) {
                    return;
                }

                auto &&r = FEEval<Integral<Expr>, Traits, LIBMESH_TAG, QUAD_DATA_NO>::apply(expr, ctx);
                assign(r, t);

                ctx.set_has_assembled(true);
            }
        }

        template <typename T1, typename T2>
        inline static void add(const T1 &val, T2 &accumulator) {
            accumulator += val;
        }

        template <typename T1, typename T2>
        inline static void add(const std::vector<T1> &vals, Number<T2> &accumulator) {
            for (const auto &v : vals) {
                accumulator += v;
            }
        }

        template <typename T1, typename T2>
        inline static void add(const Number<T1> &val, Number<T2> &accumulator) {
            accumulator += val;
        }

        template <typename T1, typename T2>
        inline static void assign(const T1 &val, T2 &accumulator) {
            accumulator = val;
        }

        template <typename T>
        inline static void assign(const std::vector<T> &vals, Number<T> &accumulator) {
            accumulator = 0.;
            for (const auto &v : vals) {
                accumulator += v;
            }
        }

        template <typename T1, typename T2>
        inline static void assign(const Number<T1> &val, Number<T2> &accumulator) {
            accumulator = val;
        }

        template <class Left, class Right, class MatT, class VecT>
        static void apply(const Equality<Left, Right> &expr,
                          Tensor<MatT, 2> &mat,
                          Tensor<VecT, 1> &vec,
                          AssemblyContext<LIBMESH_TAG> &ctx) {
            apply(expr.left(), mat.derived(), ctx);
            apply(expr.right(), vec.derived(), ctx);
        }

        template <class Matrix, class Vector>
        class EquationAssembler {
        public:
            EquationAssembler(Matrix &mat, Vector &vec, AssemblyContext<LIBMESH_TAG> &ctx)
                : mat(mat), vec(vec), ctx(ctx) {}

            template <class Eq>
            void operator()(const int index, const Eq &eq) {
                FormEval<Form, LIBMESH_TAG>::apply(eq, mat_buff, vec_buff, ctx);
                if (empty(mat) || size(mat) != size(mat_buff)) {
                    mat = mat_buff;
                    vec = vec_buff;
                } else {
                    mat += mat_buff;
                    vec += vec_buff;
                }
            }

            Matrix &mat;
            Vector &vec;

            Matrix mat_buff;
            Vector vec_buff;
            AssemblyContext<LIBMESH_TAG> &ctx;
        };

        template <class... Eq, class MatT, class VecT>
        static void apply(const Equations<Eq...> &eqs,
                          Tensor<MatT, 2> &mat,
                          Tensor<VecT, 1> &vec,
                          AssemblyContext<LIBMESH_TAG> &ctx) {
            // FEBackend<LIBMESH_TAG>::init_context(eqs, ctx);
            EquationAssembler<MatT, VecT> eq_ass(mat.derived(), vec.derived(), ctx);
            eqs.each(eq_ass);
        }

        template <class Left, class Right, class Tensor>
        static void apply(const Binary<Left, Right, Plus> &expr, Tensor &result, AssemblyContext<LIBMESH_TAG> &ctx) {
            apply(expr.left(), result, ctx);

            Tensor right;
            right.zeros(layout(result));
            apply(expr.right(), right, ctx);
            result += right;
        }

        template <class Left, class Right, class Tensor>
        static void apply(const Binary<Number<Left>, Right, Multiplies> &expr,
                          Tensor &result,
                          AssemblyContext<LIBMESH_TAG> &ctx) {
            apply(expr.right(), result, ctx);
            result *= expr.left();
        }

        template <class Left, class Right, class Tensor>
        static void apply(const Binary<Left, Right, Minus> &expr, Tensor &result, AssemblyContext<LIBMESH_TAG> &ctx) {
            apply(expr.left(), result, ctx);

            Tensor right;
            right.zeros(layout(result));
            apply(expr.right(), right, ctx);
            result -= right;
        }

        template <class Left, class Right, class Op, class Tensor>
        static void apply(const Binary<Left, Right, Op> &expr, Tensor &result, AssemblyContext<LIBMESH_TAG> &ctx) {
            apply(expr.left(), result, ctx);

            Tensor right;
            right.zeros(layout(result));
            apply(expr.right(), result, ctx);

            result = expr.operation().apply(result, right);
        }

        template <class Expr, class Tensor>
        static void apply(const Unary<Expr, Minus> &expr, Tensor &result, AssemblyContext<LIBMESH_TAG> &ctx) {
            apply(expr.expr(), result, ctx);
            result = -result;
        }

        template <class Expr, class Tensor>
        static void apply(const Negate<Expr> &expr, Tensor &result, AssemblyContext<LIBMESH_TAG> &ctx) {
            apply(expr.expr(), result, ctx);
            result = -result;
        }

        template <class Expr, class Tensor>
        static void apply(const Unary<Expr, Abs> &expr, Tensor &result, AssemblyContext<LIBMESH_TAG> &ctx) {
            apply(expr.expr(), result, ctx);
            result = abs(result);
        }

        template <class Expr, class Tensor>
        static void apply(const Unary<Expr, Sqrt> &expr, Tensor &result, AssemblyContext<LIBMESH_TAG> &ctx) {
            apply(expr.expr(), result, ctx);
            result = sqrt(result);
        }
    };

    // template<class Derived, int Backend>
    // static auto eval(const Expression<Derived> &expr, AssemblyContext<Backend> &ctx)
    // -> decltype( FEEval<Derived, Traits<Derived>, Backend, QUAD_DATA_NO>::apply(expr.derived(), ctx) )
    // {
    // 	return FEEval<Derived, Traits<Derived>, Backend, QUAD_DATA_NO>::apply(expr.derived(), ctx);
    // }

    // template<class Derived, int Backend>
    // static auto quad_eval(const Expression<Derived> &expr, AssemblyContext<Backend> &ctx)
    // -> decltype( FEEval<Derived, Traits<Derived>, Backend, QUAD_DATA_YES>::apply(expr.derived(), ctx) )
    // {
    // 	return FEEval<Derived, Traits<Derived>, Backend, QUAD_DATA_YES>::apply(expr.derived(), ctx);
    // }

    template <class Derived>
    static auto eval(const Expression<Derived> &expr, AssemblyContext<LIBMESH_TAG> &ctx) -> decltype(
        FEEval<Derived, utopia::Traits<LibMeshFunctionSpace>, LIBMESH_TAG, QUAD_DATA_NO>::apply(expr.derived(), ctx)) {
        return FEEval<Derived, utopia::Traits<LibMeshFunctionSpace>, LIBMESH_TAG, QUAD_DATA_NO>::apply(expr.derived(),
                                                                                                       ctx);
    }

    template <class Derived>
    static auto eval(Expression<Derived> &expr, AssemblyContext<LIBMESH_TAG> &ctx) -> decltype(
        FEEval<Derived, utopia::Traits<LibMeshFunctionSpace>, LIBMESH_TAG, QUAD_DATA_NO>::apply(expr.derived(), ctx)) {
        return FEEval<Derived, utopia::Traits<LibMeshFunctionSpace>, LIBMESH_TAG, QUAD_DATA_NO>::apply(expr.derived(),
                                                                                                       ctx);
    }

    template <class Derived>
    static auto quad_eval(const Expression<Derived> &expr, AssemblyContext<LIBMESH_TAG> &ctx) -> decltype(
        FEEval<Derived, utopia::Traits<LibMeshFunctionSpace>, LIBMESH_TAG, QUAD_DATA_YES>::apply(expr.derived(), ctx)) {
        return FEEval<Derived, utopia::Traits<LibMeshFunctionSpace>, LIBMESH_TAG, QUAD_DATA_YES>::apply(expr.derived(),
                                                                                                        ctx);
    }
}  // namespace utopia

#endif  // UTOPIA_LIBMESH_LINEAR_FORM_EVAL_HPP
