#ifndef UTOPIA_LIBMESH_FE_BACKEND_HPP
#define UTOPIA_LIBMESH_FE_BACKEND_HPP

#include "utopia_Reduce.hpp"
#include "utopia_libmesh_AssemblyContext.hpp"
#include "utopia_libmesh_FEForwardDeclarations.hpp"
#include "utopia_libmesh_FunctionSpace.hpp"
#include "utopia_libmesh_LambdaFunction.hpp"

#include "utopia_FEBackend.hpp"
#include "utopia_FEBasicOverloads.hpp"

#include <array>
#include <functional>
#include <type_traits>
#include "libmesh/const_function.h"

namespace utopia {

    // template<typename T>
    // void disp(const std::vector<std::vector<T>> &data, std::ostream &os = std::cout) {
    //     for(std::size_t i = 0; i < data.size(); ++i) {
    //         for(std::size_t qp = 0; qp < data[i].size(); ++qp) {
    //             std::cout << "[" << i << ", " << qp << "]:\n";
    //             disp(data[i][qp], os);
    //         }
    //     }
    // }

    // template<typename T>
    // void disp(const std::vector<T> &data, std::ostream &os = std::cout) {
    //     for(std::size_t qp = 0; qp < data.size(); ++qp) {
    //         std::cout << "[" << qp << "]:\n";
    //         disp(data[qp], os);
    //     }
    // }

    inline static double inner(const double &left, const LMDenseVector &right) { return right.sum() * left; }

    inline static double inner(const LMDenseVector &right, const double &left) { return right.sum() * left; }

    template <class Left, typename T>
    inline static T inner(const Tensor<Left, 2> &left, const libMesh::VectorValue<T> &right) {
        return inner(left.derived(), right);
    }

    template <class Left, typename T>
    inline static T inner(const Tensor<Left, 1> &left, const libMesh::VectorValue<T> &right) {
        return inner(left.derived(), right);
    }

    template <class Left, typename T>
    inline static T inner(const Tensor<Left, 2> &left, const LMDenseMatrix &right) {
        return inner(left, right);
    }

    template <typename T, class Right>
    inline static T inner(const LMDenseMatrix &left, const Tensor<Right, 2> &right) {
        return inner(right, left);
    }

    ////////////////////////////////////////////////////////////////////////////////////////////

    template <>
    class FEBackend<LIBMESH_TAG> {
    public:
        typedef utopia::Traits<LibMeshFunctionSpace> TraitsT;

        static const int Backend = LIBMESH_TAG;
        static const int Order = 1;
        static const int FILL_TYPE = FillType::DENSE;

        typedef TraitsT::Scalar Scalar;
        typedef TraitsT::Vector Vector;
        typedef TraitsT::Matrix Matrix;
        typedef TraitsT::FE FE;
        typedef TraitsT::FunctionType FunctionType;
        typedef TraitsT::GradientType GradientType;
        typedef TraitsT::DivergenceType DivergenceType;
        typedef TraitsT::JacobianType JacobianType;
        typedef TraitsT::CurlType CurlType;
        typedef TraitsT::DXType DXType;

        template <typename T>
        using FQValues = std::vector<std::vector<T>>;

        template <typename T>
        using QValues = std::vector<T>;

        typedef FQValues<LMDenseVector> VectorFunctionType;

        //////////////////////////////////////////////////////////////////////////////////////////

        // system of equations
        template <class... Eq>
        static void init_context(const Equations<Eq...> &eqs, AssemblyContext<LIBMESH_TAG> &ctx) {}

        // single equations
        template <class Left, class Right>
        static void init_context(const Equality<Left, Right> &eq, AssemblyContext<LIBMESH_TAG> &ctx) {}

        // general expressions that might contain multilinear forms
        template <class Expr>
        static void init_context(const Expression<Expr> &expr, AssemblyContext<LIBMESH_TAG> &ctx) {}

        template <class Derived, int Order>
        static void init_tensor(Tensor<Derived, Order> &t, AssemblyContext<LIBMESH_TAG> &ctx) {
            ctx.init_tensor(t.derived(), true);
        }

        class ConstraintsInitializer {
        public:
            template <class Constr>
            inline void operator()(const int, const Constr &constr) const {
                std::cout << "unimplemented constraint expression" << std::endl;
                assert(false);
            }

            template <typename T, std::size_t N>
            static STL2LibMeshLambdaFunction<T, N> make_lambda(
                std::function<std::array<T, N>(const std::array<T, N> &)> f) {
                return STL2LibMeshLambdaFunction<T, N>(f);
            }

            template <typename T, std::size_t N>
            static STL2LibMeshLambdaFunction<T, N> make_lambda(std::function<T(const std::array<T, N> &)> f) {
                return STL2LibMeshLambdaFunction<T, N>(f);
            }

            static void make_vars(const TrialFunction<LibMeshFunctionSpace> &f, std::vector<unsigned int> &vars) {
                vars.resize(1);
                vars[0] = f.space_ptr()->subspace_id();
            }

            static void make_vars(const TrialFunction<ProductFunctionSpace<LibMeshFunctionSpace>> &f,
                                  std::vector<unsigned int> &vars) {
                vars.resize(f.space_ptr()->n_subspaces());

                f.space_ptr()->each(
                    [&vars](const int i, const LibMeshFunctionSpace &ss) { vars[i] = ss.subspace_id(); });
            }

            static libMesh::DofMap &get_dof_map(LibMeshFunctionSpace &s) { return s.dof_map(); }

            static libMesh::DofMap &get_dof_map(ProductFunctionSpace<LibMeshFunctionSpace> &s) {
                return s[0].dof_map();
            }

            template <class F, class S, typename T, int Order>
            inline void operator()(
                const int,
                const DirichletBoundaryCondition<Equality<TrialFunction<S>, FunctionCoefficient<F, T, Order>>> &cond)
                const {
                std::vector<unsigned int> vars;
                make_vars(cond.expr().left(), vars);

                auto lambda = make_lambda(cond.expr().right().fun());
                // lambda.set_is_time_dependent(true);

                std::set<libMesh::boundary_id_type> bt;
                bt.insert(cond.boundary_tags().begin(), cond.boundary_tags().end());
                libMesh::DirichletBoundary d_bc(bt, vars, lambda);
                get_dof_map(*cond.expr().left().space_ptr()).add_dirichlet_boundary(d_bc);
            }

            template <class F, typename T>
            inline void operator()(
                const int,
                const DirichletBoundaryCondition<
                    Equality<TrialFunction<LibMeshFunctionSpace>, FunctionCoefficient<F, T, 0>>> &cond) const {
                std::vector<unsigned int> vars(1);
                vars[0] = cond.expr().left().space_ptr()->subspace_id();
                std::set<libMesh::boundary_id_type> bt;
                bt.insert(cond.boundary_tags().begin(), cond.boundary_tags().end());
                libMesh::DirichletBoundary d_bc(
                    bt, vars, LibMeshLambdaFunction<libMesh::Real>(cond.expr().right().fun()));
                cond.expr().left().space_ptr()->dof_map().add_dirichlet_boundary(d_bc);
            }

            template <class T>
            inline void operator()(const int,
                                   const DirichletBoundaryCondition<Equality<TrialFunction<LibMeshFunctionSpace>,
                                                                             ConstantCoefficient<T, 0>>> &&cond) {
                std::vector<unsigned int> vars(1);
                vars[0] = cond.expr().left().space_ptr()->subspace_id();
                std::set<libMesh::boundary_id_type> bt;
                bt.insert(cond.boundary_tags().begin(), cond.boundary_tags().end());
                libMesh::DirichletBoundary d_bc(bt, vars, libMesh::ConstFunction<libMesh::Real>(cond.expr().right()));
                cond.expr().left().space_ptr()->dof_map().add_dirichlet_boundary(d_bc);
            }

#ifdef UTOPIA_WITH_TINY_EXPR
            template <class S>
            inline void operator()(
                const int,
                const DirichletBoundaryCondition<Equality<TrialFunction<S>, SymbolicFunction>> &cond) const {
                auto &f = cond.expr().right();

                std::function<libMesh::Real(const libMesh::Point &)> fun =
                    [f](const libMesh::Point &xyz) -> libMesh::Real {
                    // FIXME
                    auto f_copy = f;
                    return f_copy.eval(xyz(0), xyz(1), xyz(2));
                };

                LibMeshLambdaFunction<libMesh::Real> lambda(fun);
                // lambda.set_is_time_dependent(true);

                std::vector<unsigned int> vars(1);
                vars[0] = cond.expr().left().space_ptr()->subspace_id();
                std::set<libMesh::boundary_id_type> bt;
                bt.insert(cond.boundary_tags().begin(), cond.boundary_tags().end());
                libMesh::DirichletBoundary d_bc(bt, vars, lambda);
                cond.expr().left().space_ptr()->dof_map().add_dirichlet_boundary(d_bc);
            }
#endif  // UTOPIA_WITH_TINY_EXPR

            template <class T>
            void strong_enforce(
                DirichletBoundaryCondition<Equality<TrialFunction<ProductFunctionSpace<LibMeshFunctionSpace>>,
                                                    ConstantCoefficient<T, 1>>> &&cond) {
                std::cout << "unimplemented constraint exprassion" << std::endl;
                assert(false);

                // std::vector<unsigned int> vars(1); vars[0] = cond.expr().left().space_ptr()->subspace_id();
                // std::set<libMesh::boundary_id_type> bt;
                // bt.insert(cond.boundary_tags().begin(), cond.boundary_tags().end());
                // libMesh::DirichletBoundary d_bc(bt, vars, libMesh::ConstFunction<T>(cond.expr().right()) );
                // cond.expr().left().space_ptr()->dof_map().add_dirichlet_boundary(d_bc);
            }
        };

        template <class... Constr>
        static void init_constraints(const FEConstraints<Constr...> &constr) {
            ConstraintsInitializer c_init;
            constr.each(c_init);
            /*
                template<class Right, typename T>
                void strong_enforce(const DirichletBoundaryCondition<
                                    Equality<LibMeshFEFunction,
                                    FunctionCoefficient<Right, T, 0> >
                                    > &cond)
                {

                }

            */
        }

        template <typename T>
        inline static const T &get(const std::vector<std::vector<T>> &v, const std::size_t qp, const std::size_t i) {
            return v[i][qp];
        }

        template <typename T>
        inline static const T &get(const std::vector<T> &v, const std::size_t qp, const std::size_t) {
            return v[qp];
        }

        template <typename T, int Order>
        inline static const T get(const ConstantCoefficient<T, Order> &c, const std::size_t qp, const std::size_t i) {
            return c[i];
        }

        template <typename T, int Order>
        inline static T get(Tensor<T, Order> &&c, const std::size_t qp, const std::size_t i) {
            return std::move(c.derived());
        }

        template <typename T, int Order>
        inline static const T &get(const Tensor<T, Order> &c, const std::size_t qp, const std::size_t i) {
            return c.derived();
        }

        inline static const double get(const double &c, const std::size_t qp, const std::size_t i) { return c; }

        template <typename T>
        inline static void add(Matrix &mat, const int i, const int j, const T value) {
            mat.add(i, j, value);
        }

        template <typename T>
        inline static void add(Vector &vec, const int i, const int j, const T value) {
            assert(j == 0);
            vec.add(i, value);
        }

        //////////////////////////////////////////////////////////////////////////////////////////
        static inline unsigned int offset(const LibMeshFunctionSpace &space, const AssemblyContext<LIBMESH_TAG> &ctx) {
            return ctx.offset[space.subspace_id()];
        }

        static inline unsigned int offset(const ProductFunctionSpace<LibMeshFunctionSpace> &space,
                                          const AssemblyContext<LIBMESH_TAG> &ctx) {
            return offset(space[0], ctx);
        }

        static inline unsigned int n_shape_functions(const LibMeshFunctionSpace &space,
                                                     const AssemblyContext<LIBMESH_TAG> &ctx) {
            return ctx.fe()[space.subspace_id()]->n_shape_functions();
        }

        static inline Range range(const LibMeshFunctionSpace &space, const AssemblyContext<LIBMESH_TAG> &ctx) {
            auto off = offset(space, ctx);
            return {off, off + n_shape_functions(space, ctx)};
        }

        static inline Range range(const ProductFunctionSpace<LibMeshFunctionSpace> &space,
                                  const AssemblyContext<LIBMESH_TAG> &ctx) {
            auto off = offset(space, ctx);
            return {off, off + n_shape_functions(space, ctx)};
        }

        static inline unsigned int n_shape_functions(const ProductFunctionSpace<LibMeshFunctionSpace> &space,
                                                     const AssemblyContext<LIBMESH_TAG> &ctx) {
            unsigned int ret = 0;

            space.each([&](const int, const LibMeshFunctionSpace &subspace) {
                ret += ctx.fe()[subspace.subspace_id()]->n_shape_functions();
            });

            return ret;
        }
        //////////////////////////////////////////////////////////////////////////////////////////

        // Function
        // scalar fe functions
        static const FunctionType &fun(const TrialFunction<LibMeshFunctionSpace> &fun,
                                       AssemblyContext<LIBMESH_TAG> &ctx) {
            assert(!ctx.trial().empty());
            return ctx.trial()[fun.space_ptr()->subspace_id()]->get_phi();
        }

        static const FunctionType &fun(const TestFunction<LibMeshFunctionSpace> &fun,
                                       AssemblyContext<LIBMESH_TAG> &ctx) {
            assert(!ctx.test().empty());
            return ctx.test()[fun.space_ptr()->subspace_id()]->get_phi();
        }

        template <typename T>
        static T apply(const BlockVar<T> &var, const AssemblyContext<LIBMESH_TAG> &ctx) {
            return var.get(ctx.block_id());
        }

#ifdef UTOPIA_WITH_TINY_EXPR
        static QValues<double> apply(const SymbolicFunction &fun, const AssemblyContext<LIBMESH_TAG> &ctx) {
            const auto &xyz = ctx.test()[0]->get_xyz();
            auto n_quad_points = xyz.size();
            QValues<double> ret;
            ret.resize(n_quad_points);

            // FIXME
            auto f_copy = fun;

            for (std::size_t qp = 0; qp < n_quad_points; ++qp) {
                ret[qp] = f_copy.eval(xyz[qp](0), xyz[qp](1), xyz[qp](2));
            }

            return ret;
        }
#endif  // UTOPIA_WITH_TINY_EXPR

        static auto determinant(const QValues<LMDenseMatrix> &mats, const AssemblyContext<LIBMESH_TAG> &ctx)
            -> QValues<double> {
            const auto n = mats.size();
            std::vector<double> dets(n);
            for (std::size_t i = 0; i < n; ++i) {
                dets[i] = utopia::det(mats[i]);
            }

            return dets;
        }

        static Matrix build(const Size &s, const Identity &, const AssemblyContext<LIBMESH_TAG> &) {
            return identity(s);
        }

        static auto inverse(const QValues<LMDenseMatrix> &mats, const AssemblyContext<LIBMESH_TAG> &)
            -> QValues<LMDenseMatrix> {
            const auto n = mats.size();
            QValues<LMDenseMatrix> ret(n);
            for (std::size_t i = 0; i < n; ++i) {
                ret[i] = utopia::inv(mats[i]);
            }

            return ret;
        }

        // template<typename Tensor>
        static auto trace(const QValues<LMDenseMatrix> &mats, const AssemblyContext<LIBMESH_TAG> &ctx)
            -> QValues<double> {
            auto n_quad_points = mats.size();
            QValues<double> ret(n_quad_points);
            for (std::size_t qp = 0; qp < n_quad_points; ++qp) {
                ret[qp] = utopia::trace(mats[qp]);
            }

            return ret;
        }

        static auto trace_times_identity(const QValues<LMDenseMatrix> &mats, const AssemblyContext<LIBMESH_TAG> &ctx)
            -> QValues<LMDenseMatrix> {
            auto n_quad_points = mats.size();
            QValues<LMDenseMatrix> ret(mats.size());
            QValues<double> traces = trace(mats, ctx);

            for (std::size_t qp = 0; qp < n_quad_points; ++qp) {
                auto s = size(mats[qp]);
                ret[qp].identity(serial_layout(s.get(0), s.get(1)), traces[qp]);
                // ret[qp] = traces[qp] * identity(size(mats[qp]));
            }

            return ret;
        }

        template <typename Tensor>
        static auto trace_times_identity(const FQValues<Tensor> &mats, const AssemblyContext<LIBMESH_TAG> &ctx)
            -> FQValues<Tensor> {
            auto n_funs = mats.size();
            FQValues<Tensor> ret(n_funs);
            for (std::size_t i = 0; i < n_funs; ++i) {
                ret[i] = trace_times_identity(mats[i], ctx);
            }

            return ret;
        }

        static auto trace(const LMDenseMatrix &mat, const AssemblyContext<LIBMESH_TAG> &ctx) -> double {
            auto s = size(mat);
            auto n = s.get(0) < s.get(1) ? s.get(0) : s.get(1);

            Read<LMDenseMatrix> r_m(mat);

            auto ret = mat.get(0, 0);
            for (std::size_t i = 1; i < n; ++i) {
                ret += mat.get(i, i);
            }

            return ret;
        }

        static auto transpose(const LMDenseMatrix &mat, const AssemblyContext<LIBMESH_TAG> &ctx) -> LMDenseMatrix {
            // forward to algebra backend
            return utopia::transpose(mat);
        }

        static auto transpose(const QValues<LMDenseMatrix> &mats, const AssemblyContext<LIBMESH_TAG> &ctx)
            -> QValues<LMDenseMatrix> {
            const auto n = mats.size();
            QValues<LMDenseMatrix> ret(n);
            for (std::size_t i = 0; i < n; ++i) {
                ret[i] = utopia::transpose(mats[i]);
            }

            return ret;
        }

        static auto transpose(const FQValues<LMDenseMatrix> &mats, const AssemblyContext<LIBMESH_TAG> &ctx)
            -> FQValues<LMDenseMatrix> {
            const auto n = mats.size();
            const auto n_qp = mats[0].size();

            FQValues<LMDenseMatrix> ret(n);
            for (std::size_t i = 0; i < n; ++i) {
                ret[i].resize(n_qp);
                for (std::size_t k = 0; k < n_qp; ++k) {
                    ret[i][k] = utopia::transpose(mats[i][k]);
                }
            }

            return ret;
        }

        static auto transpose(FQValues<LMDenseMatrix> &&mats, const AssemblyContext<LIBMESH_TAG> &ctx)
            -> FQValues<LMDenseMatrix> {
            const auto n = mats.size();
            const auto n_qp = mats[0].size();

            for (std::size_t i = 0; i < n; ++i) {
                for (std::size_t k = 0; k < n_qp; ++k) {
                    mats[i][k] = utopia::transpose(mats[i][k]);
                }
            }

            return std::move(mats);
        }

        static auto apply_binary(const SymbolicTensor<Identity, 2> &,
                                 QValues<LMDenseMatrix> &&mats,
                                 const Plus &,
                                 const AssemblyContext<LIBMESH_TAG> &) -> QValues<LMDenseMatrix> {
            auto s = size(mats[0]);
            for (auto &m : mats) {
                // m += identity(s);
                m.shift_diag(1.0);
            }

            // std::cout << "-----------------------" << std::endl;
            // disp(mats);

            return std::move(mats);
        }

        static auto apply_binary(const SymbolicTensor<Identity, 2> &,
                                 const QValues<LMDenseMatrix> &mats,
                                 const Plus &,
                                 const AssemblyContext<LIBMESH_TAG> &) -> QValues<LMDenseMatrix> {
            QValues<LMDenseMatrix> ret = mats;
            auto s = size(mats[0]);
            for (auto &m : ret) {
                // m += identity(s);
                m.shift_diag(1.0);
            }

            // std::cout << "-----------------------" << std::endl;
            // disp(mats);

            return ret;
        }

        static auto apply_binary(const double factor,
                                 QValues<double> &&vals,
                                 const Multiplies &,
                                 const AssemblyContext<LIBMESH_TAG> &) -> QValues<double> {
            for (auto &v : vals) {
                v *= factor;
            }

            return std::move(vals);
        }

        template <typename T>
        static auto apply_binary(const ConstantCoefficient<T, 0> &factor,
                                 QValues<double> &&vals,
                                 const Multiplies &,
                                 const AssemblyContext<LIBMESH_TAG> &) -> QValues<double> {
            for (auto &v : vals) {
                v *= factor.expr();
            }

            return std::move(vals);
        }

        static auto apply_binary(QValues<double> &&left,
                                 const QValues<double> &right,
                                 const Multiplies &,
                                 const AssemblyContext<LIBMESH_TAG> &) -> QValues<double> {
            std::size_t n = left.size();

            for (std::size_t i = 0; i < n; ++i) {
                left[i] *= right[i];
            }

            return std::move(left);
        }

        static auto apply_binary(const double factor,
                                 const QValues<double> &vals_in,
                                 const Multiplies &,
                                 const AssemblyContext<LIBMESH_TAG> &) -> QValues<double> {
            QValues<double> vals;
            for (auto &v : vals) {
                v *= factor;
            }

            return vals;
        }

        template <class T, int Order>
        static auto apply_binary(const double factor,
                                 Tensor<T, Order> &&v,
                                 const Multiplies &,
                                 const AssemblyContext<LIBMESH_TAG> &) -> T {
            v.derived() *= factor;
            return std::move(v.derived());
        }

        static auto apply_binary(const LMDenseMatrix &mat,
                                 QValues<LMDenseMatrix> &&mats,
                                 const Minus &,
                                 const AssemblyContext<LIBMESH_TAG> &) -> QValues<LMDenseMatrix> {
            auto s = size(mats[0]);
            for (auto &m : mats) {
                m = mat - m;
            }

            return std::move(mats);
        }

        static auto apply_binary(const LMDenseMatrix &mat,
                                 QValues<LMDenseMatrix> &&mats,
                                 const Plus &,
                                 const AssemblyContext<LIBMESH_TAG> &) -> QValues<LMDenseMatrix> {
            auto s = size(mats[0]);
            for (auto &m : mats) {
                m += mat;
            }

            return std::move(mats);
        }

        static auto apply_binary(QValues<LMDenseMatrix> &&mats,
                                 const SymbolicTensor<Identity, 2> &,
                                 const Minus &,
                                 const AssemblyContext<LIBMESH_TAG> &) -> QValues<LMDenseMatrix> {
            auto s = size(mats[0]);
            for (auto &m : mats) {
                // m -= identity(s);
                m.shift_diag(-1.0);
            }

            return std::move(mats);
        }

        static auto apply_binary(const double val,
                                 LMDenseMatrix &&mat,
                                 const Multiplies &,
                                 const AssemblyContext<LIBMESH_TAG> &) -> LMDenseMatrix {
            mat *= val;
            return std::move(mat);
        }

        static auto apply_binary(const double val,
                                 QValues<LMDenseMatrix> &&mats,
                                 const Plus &,
                                 const AssemblyContext<LIBMESH_TAG> &) -> QValues<Matrix> {
            for (auto &m : mats) {
                m *= val;
            }

            return std::move(mats);
        }

        template <typename T>
        static auto apply_binary(QValues<T> &&vals,
                                 const double val,
                                 const Multiplies &op,
                                 const AssemblyContext<LIBMESH_TAG> &ctx) -> QValues<T> {
            return apply_binary(val, std::forward<QValues<T>>(vals), op, ctx);
        }

        template <typename T>
        static auto apply_binary(const double val,
                                 QValues<T> &&vals,
                                 const Multiplies &,
                                 const AssemblyContext<LIBMESH_TAG> &) -> QValues<T> {
            for (auto &v : vals) {
                v *= val;
            }

            return std::move(vals);
        }

        template <typename T>
        static auto apply_binary(const double val,
                                 FQValues<T> &&vals,
                                 const Multiplies &,
                                 const AssemblyContext<LIBMESH_TAG> &) -> FQValues<T> {
            for (auto &v : vals) {
                for (auto &vk : v) {
                    vk *= val;
                }
            }

            return std::move(vals);
        }

        template <typename T>
        static auto apply_binary(FQValues<T> &&left,
                                 const FQValues<T> &right,
                                 const Divides &op,
                                 const AssemblyContext<LIBMESH_TAG> &ctx) -> FQValues<T> {
            std::size_t n = left.size();
            std::size_t n_quad_points = left[0].size();

            for (std::size_t i = 0; i != n; ++i) {
                for (std::size_t qp = 0; qp < n_quad_points; ++qp) {
                    left[i][qp] /= right[i][qp];
                }
            }

            return std::move(left);
        }

        template <typename T>
        static auto apply_binary(QValues<T> &&left,
                                 const QValues<T> &right,
                                 const Divides &op,
                                 const AssemblyContext<LIBMESH_TAG> &ctx) -> QValues<T> {
            std::size_t n = left.size();

            for (std::size_t i = 0; i != n; ++i) {
                left[i] /= right[i];
            }

            return std::move(left);
        }

        template <typename T>
        static auto apply_binary(const QValues<T> &left,
                                 const QValues<T> &right,
                                 const Divides &op,
                                 const AssemblyContext<LIBMESH_TAG> &ctx) -> QValues<T> {
            std::size_t n = left.size();
            QValues<T> ret = left;

            for (std::size_t i = 0; i != n; ++i) {
                ret[i] /= right[i];
            }

            return std::move(ret);
        }

        template <typename T>
        static auto apply_binary(QValues<T> &&left,
                                 const QValues<T> &right,
                                 const Plus &,
                                 const AssemblyContext<LIBMESH_TAG> &) -> QValues<T> {
            std::size_t n = left.size();

            for (std::size_t i = 0; i != n; ++i) {
                left[i] += right[i];
            }

            return std::move(left);
        }

        template <typename T>
        static auto apply_binary(FQValues<T> &&left,
                                 const FQValues<T> &right,
                                 const Plus &,
                                 const AssemblyContext<LIBMESH_TAG> &) -> FQValues<T> {
            std::size_t n = left.size();
            std::size_t n_qp = left[0].size();

            for (std::size_t i = 0; i != n; ++i) {
                for (std::size_t j = 0; j != n_qp; ++j) {
                    left[i][j] += right[i][j];
                }
            }

            return std::move(left);
        }

        template <typename T>
        static auto apply_binary(QValues<T> &&left,
                                 const QValues<T> &right,
                                 const Minus &,
                                 const AssemblyContext<LIBMESH_TAG> &) -> QValues<T> {
            std::size_t n = left.size();

            for (std::size_t i = 0; i != n; ++i) {
                left[i] -= right[i];
            }

            return std::move(left);
        }

        template <typename T>
        static auto apply_binary(FQValues<T> &&left,
                                 const FQValues<T> &right,
                                 const Minus &,
                                 const AssemblyContext<LIBMESH_TAG> &) -> FQValues<T> {
            std::size_t n = left.size();
            std::size_t n_qp = left[0].size();

            for (std::size_t i = 0; i != n; ++i) {
                for (std::size_t j = 0; j != n_qp; ++j) {
                    left[i][j] -= right[i][j];
                }
            }

            return std::move(left);
        }

        template <typename T1, typename T2, class Op>
        static auto apply_binary(const FQValues<T1> &left,
                                 const QValues<T2> &right,
                                 const Op &op,
                                 const AssemblyContext<LIBMESH_TAG> &ctx) -> FQValues<T1> {
            auto ret = left;
            std::size_t n_funs = left.size();
            std::size_t n_quad_points = left[0].size();

            for (std::size_t i = 0; i != n_funs; ++i) {
                for (std::size_t qp = 0; qp != n_quad_points; ++qp) {
                    ret[i][qp] = apply_binary(left[i][qp], right[qp], op, ctx);
                }
            }

            return ret;
        }

        template <typename T1, typename T2, class Op>
        static auto apply_binary(const FQValues<T1> &left,
                                 const FQValues<T2> &right,
                                 const Op &op,
                                 const AssemblyContext<LIBMESH_TAG> &ctx) -> FQValues<T1> {
            auto ret = left;
            std::size_t n_funs = left.size();
            std::size_t n_quad_points = left[0].size();

            for (std::size_t i = 0; i != n_funs; ++i) {
                for (std::size_t qp = 0; qp != n_quad_points; ++qp) {
                    ret[i][qp] = apply_binary(left[i][qp], right[i][qp], op, ctx);
                }
            }

            return ret;
        }

        template <typename T, class Op>
        static auto apply_binary(const double &left,
                                 FQValues<T> &&right,
                                 const Op &op,
                                 const AssemblyContext<LIBMESH_TAG> &ctx) -> FQValues<T> {
            std::size_t n_funs = right.size();
            std::size_t n_quad_points = right[0].size();

            for (std::size_t i = 0; i != n_funs; ++i) {
                for (std::size_t qp = 0; qp != n_quad_points; ++qp) {
                    right[i][qp] = apply_binary(left, right[i][qp], op, ctx);
                }
            }

            return std::move(right);
        }

        template <typename T, class Op>
        static auto apply_binary(const double &left,
                                 QValues<T> &&right,
                                 const Op &op,
                                 const AssemblyContext<LIBMESH_TAG> &ctx) -> QValues<T> {
            std::size_t n_quad_points = right.size();

            for (std::size_t qp = 0; qp != n_quad_points; ++qp) {
                right[qp] = apply_binary(left, right[qp], op, ctx);
            }

            return std::move(right);
        }

        static auto apply_binary(const FQValues<LMDenseMatrix> &left,
                                 const FQValues<LMDenseMatrix> &right,
                                 const Plus &op,
                                 const AssemblyContext<LIBMESH_TAG> &ctx) -> FQValues<LMDenseMatrix> {
            auto ret = left;
            std::size_t n_funs = left.size();
            std::size_t n_quad_points = left[0].size();

            for (std::size_t i = 0; i != n_funs; ++i) {
                for (std::size_t qp = 0; qp != n_quad_points; ++qp) {
                    ret[i][qp] = apply_binary(left[i][qp], right[i][qp], op, ctx);
                }
            }

            return ret;
        }

        template <typename T1, typename T2, class Op>
        static auto apply_binary(const QValues<T1> &left,
                                 const FQValues<T2> &right,
                                 const Op &op,
                                 const AssemblyContext<LIBMESH_TAG> &ctx) -> FQValues<T1> {
            FQValues<T1> ret;
            std::size_t n_funs = right.size();
            std::size_t n_quad_points = right[0].size();

            ret.resize(n_funs);

            for (std::size_t i = 0; i != n_funs; ++i) {
                ret[i].resize(n_quad_points);

                for (std::size_t qp = 0; qp != n_quad_points; ++qp) {
                    ret[i][qp] = apply_binary(left[qp], right[i][qp], op, ctx);
                }
            }

            return ret;
        }

        template <class Op>
        static LMDenseMatrix apply_binary(const LMDenseMatrix &left,
                                          const LMDenseMatrix &right,
                                          const Op &op,
                                          const AssemblyContext<LIBMESH_TAG> &ctx) {
            LMDenseMatrix ret = left;

            Write<LMDenseMatrix> w_(ret);
            Read<LMDenseMatrix> r_(right);
            each_read(left, [&right, &op, &ret](const SizeType i, const SizeType j, const double value) {
                ret.set(i, j, op.apply(value, right.get(i, j)));
            });

            return ret;
        }

        template <class Op>
        static auto apply_binary(const double &left,
                                 const double &right,
                                 const Op &op,
                                 const AssemblyContext<LIBMESH_TAG> &) -> double {
            return op.apply(left, right);
        }

        template <typename T>
        static auto apply_binary(FQValues<T> &&left, const T &val, const Minus &, const AssemblyContext<LIBMESH_TAG> &)
            -> FQValues<T> {
            for (auto &ll : left) {
                for (auto &l : ll) {
                    left -= val;
                }
            }

            return std::move(left);
        }

        template <typename T>
        static auto apply_binary(QValues<T> &&left, const T &val, const Minus &, const AssemblyContext<LIBMESH_TAG> &)
            -> QValues<T> {
            std::size_t n = left.size();

            for (std::size_t i = 0; i != n; ++i) {
                left[i] -= val;
            }

            return std::move(left);
        }

        template <typename T>
        static auto apply_binary(FQValues<T> &&left,
                                 const QValues<T> &val,
                                 const Minus &,
                                 const AssemblyContext<LIBMESH_TAG> &) -> FQValues<T> {
            const std::size_t n = left.size();
            const std::size_t n_qp = left[0].size();

            for (std::size_t i = 0; i != n; ++i) {
                for (std::size_t j = 0; j != n_qp; ++j) {
                    left[i][j] -= val[i];
                }
            }

            return std::move(left);
        }

        template <typename T>
        static auto apply_binary(const double val,
                                 const std::vector<T> &vals,
                                 const Multiplies &,
                                 const AssemblyContext<LIBMESH_TAG> &) -> std::vector<T> {
            auto ret = vals;
            for (auto &v : ret) {
                v *= val;
            }

            return ret;
        }

        template <typename T>
        static auto multiply(const std::vector<double> &scale,
                             std::vector<T> &&vals,
                             const AssemblyContext<LIBMESH_TAG> &) -> std::vector<T> {
            std::size_t n = scale.size();

            for (std::size_t i = 0; i != n; ++i) {
                vals[i] *= scale[i];
            }

            return std::move(vals);
        }

        template <typename T>
        static auto multiply(const std::vector<double> &scale,
                             std::vector<std::vector<T>> &&vals,
                             const AssemblyContext<LIBMESH_TAG> &) -> std::vector<std::vector<T>> {
            std::size_t n_quad_points = scale.size();
            std::size_t n_funs = vals.size();

            for (std::size_t i = 0; i != n_funs; ++i) {
                for (std::size_t qp = 0; qp != n_quad_points; ++qp) {
                    vals[i][qp] *= scale[qp];
                }
            }

            return std::move(vals);
        }

        template <typename T>
        static auto multiply(const QValues<T> &left,
                             const FQValues<double> &right,
                             const AssemblyContext<LIBMESH_TAG> &) -> FQValues<T> {
            FQValues<T> ret;
            std::size_t n_quad_points = left.size();
            std::size_t n_funs = right.size();

            assert(n_quad_points == right[0].size());

            ret.resize(n_funs);

            for (std::size_t i = 0; i != n_funs; ++i) {
                ret[i].resize(n_quad_points);

                for (std::size_t qp = 0; qp != n_quad_points; ++qp) {
                    ret[i][qp] = left[qp];
                    ret[i][qp] *= right[i][qp];
                }
            }

            return ret;
        }

        template <typename T>
        static auto multiply(const FQValues<double> &left,
                             const QValues<T> &right,
                             const AssemblyContext<LIBMESH_TAG> &) -> FQValues<T> {
            FQValues<T> ret;
            std::size_t n_quad_points = right.size();
            std::size_t n_funs = left.size();

            assert(n_quad_points == left[0].size());

            ret.resize(n_funs);

            for (std::size_t i = 0; i != n_funs; ++i) {
                ret[i].resize(n_quad_points);

                for (std::size_t qp = 0; qp != n_quad_points; ++qp) {
                    ret[i][qp] = right[qp];
                    ret[i][qp] *= left[i][qp];
                }
            }

            return ret;
        }

        template <typename T>
        static auto multiply(const double val, const std::vector<T> &vals, const AssemblyContext<LIBMESH_TAG> &)
            -> std::vector<T> {
            auto ret = vals;
            for (auto &v : ret) {
                v *= val;
            }

            return ret;
        }

        template <typename T>
        static auto multiply(const Matrix &left, const QValues<T> &vals, const AssemblyContext<LIBMESH_TAG> &)
            -> QValues<T> {
            auto ret = vals;
            auto n = vals.size();
            for (std::size_t i = 0; i < n; ++i) {
                multiply(left, vals[i], ret[i]);
            }

            return ret;
        }

        inline static auto multiply(const QValues<LMDenseMatrix> &left,
                                    const QValues<LMDenseVector> &right,
                                    const AssemblyContext<LIBMESH_TAG> &) -> QValues<LMDenseVector> {
            auto ret = right;
            auto n = right.size();
            for (std::size_t i = 0; i < n; ++i) {
                ret[i] = left[i] * right[i];
            }

            return ret;
        }

        template <typename T>
        static auto multiply(const double val, FQValues<T> &&vals, const AssemblyContext<LIBMESH_TAG> &)
            -> FQValues<T> {
            for (auto &v : vals) {
                for (auto &v_i : v) {
                    v_i *= val;
                }
            }

            return std::move(vals);
        }

        template <typename T>
        static auto apply_binary(const double val,
                                 const FQValues<T> &vals,
                                 const Multiplies &,
                                 const AssemblyContext<LIBMESH_TAG> &) -> FQValues<T> {
            auto ret = vals;
            for (auto &v_i : ret) {
                for (auto &v_j : v_i) {
                    v_j *= val;
                }
            }

            return ret;
        }

        template <class Op>
        static auto apply_unary(std::vector<double> &&vals, const Op &op, const AssemblyContext<LIBMESH_TAG> &)
            -> std::vector<double> {
            for (auto &v : vals) {
                v = op.apply(v);
            }

            return std::move(vals);
        }

        static auto norm2(const std::vector<LMDenseVector> &vals, const AssemblyContext<LIBMESH_TAG> &)
            -> std::vector<double> {
            auto n = vals.size();
            std::vector<double> ret(n);

            for (std::size_t i = 0; i < n; ++i) {
                ret[i] = utopia::norm2(vals[i]);
            }

            return ret;
        }

        static auto norm2(const std::vector<LMDenseMatrix> &vals, const AssemblyContext<LIBMESH_TAG> &)
            -> std::vector<double> {
            auto n = vals.size();
            std::vector<double> ret(n);

            for (std::size_t i = 0; i < n; ++i) {
                ret[i] = utopia::norm2(vals[i]);
            }

            return ret;
        }

        template <class Op>
        static auto apply_unary(FQValues<double> &&vals, const Op &op, const AssemblyContext<LIBMESH_TAG> &)
            -> FQValues<double> {
            for (auto &v : vals) {
                for (auto &vv : v) {
                    vv = op.apply(vv);
                }
            }

            return std::move(vals);
        }

        template <class Op>
        static auto apply_unary(FQValues<LMDenseVector> &&vals, const Op &op, const AssemblyContext<LIBMESH_TAG> &)
            -> FQValues<LMDenseVector> {
            for (auto &v : vals) {
                for (auto &vv : v) {
                    // for(int i = 0; i < LIBMESH_DIM; ++i) {
                    //     vv(i) = op.apply(vv(i));
                    // }

                    vv.transform(op);
                }
            }

            return std::move(vals);
        }

        template <class Op>
        static auto apply_unary(LMDenseVector &&vals, const Op &op, const AssemblyContext<LIBMESH_TAG> &)
            -> LMDenseVector {
            vals.transform(op);
            return std::move(vals);
        }

        // vector fe functions
        static void fun_aux(ProductFunctionSpace<LibMeshFunctionSpace> &space,
                            std::vector<std::unique_ptr<FE>> &fe_object,
                            AssemblyContext<LIBMESH_TAG> &ctx,
                            VectorFunctionType &ret) {
            assert(!fe_object.empty());

            const auto &sub_0 = space[0];

            const std::size_t n_quad_points = fe_object[sub_0.subspace_id()]->get_phi()[0].size();

            unsigned int n_shape_functions = 0;
            space.each([&n_shape_functions, &fe_object](const int index, const LibMeshFunctionSpace &subspace) {
                n_shape_functions += fe_object[subspace.subspace_id()]->n_shape_functions();
            });

            ret.resize(n_shape_functions);
            for (std::size_t i = 0; i < n_shape_functions; ++i) {
                ret[i].resize(n_quad_points);
                for (auto &r : ret[i]) {
                    r.resize(space.n_subspaces());
                    r.set(0.0);
                }
            }

            for (std::size_t qp = 0; qp < n_quad_points; ++qp) {
                int offset = 0;
                space.each([&offset, &fe_object, &ret, qp](const int sub_index, const LibMeshFunctionSpace &s) {
                    const auto &fe = fe_object[s.subspace_id()];
                    assert((static_cast<bool>(fe)));
                    const auto &fun = fe->get_phi();
                    uint n_shape_i = fe->n_shape_functions();

                    for (uint j = 0; j < n_shape_i; ++j) {
                        ret[offset++][qp].set(sub_index, fun[j][qp]);
                    }
                });
            }
        }

        static VectorFunctionType fun(const TestFunction<ProductFunctionSpace<LibMeshFunctionSpace>> &fun,
                                      AssemblyContext<LIBMESH_TAG> &ctx) {
            auto space_ptr = fun.space_ptr();
            VectorFunctionType ret;

            fun_aux(*space_ptr, ctx.test(), ctx, ret);
            return ret;
        }

        static VectorFunctionType fun(const TrialFunction<ProductFunctionSpace<LibMeshFunctionSpace>> &fun,
                                      AssemblyContext<LIBMESH_TAG> &ctx) {
            auto space_ptr = fun.space_ptr();
            VectorFunctionType ret;
            fun_aux(*space_ptr, ctx.trial(), ctx, ret);
            return ret;
        }

        //////////////////////////////////////////////////////////////////////////////////////////

        // Gradient
        // scalar fe functions
        static const GradientType &grad(const TrialFunction<LibMeshFunctionSpace> &fun,
                                        AssemblyContext<LIBMESH_TAG> &ctx) {
            // return ctx.trial()[fun.space_ptr()->subspace_id()]->get_dphi();
            assert(!ctx.grad(fun.space_ptr()->subspace_id()).empty());
            return ctx.grad(fun.space_ptr()->subspace_id());
        }

        static const GradientType &grad(const TestFunction<LibMeshFunctionSpace> &fun,
                                        AssemblyContext<LIBMESH_TAG> &ctx) {
            // return ctx.test()[fun.space_ptr()->subspace_id()]->get_dphi();
            assert(!ctx.grad(fun.space_ptr()->subspace_id()).empty());
            return ctx.grad(fun.space_ptr()->subspace_id());
        }

        // vector fe functions
        static void grad_aux(ProductFunctionSpace<LibMeshFunctionSpace> &space,
                             std::vector<std::unique_ptr<FE>> &fe_object,
                             AssemblyContext<LIBMESH_TAG> &ctx,
                             JacobianType &ret) {
            assert(!fe_object.empty());
            const auto &sub_0 = space[0];
            const std::size_t n_quad_points = fe_object[sub_0.subspace_id()]->get_phi()[0].size();
            const uint dim = space[0].mesh().spatial_dimension();

            uint n_shape_functions = 0;
            space.each([&fe_object, &n_shape_functions](const int, const LibMeshFunctionSpace &subspace) {
                n_shape_functions += fe_object[subspace.subspace_id()]->n_shape_functions();
            });

            ret.resize(n_shape_functions);
            for (std::size_t i = 0; i < n_shape_functions; ++i) {
                ret[i].resize(n_quad_points);
                // TensorValue is by default initialized to 0s

                for (std::size_t qp = 0; qp < n_quad_points; ++qp) {
                    ret[i][qp].resize(space.n_subspaces(), dim);
                    ret[i][qp].set(0.0);
                }
            }

            for (std::size_t qp = 0; qp < n_quad_points; ++qp) {
                int offset = 0;
                space.each(
                    [&offset, &fe_object, &ret, &space, qp, dim](const int sub_index, const LibMeshFunctionSpace &s) {
                        const auto &fe = fe_object[s.subspace_id()];
                        const uint n_shape_i = fe->n_shape_functions();

                        for (uint j = 0; j < n_shape_i; ++j, offset++) {
                            const auto &grad = fe->get_dphi()[j][qp];

                            for (uint d = 0; d < dim; ++d) {
                                ret[offset][qp].set(sub_index, d, grad(d));
                            }
                        }
                    });
            }
        }

        // vfe
        static const JacobianType &grad(const TrialFunction<ProductFunctionSpace<LibMeshFunctionSpace>> &fun,
                                        AssemblyContext<LIBMESH_TAG> &ctx) {
            return ctx.vector_fe()[fun.space_ptr()->subspace(0).subspace_id()]->grad;
        }

        // vfe
        static const JacobianType &grad(const TestFunction<ProductFunctionSpace<LibMeshFunctionSpace>> &fun,
                                        AssemblyContext<LIBMESH_TAG> &ctx) {
            return ctx.vector_fe()[fun.space_ptr()->subspace(0).subspace_id()]->grad;
        }

        //////////////////////////////////////////////////////////////////////////////////////////

        // Divergence
        // vector fe functions
        static void div_aux(ProductFunctionSpace<LibMeshFunctionSpace> &space,
                            std::vector<std::unique_ptr<FE>> &fe_object,
                            AssemblyContext<LIBMESH_TAG> &ctx,
                            FunctionType &ret) {
            const auto &sub_0 = space[0];
            ret.resize(fe_object[sub_0.subspace_id()]->get_dphi().size());

            assert(!fe_object.empty());
            const std::size_t n_quad_points = fe_object[sub_0.subspace_id()]->get_dphi()[0].size();
            const uint dim = space[0].mesh().spatial_dimension();

            uint n_shape_functions = 0;
            space.each([&fe_object, &n_shape_functions](const int, const LibMeshFunctionSpace &subspace) {
                n_shape_functions += fe_object[subspace.subspace_id()]->n_shape_functions();
            });

            ret.resize(n_shape_functions);
            for (std::size_t i = 0; i < n_shape_functions; ++i) {
                ret[i].resize(n_quad_points);
                // TensorValue is by default initialized to 0s
            }

            for (std::size_t qp = 0; qp < n_quad_points; ++qp) {
                int offset = 0;
                space.each(
                    [&offset, &fe_object, &ret, &space, qp, dim](const int sub_index, const LibMeshFunctionSpace &s) {
                        const auto &fe = fe_object[s.subspace_id()];
                        const uint n_shape_i = fe->n_shape_functions();

                        for (uint j = 0; j < n_shape_i; ++j, offset++) {
                            const auto &grad = fe->get_dphi()[j][qp];
                            ret[offset][qp] = grad(sub_index);
                        }
                    });
            }
        }

        static FunctionType div(const TrialFunction<ProductFunctionSpace<LibMeshFunctionSpace>> &fun,
                                AssemblyContext<LIBMESH_TAG> &ctx) {
            auto space_ptr = fun.space_ptr();
            FunctionType ret;
            div_aux(*space_ptr, ctx.trial(), ctx, ret);
            return ret;
        }

        static FunctionType div(const TestFunction<ProductFunctionSpace<LibMeshFunctionSpace>> &fun,
                                AssemblyContext<LIBMESH_TAG> &ctx) {
            auto space_ptr = fun.space_ptr();
            FunctionType ret;
            div_aux(*space_ptr, ctx.test(), ctx, ret);
            return ret;
        }

        //////////////////////////////////////////////////////////////////////////////////////////
        // Curl
        static void eval_curl(const uint dim, const LMDenseMatrix &deriv, LMDenseVector &result) {
            result.resize(3);

            switch (dim) {
                case 2: {
                    result.set(0.0);
                    result.set(2, deriv.get(1, 0) - deriv.get(0, 1));
                    break;
                }

                case 3: {
                    result.set(0, deriv.get(2, 1) - deriv.get(1, 2));
                    result.set(1, deriv.get(0, 2) - deriv.get(2, 0));
                    result.set(2, deriv.get(1, 0) - deriv.get(0, 1));
                    break;
                }

                default: {
                    assert(false);
                }
            }
        }

        static void curl_aux(ProductFunctionSpace<LibMeshFunctionSpace> &space,
                             std::vector<std::unique_ptr<FE>> &fe_object,
                             AssemblyContext<LIBMESH_TAG> &ctx,
                             CurlType &ret) {
            const auto &sub_0 = space[0];
            JacobianType grads;

            grad_aux(space, fe_object, ctx, grads);

            const uint n_shape_x = fe_object[sub_0.subspace_id()]->n_shape_functions();
            uint n_shape_functions = 0;

            space.each([&n_shape_functions, &fe_object, &n_shape_x](const int, LibMeshFunctionSpace &subspace) {
                assert(n_shape_x == fe_object[subspace.subspace_id()]->n_shape_functions());
                n_shape_functions += fe_object[subspace.subspace_id()]->n_shape_functions();
            });

            const std::size_t n_subspaces = space.n_subspaces();
            const std::size_t n_quad_points = grads[0].size();
            const std::size_t dim = sub_0.mesh().spatial_dimension();

            ret.resize(n_shape_functions);
            for (auto &r : ret) {
                r.resize(n_quad_points);
            }

            for (std::size_t qp = 0; qp < n_quad_points; ++qp) {
                uint offset = 0;
                for (uint i = 0; i < n_subspaces; ++i) {
                    for (uint j = 0; j < n_shape_x; ++j, ++offset) {
                        auto &ret_j = ret[offset][qp];
                        auto &grad_j = grads[offset][qp];
                        eval_curl(dim, grad_j, ret_j);
                    }
                }
            }
        }

        static CurlType curl(const TrialFunction<ProductFunctionSpace<LibMeshFunctionSpace>> &fun,
                             AssemblyContext<LIBMESH_TAG> &ctx) {
            auto space_ptr = fun.space_ptr();
            CurlType ret;
            curl_aux(*space_ptr, ctx.trial(), ctx, ret);
            return ret;
        }

        static CurlType curl(const TestFunction<ProductFunctionSpace<LibMeshFunctionSpace>> &fun,
                             AssemblyContext<LIBMESH_TAG> &ctx) {
            auto space_ptr = fun.space_ptr();
            CurlType ret;
            curl_aux(*space_ptr, ctx.test(), ctx, ret);
            return ret;
        }

        //////////////////////////////////////////////////////////////////////////////////////////
        // Interpolate
        static void gather_interp_values(const UVector &c,
                                         const TrialFunction<LibMeshFunctionSpace> &f,
                                         Vector &element_values,
                                         AssemblyContext<LIBMESH_TAG> &ctx) {
            // auto &c   = interp.coefficient();
            // auto &f   = interp.fun();

            auto space_ptr = f.space_ptr();
            const auto &mesh = space_ptr->mesh();
            const auto &elem_ptr = utopia::elem_ptr(mesh, ctx.current_element());
            const auto &dof_map = space_ptr->dof_map();

            std::vector<libMesh::dof_id_type> indices;
            dof_map.dof_indices(elem_ptr, indices, space_ptr->subspace_id());

            element_values = zeros(indices.size());

            Write<Vector> w(element_values);
            Read<UVector> r(c);

            typename Traits<UVector>::IndexSet u_index;
            u_index.insert(u_index.end(), indices.begin(), indices.end());

            // for(std::size_t i = 0; i < indices.size(); ++i) {
            // 	element_values.set(i, c.get(indices[i]));
            // }
            // std::cout << raw_type(c) << std::endl;
            assert(c.has_ghosts() || c.comm().size() == 1);
            c.get(u_index, element_values.entries());
        }

        template <class Space>
        static std::vector<Scalar> fun(const UVector &c,
                                       const TrialFunction<Space> &f,
                                       AssemblyContext<LIBMESH_TAG> &ctx) {
            // auto &c  = interp.coefficient();
            // auto &f  = interp.fun();
            auto &&g = fun(f, ctx);

            Vector element_values;
            gather_interp_values(c, f, element_values, ctx);

            const std::size_t n_shape_functions = g.size();
            const std::size_t n_quad_points = g[0].size();

            std::vector<Scalar> ret(n_quad_points, 0.);

            for (std::size_t i = 0; i < n_shape_functions; ++i) {
                for (std::size_t qp = 0; qp < n_quad_points; ++qp) {
                    ret[qp] += g[i][qp] * element_values.get(i);
                }
            }

            return ret;
        }

        template <class C>
        static auto fun(const Interpolate<C, TrialFunction<ProductFunctionSpace<LibMeshFunctionSpace>>> &f,
                        AssemblyContext<LIBMESH_TAG> &ctx) -> QValues<LMDenseVector> {
            return fun(f.coefficient(), f.fun(), ctx);
        }

        template <class Space>
        static std::vector<LMDenseVector> fun(const UVector &c,
                                              const TrialFunction<ProductFunctionSpace<Space>> &f,
                                              AssemblyContext<LIBMESH_TAG> &ctx) {
            // auto &c  = interp.coefficient();
            // auto &f  = interp.fun();
            auto &&g = fun(f, ctx);

            Vector element_values;
            gather_interp_values(c.derived(), f, element_values, ctx);

            const std::size_t n_shape_functions = g.size();
            const std::size_t n_quad_points = g[0].size();

            std::vector<LMDenseVector> ret(n_quad_points);

            for (std::size_t qp = 0; qp < n_quad_points; ++qp) {
                ret[qp].resize(f.codim());
                ret[qp].set(0.0);

                for (std::size_t i = 0; i < n_shape_functions; ++i) {
                    if (std::is_rvalue_reference<decltype(g)>::value) {
                        g[i][qp] *= element_values.get(i);
                        ret[qp] += g[i][qp];
                    } else {
                        auto temp = g[i][qp];
                        temp *= element_values.get(i);
                        ret[qp] += temp;
                    }
                }
            }

            return ret;
        }

        static std::vector<Vector> grad(const UVector &c,
                                        const TrialFunction<LibMeshFunctionSpace> &f,
                                        AssemblyContext<LIBMESH_TAG> &ctx) {
            // auto &c  = interp.coefficient();
            // auto &f  = interp.fun();
            auto &&g = grad(f, ctx);

            Vector element_values;
            gather_interp_values(c, f, element_values, ctx);

            const std::size_t n_shape_functions = g.size();
            const std::size_t n_quad_points = g[0].size();

            std::vector<Vector> ret(n_quad_points);

            // const uint dim = f.space_ptr()->mesh().mesh_dimension();
            const uint dim = f.space_ptr()->mesh().spatial_dimension();

            for (auto &r : ret) {
                r = zeros(dim);
            }

            for (std::size_t i = 0; i < n_shape_functions; ++i) {
                for (std::size_t qp = 0; qp < n_quad_points; ++qp) {
                    ret[qp] += element_values.get(i) * g[i][qp];
                    // add(ret[qp], element_values.get(i) * g[i][qp]);
                }
            }

            return ret;
        }

        static auto grad(const Interpolate<const UVector, TrialFunction<ProductFunctionSpace<LibMeshFunctionSpace>>> &f,
                         AssemblyContext<LIBMESH_TAG> &ctx) -> std::vector<Matrix>
        // std::vector<TensorValueT>
        {
            return grad(f.coefficient(), f.fun(), ctx);
        }

        static auto div(const Interpolate<const UVector, TrialFunction<ProductFunctionSpace<LibMeshFunctionSpace>>> &f,
                        AssemblyContext<LIBMESH_TAG> &ctx) -> std::vector<double>
        // std::vector<TensorValueT>
        {
            return div(f.coefficient(), f.fun(), ctx);
        }

        static auto div(const Interpolate<UVector, TrialFunction<ProductFunctionSpace<LibMeshFunctionSpace>>> &f,
                        AssemblyContext<LIBMESH_TAG> &ctx) -> std::vector<double>
        // std::vector<TensorValueT>
        {
            return div(f.coefficient(), f.fun(), ctx);
        }

        static std::vector<Matrix> grad(const UVector &c,
                                        const TrialFunction<ProductFunctionSpace<LibMeshFunctionSpace>> &f,
                                        AssemblyContext<LIBMESH_TAG> &ctx) {
            // auto &f   = interp.fun();
            auto space_ptr = f.space_ptr();
            auto &&g = grad(f, ctx);

            Vector element_values;
            gather_interp_values(c, f, element_values, ctx);

            const SizeType rows = space_ptr->n_subspaces();
            const SizeType cols = space_ptr->subspace(0).mesh().spatial_dimension();
            Size s{rows, cols};

            const std::size_t n_shape_functions = g.size();
            const std::size_t n_quad_points = g[0].size();

            std::vector<Matrix> ret(n_quad_points);

            auto ml = serial_layout(rows, cols);
            for (auto &r : ret) {
                r.dense(ml, 0.0);
            }

            for (std::size_t i = 0; i < n_shape_functions; ++i) {
                for (std::size_t qp = 0; qp < n_quad_points; ++qp) {
                    // add(ret[qp], element_values.get(i) * g[i][qp]);
                    ret[qp] += element_values.get(i) * g[i][qp];
                }
            }

            return ret;
        }

        static std::vector<double> div(const UVector &c,
                                       const TrialFunction<ProductFunctionSpace<LibMeshFunctionSpace>> &f,
                                       AssemblyContext<LIBMESH_TAG> &ctx) {
            auto space_ptr = f.space_ptr();
            auto &&g = div(f, ctx);

            Vector element_values;
            gather_interp_values(c, f, element_values, ctx);

            const SizeType rows = space_ptr->n_subspaces();
            const SizeType cols = space_ptr->subspace(0).mesh().spatial_dimension();
            Size s{rows, cols};

            const std::size_t n_shape_functions = g.size();
            const std::size_t n_quad_points = g[0].size();

            std::vector<double> ret(n_quad_points, 0.0);

            for (std::size_t i = 0; i < n_shape_functions; ++i) {
                for (std::size_t qp = 0; qp < n_quad_points; ++qp) {
                    ret[qp] += g[i][qp] * element_values.get(i);
                }
            }

            return ret;
        }

        static void gather_interp_values(const UVector &c,
                                         const TrialFunction<ProductFunctionSpace<LibMeshFunctionSpace>> &f,
                                         Vector &element_values,
                                         AssemblyContext<LIBMESH_TAG> &ctx) {
            using IndexArray = typename Traits<UVector>::IndexArray;

            // auto &f   = interp.fun();

            auto space_ptr = f.space_ptr();
            const auto &sub_0 = space_ptr->subspace(0);
            const auto &mesh = sub_0.mesh();
            const auto &dof_map = sub_0.dof_map();

            const auto &elem_ptr = utopia::elem_ptr(mesh, ctx.current_element());

            IndexArray prod_indices;
            std::vector<libMesh::dof_id_type> indices;

            space_ptr->each([&](const int sub_index, const LibMeshFunctionSpace &space) {
                dof_map.dof_indices(elem_ptr, indices, space.subspace_id());
                prod_indices.insert(prod_indices.end(), indices.begin(), indices.end());
            });

            const std::size_t n_indices = prod_indices.size();
            element_values = zeros(n_indices);

            Write<Vector> w(element_values);
            Read<UVector> r(c);
            assert(c.has_ghosts() || c.comm().size() == 1);
            c.get(prod_indices, element_values.entries());
        }

        //////////////////////////////////////////////////////////////////////////////////////////

        inline static auto apply_binary(const Number<double> &left,
                                        const Gradient<TrialFunction<LibMeshFunctionSpace>> &right,
                                        const Multiplies &op,
                                        AssemblyContext<LIBMESH_TAG> &ctx) -> GradientType {
            return multiply(static_cast<double>(left), right, ctx);
        }

        template <class Space, class Op>
        inline static auto apply_binary(const Number<double> &left,
                                        const TrialFunction<Space> &right,
                                        const Op &op,
                                        AssemblyContext<LIBMESH_TAG> &ctx) ->
            typename remove_ref_and_const<decltype(fun(right, ctx))>::type {
            typename remove_ref_and_const<decltype(fun(right, ctx))>::type ret = fun(right, ctx);

            const double num = left;
            for (auto &v : ret) {
                for (auto &s : v) {
                    s = op.apply(num, s);
                }
            }

            return ret;
        }

        template <class Space, class Op>
        inline static auto apply_binary(const Number<double> &left,
                                        const TestFunction<Space> &right,
                                        const Op &op,
                                        AssemblyContext<LIBMESH_TAG> &ctx) ->
            typename remove_ref_and_const<decltype(fun(right, ctx))>::type {
            typename remove_ref_and_const<decltype(fun(right, ctx))>::type ret = fun(right, ctx);

            const double num = left;
            for (auto &v : ret) {
                for (auto &s : v) {
                    s = op.apply(num, s);
                }
            }

            return ret;
        }

        template <typename T, typename C>
        inline static auto inner(const Interpolate<C, TrialFunction<LibMeshFunctionSpace>> &right,
                                 const ConstantCoefficient<T, 0> &left,
                                 AssemblyContext<LIBMESH_TAG> &ctx) -> QValues<double> {
            return apply_binary(left, right, Multiplies(), ctx);
        }

        template <typename T, typename C>
        inline static auto inner(const Interpolate<C, TrialFunction<LibMeshFunctionSpace>> &left,
                                 const Interpolate<C, TrialFunction<LibMeshFunctionSpace>> &right,
                                 AssemblyContext<LIBMESH_TAG> &ctx) -> QValues<double> {
            auto f_left = fun(left.coefficient(), left.fun(), ctx);
            auto &&f_right = fun(right.coefficient(), right.fun(), ctx);

            std::size_t n = f_left.size();

            for (std::size_t i = 0; i < n; ++i) {
                f_left[i] *= f_right[i];
            }

            return f_left;
        }

        template <typename T>
        inline static auto inner(QValues<double> &&left,
                                 const ConstantCoefficient<T, 0> &right,
                                 AssemblyContext<LIBMESH_TAG> &ctx) -> std::vector<double> {
            for (auto &l : left) {
                l *= right.expr();
            }

            return std::move(left);
        }

        template <typename T>
        inline static auto inner(const ConstantCoefficient<T, 0> &left,
                                 QValues<double> &&right,
                                 AssemblyContext<LIBMESH_TAG> &ctx) -> std::vector<double> {
            for (auto &r : right) {
                r *= left.expr();
            }

            return std::move(right);
        }

        template <typename T, typename C>
        inline static auto inner(const ConstantCoefficient<T, 0> &left,
                                 const Interpolate<C, TrialFunction<LibMeshFunctionSpace>> &right,
                                 AssemblyContext<LIBMESH_TAG> &ctx) -> std::vector<double> {
            return apply_binary(left, right, Multiplies(), ctx);
        }

        template <typename T, typename C>
        inline static auto apply_binary(const ConstantCoefficient<T, 0> &left,
                                        const Interpolate<C, TrialFunction<LibMeshFunctionSpace>> &right,
                                        const Multiplies &,
                                        AssemblyContext<LIBMESH_TAG> &ctx) -> std::vector<double> {
            auto &&f = fun(right, ctx);

            for (auto &f_v : f) {
                f_v *= left.expr();
            }

            return std::move(f);
        }

        template <typename T, typename C>
        inline static auto apply_binary(const ConstantCoefficient<T, 0> &left,
                                        const Interpolate<C, TrialFunction<LibMeshFunctionSpace>> &right,
                                        const Minus &,
                                        AssemblyContext<LIBMESH_TAG> &ctx) -> std::vector<double> {
            auto &&f = fun(right.coefficient(), right.fun(), ctx);

            for (auto &f_v : f) {
                f_v = left.expr() - f_v;
            }

            return std::move(f);
        }

        template <typename T, typename C>
        inline static auto apply_binary(const Interpolate<C, TrialFunction<LibMeshFunctionSpace>> &left,
                                        const ConstantCoefficient<T, 0> &right,
                                        const Minus &,
                                        AssemblyContext<LIBMESH_TAG> &ctx) -> std::vector<double> {
            auto &&f = fun(left.coefficient(), left.fun(), ctx);

            for (auto &f_v : f) {
                f_v -= right.expr();
            }

            return std::move(f);
        }

        template <typename T>
        inline static auto apply_binary(const ConstantCoefficient<T, 0> &left,
                                        std::vector<T> &&values,
                                        const Plus &,
                                        AssemblyContext<LIBMESH_TAG> &ctx) -> std::vector<T> {
            for (auto &v : values) {
                v += left.expr();
            }

            return std::move(values);
        }

        template <typename T, typename C>
        inline static auto apply_binary(const ConstantCoefficient<T, 0> &left,
                                        const Interpolate<C, TrialFunction<LibMeshFunctionSpace>> &right,
                                        const Plus &,
                                        AssemblyContext<LIBMESH_TAG> &ctx) -> std::vector<double> {
            auto &&f = fun(right.coefficient(), right.fun(), ctx);

            for (auto &f_v : f) {
                f_v += left.expr();
            }

            return std::move(f);
        }

        template <typename T, typename C>
        inline static auto apply_binary(const ConstantCoefficient<T, 0> &left,
                                        const Interpolate<C, TrialFunction<LibMeshFunctionSpace>> &right,
                                        const Divides &,
                                        AssemblyContext<LIBMESH_TAG> &ctx) -> std::vector<double> {
            auto &&f = fun(right.coefficient(), right.fun(), ctx);

            for (auto &f_v : f) {
                f_v = left.expr() / f_v;
            }

            return std::move(f);
        }

        template <typename T>
        inline static auto apply_binary(QValues<T> &&left,
                                        const QValues<T> &right,
                                        const Divides &,
                                        AssemblyContext<LIBMESH_TAG> &ctx) -> QValues<T> {
            auto n = right.size();
            assert(n == left.size());

            for (std::size_t i = 0; i < n; ++i) {
                left[i] /= right[i];
            }

            return std::move(left);
        }

        template <typename T, int Order>
        inline static auto apply_binary(Tensor<T, Order> &&left,
                                        const Tensor<T, Order> &right,
                                        const Divides &,
                                        AssemblyContext<LIBMESH_TAG> &ctx) -> T {
            left.derived() /= right.derived();
            return std::move(left.derived());
        }

        template <typename T>
        inline static auto apply_binary(FQValues<T> &&left,
                                        const FQValues<T> &right,
                                        const Divides &op,
                                        AssemblyContext<LIBMESH_TAG> &ctx) -> FQValues<T> {
            auto n = right.size();
            auto n_quad_points = right[0].size();

            assert(n == left.size());
            assert(n_quad_points == left[0].size());

            for (std::size_t i = 0; i < n; ++i) {
                for (std::size_t qp = 0; i < n_quad_points; ++i) {
                    // left[i][qp] = apply_binary(left[i][qp], right[i][qp], op, ctx);
                    left[i][qp] /= right[i][qp];
                }
            }

            return std::move(left);
        }

        template <typename T>
        inline static auto apply_binary(const QValues<T> &left,
                                        const FQValues<double> &right,
                                        const Divides &op,
                                        AssemblyContext<LIBMESH_TAG> &ctx) -> FQValues<T> {
            auto n = right.size();
            auto n_quad_points = right[0].size();

            FQValues<T> ret;
            ret.resize(n);

            assert(n == left.size());
            assert(n_quad_points == left[0].size());

            for (std::size_t i = 0; i < n; ++i) {
                ret[i].resize(n_quad_points);

                for (std::size_t qp = 0; i < n_quad_points; ++i) {
                    ret[i][qp] = left[qp] / right[i][qp];
                }
            }

            return ret;
        }

        inline static auto apply_binary(QValues<LMDenseMatrix> &&left,
                                        const QValues<double> &right,
                                        const Divides &op,
                                        AssemblyContext<LIBMESH_TAG> &ctx) -> QValues<LMDenseMatrix> {
            auto n = right.size();
            assert(n == left.size());
            for (std::size_t i = 0; i < n; ++i) {
                left[i] /= right[i];
            }

            return std::move(left);
        }

        template <class Left, class Right>
        inline static auto grad_t_plus_grad(const double scaling,
                                            const Left &left,
                                            const Right &right,
                                            AssemblyContext<LIBMESH_TAG> &ctx) ->
            typename remove_ref_and_const<decltype(grad(right, ctx))>::type {
            auto ret = grad(right, ctx);
            auto &&grad_left = grad(left, ctx);

            grad_t_plus_grad_aux(scaling, grad_left, ret);
            return std::move(ret);
        }

        template <class C>
        inline static auto grad_t_plus_grad(
            const double scaling,
            const GradInterpolate<C, TrialFunction<ProductFunctionSpace<LibMeshFunctionSpace>>> &left,
            const GradInterpolate<C, TrialFunction<ProductFunctionSpace<LibMeshFunctionSpace>>> &right,
            AssemblyContext<LIBMESH_TAG> &ctx) -> QValues<LMDenseMatrix> {
            const auto &l = left.expr();
            const auto &r = right.expr();

            if (l.fun().equals(r.fun())) {
                auto &&grad_left = grad(l, ctx);
                auto ret = grad_left;
                grad_t_plus_grad_aux(scaling, grad_left, ret);
                return std::move(ret);
            } else {
                auto ret = grad(r, ctx);
                auto &&grad_left = grad(l, ctx);
                grad_t_plus_grad_aux(scaling, grad_left, ret);
                return std::move(ret);
            }
        }

        static void grad_t_plus_grad_aux(const double scaling,
                                         const QValues<LMDenseMatrix> &left,
                                         QValues<LMDenseMatrix> &right) {
            for (std::size_t qp = 0; qp < left.size(); ++qp) {
                right[qp] += utopia::transpose(left[qp]);
                right[qp] *= scaling;
            }
        }

        static void grad_t_plus_grad_aux(const double scaling,
                                         const FQValues<LMDenseMatrix> &left,
                                         FQValues<LMDenseMatrix> &right) {
            for (std::size_t i = 0; i < left.size(); ++i) {
                for (std::size_t qp = 0; qp < left[i].size(); ++qp) {
                    right[i][qp] += utopia::transpose(left[i][qp]);
                    right[i][qp] *= scaling;
                }
            }
        }

        //////////////////////////////////////////////////////////////////////////////////////////////////////
        template <typename C>
        inline static auto multiply(const QValues<double> &left,
                                    const GradInterpolate<C, TrialFunction<LibMeshFunctionSpace>> &right,
                                    AssemblyContext<LIBMESH_TAG> &ctx)
            -> decltype(grad(right.expr().coefficient(), right.expr().fun(), ctx)) {
            auto &&g = grad(right.expr().coefficient(), right.expr().fun(), ctx);

            std::size_t n_quad_points = g.size();

            for (std::size_t qp = 0; qp < n_quad_points; ++qp) {
                g[qp] *= left[qp];
            }

            return std::move(g);
        }

        template <typename C>
        inline static auto multiply(const QValues<LMDenseMatrix> &left,
                                    const GradInterpolate<C, TrialFunction<LibMeshFunctionSpace>> &right,
                                    AssemblyContext<LIBMESH_TAG> &ctx)
            -> decltype(grad(right.expr().coefficient(), right.expr().fun(), ctx)) {
            auto g = grad(right.expr().coefficient(), right.expr().fun(), ctx);
            auto ret = g;

            std::size_t n_quad_points = g.size();

            for (std::size_t qp = 0; qp < n_quad_points; ++qp) {
                multiply(left[qp], g[qp], ret[qp]);
            }

            return std::move(ret);
        }

        template <typename T, typename C>
        inline static auto multiply(const ConstantCoefficient<T, 0> &left,
                                    const GradInterpolate<C, TrialFunction<LibMeshFunctionSpace>> &right,
                                    AssemblyContext<LIBMESH_TAG> &ctx)
            -> decltype(grad(right.expr().coefficient(), right.expr().fun(), ctx)) {
            auto &&g = grad(right.expr().coefficient(), right.expr().fun(), ctx);

            std::size_t n_quad_points = g.size();

            for (std::size_t qp = 0; qp < n_quad_points; ++qp) {
                g[qp] *= left.expr();
            }

            return std::move(g);
        }

        template <typename C>
        inline static auto multiply(const LMDenseMatrix &mat,
                                    const GradInterpolate<C, TrialFunction<LibMeshFunctionSpace>> &right,
                                    AssemblyContext<LIBMESH_TAG> &ctx)
            -> decltype(grad(right.expr().coefficient(), right.expr().fun(), ctx)) {
            auto &&g = grad(right.expr().coefficient(), right.expr().fun(), ctx);
            auto ret = g;

            std::size_t n_quad_points = g.size();

            for (std::size_t qp = 0; qp < n_quad_points; ++qp) {
                multiply(mat, g[qp], ret[qp]);
            }

            return std::move(ret);
        }

        template <typename T>
        inline static auto multiply(const std::vector<std::vector<double>> &left,
                                    const ConstantCoefficient<T, 0> &right,
                                    const AssemblyContext<LIBMESH_TAG> &ctx) -> std::vector<std::vector<double>> {
            auto ret = left;
            for (auto &r : ret) {
                for (auto &v : r) {
                    v *= right.expr();
                }
            }

            return ret;
        }

        inline static auto multiply(const std::vector<std::vector<double>> &left,
                                    const std::vector<Scalar> &right,
                                    const AssemblyContext<LIBMESH_TAG> &ctx) -> std::vector<std::vector<double>> {
            auto ret = left;
            for (auto &r : ret) {
                for (std::size_t i = 0; i < right.size(); ++i) {
                    r[i] *= right[i];
                }
            }

            return ret;
        }

        template <class T>
        inline static auto multiply(const QValues<double> &left,
                                    const T &right,
                                    const AssemblyContext<LIBMESH_TAG> &ctx) -> QValues<T> {
            QValues<T> ret(left.size());

            auto n_qp = left.size();

            for (std::size_t qp = 0; qp < n_qp; ++qp) {
                ret[qp] = left[qp] * right;
            }

            return ret;
        }

        inline static auto multiply(const FQValues<LMDenseMatrix> &left,
                                    const QValues<LMDenseMatrix> &right,
                                    const AssemblyContext<LIBMESH_TAG> &ctx) -> FQValues<LMDenseMatrix> {
            FQValues<LMDenseMatrix> ret(left.size());

            for (std::size_t i = 0; i < left.size(); ++i) {
                ret[i].resize(left[i].size());

                for (std::size_t qp = 0; qp < right.size(); ++qp) {
                    multiply(left[i][qp], right[qp], ret[i][qp]);
                }
            }

            return ret;
        }

        inline static auto multiply(const QValues<LMDenseMatrix> &left,
                                    const FQValues<LMDenseMatrix> &right,
                                    const AssemblyContext<LIBMESH_TAG> &ctx) -> FQValues<LMDenseMatrix> {
            FQValues<LMDenseMatrix> ret(right.size());

            for (std::size_t i = 0; i < right.size(); ++i) {
                ret[i].resize(right[i].size());

                for (std::size_t qp = 0; qp < left.size(); ++qp) {
                    multiply(left[qp], right[i][qp], ret[i][qp]);
                }
            }

            return ret;
        }

        template <class Tensor>
        inline static auto multiply(const QValues<Tensor> &left,
                                    const QValues<Tensor> &right,
                                    const AssemblyContext<LIBMESH_TAG> &ctx) -> QValues<Tensor> {
            assert(left.size() == right.size());

            const auto n_quad_points = right.size();

            QValues<Tensor> ret(n_quad_points);

            for (std::size_t qp = 0; qp < n_quad_points; ++qp) {
                multiply(left[qp], right[qp], ret[qp]);
            }

            return ret;
        }

        inline static auto multiply(const QValues<double> &left,
                                    const QValues<double> &right,
                                    const AssemblyContext<LIBMESH_TAG> &ctx) -> QValues<double> {
            auto ret = left;
            for (std::size_t i = 0; i < right.size(); ++i) {
                ret[i] *= right[i];
            }

            return ret;
        }

        template <typename T>
        inline static auto multiply(const QValues<double> &left,
                                    const ConstantCoefficient<T, 0> &right,
                                    const AssemblyContext<LIBMESH_TAG> &ctx) -> QValues<double> {
            auto ret = left;
            for (auto &v : ret) {
                v *= right.expr();
            }

            return ret;
        }

        inline static auto multiply(const double scale, LMDenseMatrix &&tensor, const AssemblyContext<LIBMESH_TAG> &)
            -> LMDenseMatrix {
            tensor *= scale;
            return std::move(tensor);
        }

        template <class Space>
        inline static auto multiply(const Scalar &left,
                                    const Gradient<TrialFunction<Space>> &right,
                                    AssemblyContext<LIBMESH_TAG> &ctx) ->
            typename remove_ref_and_const<decltype(grad(right.expr(), ctx))>::type {
            typename remove_ref_and_const<decltype(grad(right.expr(), ctx))>::type ret = grad(right.expr(), ctx);

            for (auto &v : ret) {
                for (auto &s : v) {
                    s = left * s;
                }
            }

            return ret;
        }

        inline static auto multiply(const LMDenseMatrix &left,
                                    const Gradient<TrialFunction<LibMeshFunctionSpace>> &right,
                                    AssemblyContext<LIBMESH_TAG> &ctx) ->
            typename remove_ref_and_const<decltype(grad(right.expr(), ctx))>::type {
            const GradientType &g = grad(right.expr(), ctx);
            GradientType ret = g;

            const std::size_t n = g.size();
            const std::size_t n_qp = g[0].size();

            for (std::size_t i = 0; i < n; ++i) {
                for (std::size_t k = 0; k < n_qp; ++k) {
                    ret[i][k] = left * g[i][k];
                }
            }
            return ret;
        }

        inline static void multiply(const LMDenseMatrix &left, const double &right, LMDenseMatrix &out) {
            out = right * left;
        }

        inline static void multiply(LMDenseMatrix &&left, const double &right, LMDenseMatrix &out) {
            out = std::move(left);
            out *= right;
        }

        inline static LMDenseMatrix multiply(LMDenseMatrix &&left,
                                             const double &right,
                                             AssemblyContext<LIBMESH_TAG> &) {
            left *= right;
            return left;
        }

        inline static void multiply(const LMDenseMatrix &left, const LMDenseVector &right, LMDenseVector &out) {
            out = left * right;
        }

        inline static void multiply(const LMDenseMatrix &left, const LMDenseMatrix &right, LMDenseMatrix &out) {
            out = left * right;
        }

        // template<class Space>
        inline static auto multiply(const LMDenseMatrix &left,
                                    const Gradient<TrialFunction<ProductFunctionSpace<LibMeshFunctionSpace>>> &right,
                                    AssemblyContext<LIBMESH_TAG> &ctx) ->
            typename remove_ref_and_const<decltype(grad(right.expr(), ctx))>::type {
            typename remove_ref_and_const<decltype(grad(right.expr(), ctx))>::type ret = grad(right.expr(), ctx);

            for (auto &v : ret) {
                for (auto &s : v) {
                    auto s_copy = s;
                    multiply(left, s_copy, s);
                }
            }

            return ret;
        }

        template <class Space>
        inline static auto multiply(const QValues<double> &left,
                                    const Gradient<TrialFunction<Space>> &right,
                                    AssemblyContext<LIBMESH_TAG> &ctx) ->
            typename remove_ref_and_const<decltype(grad(right.expr(), ctx))>::type {
            typename remove_ref_and_const<decltype(grad(right.expr(), ctx))>::type ret = grad(right.expr(), ctx);

            const std::size_t n_functions = ret.size();
            const std::size_t n_quad_points = ret[0].size();

            assert(n_quad_points == left.size());

            for (std::size_t qp = 0; qp < n_quad_points; ++qp) {
                for (std::size_t i = 0; i < n_functions; ++i) {
                    ret[i][qp] *= left[qp];
                }
            }

            return ret;
        }

        template <class Left, class Space>
        inline static auto multiply(const Left &left,
                                    const Gradient<TestFunction<Space>> &right,
                                    AssemblyContext<LIBMESH_TAG> &ctx) ->
            typename remove_ref_and_const<decltype(grad(right.expr(), ctx))>::type {
            typename remove_ref_and_const<decltype(grad(right.expr(), ctx))>::type ret = grad(right.expr(), ctx);

            const std::size_t n_functions = ret.size();
            const std::size_t n_quad_points = ret[0].size();

            for (std::size_t qp = 0; qp < n_quad_points; ++qp) {
                for (std::size_t i = 0; i < n_functions; ++i) {
                    ret[i][qp] = left * ret[i][qp];
                }
            }

            return ret;
        }

        template <class Left, class Space>
        inline static auto multiply(const std::vector<std::vector<Left>> &left,
                                    const Gradient<TestFunction<Space>> &right,
                                    AssemblyContext<LIBMESH_TAG> &ctx) ->
            typename remove_ref_and_const<decltype(grad(right.expr(), ctx))>::type {
            typename remove_ref_and_const<decltype(grad(right.expr(), ctx))>::type ret = grad(right.expr(), ctx);

            const std::size_t n_functions = ret.size();
            const std::size_t n_quad_points = ret[0].size();

            for (std::size_t qp = 0; qp < n_quad_points; ++qp) {
                for (std::size_t i = 0; i < n_functions; ++i) {
                    ret[i][qp] = left[i][qp] * ret[i][qp];
                }
            }

            return ret;
        }

        template <class T1, typename T2, typename T3>
        void multiply(const LMDenseMatrix &left,
                      const libMesh::TensorValue<T2> &right,
                      libMesh::TensorValue<T3> &result) {
            auto s = size(left);
            result.zero();

            Read<LMDenseMatrix> r_l(left);

            for (SizeType i = 0; i < s.get(0); ++i) {
                for (SizeType j = 0; j < s.get(1); ++j) {
                    for (SizeType k = 0; k < LIBMESH_DIM; ++k) {
                        result(i, k) += left.get(i, j) * right(j, k);
                    }
                }
            }
        }

        template <class Left, class Space>
        inline static auto multiply(const std::vector<Left> &left,
                                    const Gradient<TestFunction<Space>> &right,
                                    AssemblyContext<LIBMESH_TAG> &ctx) ->
            typename remove_ref_and_const<decltype(grad(right.expr(), ctx))>::type {
            typename remove_ref_and_const<decltype(grad(right.expr(), ctx))>::type ret = grad(right.expr(), ctx);

            const std::size_t n_functions = ret.size();
            const std::size_t n_quad_points = ret[0].size();

            auto temp = ret[0][0];

            for (std::size_t qp = 0; qp < n_quad_points; ++qp) {
                for (std::size_t i = 0; i < n_functions; ++i) {
                    temp = ret[i][qp];
                    multiply(left[qp], temp, ret[i][qp]);
                }
            }

            return ret;
        }

        template <class Space>
        inline static auto multiply(const std::vector<double> &left,
                                    const Gradient<TestFunction<Space>> &right,
                                    AssemblyContext<LIBMESH_TAG> &ctx) ->
            typename remove_ref_and_const<decltype(grad(right.expr(), ctx))>::type {
            typename remove_ref_and_const<decltype(grad(right.expr(), ctx))>::type ret = grad(right.expr(), ctx);

            const std::size_t n_functions = ret.size();
            const std::size_t n_quad_points = ret[0].size();

            for (std::size_t qp = 0; qp < n_quad_points; ++qp) {
                for (std::size_t i = 0; i < n_functions; ++i) {
                    ret[i][qp] *= left[qp];
                }
            }

            return ret;
        }

        template <class Space>
        inline static auto multiply(const std::vector<double> &left,
                                    const TrialFunction<Space> &right,
                                    AssemblyContext<LIBMESH_TAG> &ctx) ->
            typename remove_ref_and_const<decltype(fun(right, ctx))>::type {
            typename remove_ref_and_const<decltype(fun(right, ctx))>::type ret = fun(right, ctx);

            const std::size_t n_quad_points = left.size();
            const std::size_t n_functions = ret.size();

            assert(n_quad_points == left.size());

            for (std::size_t qp = 0; qp < n_quad_points; ++qp) {
                for (std::size_t i = 0; i < n_functions; ++i) {
                    ret[i][qp] *= left[qp];
                }
            }

            return ret;
        }

        template <class Space>
        inline static auto multiply(const std::vector<Matrix> &left,
                                    const Gradient<TrialFunction<Space>> &right,
                                    AssemblyContext<LIBMESH_TAG> &ctx) ->
            typename remove_ref_and_const<decltype(grad(right.expr(), ctx))>::type {
            typename remove_ref_and_const<decltype(grad(right.expr(), ctx))>::type ret = grad(right.expr(), ctx);

            const std::size_t n_quad_points = left.size();
            const std::size_t n_functions = ret.size();

            for (std::size_t qp = 0; qp < n_quad_points; ++qp) {
                for (std::size_t i = 0; i < n_functions; ++i) {
                    auto g = ret[i][qp];
                    multiply(left[qp], g, ret[i][qp]);
                }
            }

            return ret;
        }

        template <typename T1, typename T2>
        inline static auto multiply(const std::vector<std::vector<T1>> &left,
                                    const std::vector<std::vector<T2>> &right,
                                    AssemblyContext<LIBMESH_TAG> &ctx)
            -> std::vector<std::vector<decltype(T1() * T2())>> {
            std::vector<std::vector<decltype(T1() * T2())>> ret(left.size());

            return ret;
        }

        template <class Space>
        inline static auto multiply(const std::vector<Matrix> &left,
                                    const TrialFunction<Space> &right,
                                    AssemblyContext<LIBMESH_TAG> &ctx) ->
            typename remove_ref_and_const<decltype(fun(right, ctx))>::type {
            typename remove_ref_and_const<decltype(fun(right, ctx))>::type ret = fun(right, ctx);

            const std::size_t n_quad_points = left.size();
            const std::size_t n_functions = ret.size();

            for (std::size_t qp = 0; qp < n_quad_points; ++qp) {
                for (std::size_t i = 0; i < n_functions; ++i) {
                    auto v = ret[i][qp];
                    multiply(left[qp], v, ret[i][qp]);
                }
            }

            return ret;
        }

        inline static auto multiply(const std::vector<Vector> &left,
                                    const TrialFunction<LibMeshFunctionSpace> &right,
                                    AssemblyContext<LIBMESH_TAG> &ctx) -> std::vector<std::vector<Vector>> {
            const auto &f = fun(right, ctx);
            std::vector<std::vector<Vector>> ret(f.size());

            const std::size_t n_quad_points = left.size();
            const std::size_t n_functions = ret.size();

            for (std::size_t i = 0; i < n_functions; ++i) {
                ret[i].resize(left.size());
                for (std::size_t qp = 0; qp < n_quad_points; ++qp) {
                    // multiply(left[qp], f[i][qp], ret[i][qp]);
                    ret[i][qp] = f[i][qp] * left[qp];
                }
            }

            return ret;
        }

        template <class Left, class Space>
        inline static auto multiply(const Left &left,
                                    const TrialFunction<Space> &right,
                                    AssemblyContext<LIBMESH_TAG> &ctx) ->
            typename remove_ref_and_const<decltype(fun(right, ctx))>::type {
            typename remove_ref_and_const<decltype(fun(right, ctx))>::type ret = fun(right, ctx);

            for (auto &v : ret) {
                for (auto &s : v) {
                    s = left * s;
                }
            }

            return ret;
        }

        template <class Left, class Space>
        inline static auto multiply(const Number<Left> &left,
                                    const TestFunction<Space> &right,
                                    AssemblyContext<LIBMESH_TAG> &ctx) ->
            typename remove_ref_and_const<decltype(fun(right, ctx))>::type {
            typename remove_ref_and_const<decltype(fun(right, ctx))>::type ret = fun(right, ctx);

            const std::size_t n_functions = ret.size();
            const std::size_t n_quad_points = ret[0].size();

            for (std::size_t qp = 0; qp < n_quad_points; ++qp) {
                for (std::size_t i = 0; i < n_functions; ++i) {
                    ret[i][qp] = left * ret[i][qp];
                }
            }

            return ret;
        }

        template <class Space>
        inline static auto multiply(const QValues<double> &left,
                                    const TestFunction<Space> &right,
                                    AssemblyContext<LIBMESH_TAG> &ctx) ->
            typename remove_ref_and_const<decltype(fun(right, ctx))>::type {
            typename remove_ref_and_const<decltype(fun(right, ctx))>::type ret = fun(right, ctx);

            const std::size_t n_functions = ret.size();
            const std::size_t n_quad_points = ret[0].size();

            for (std::size_t qp = 0; qp < n_quad_points; ++qp) {
                for (std::size_t i = 0; i < n_functions; ++i) {
                    ret[i][qp] *= left[qp];
                }
            }

            return ret;
        }

        template <class Space>
        inline static auto multiply(const Matrix &left,
                                    const TestFunction<ProductFunctionSpace<Space>> &right,
                                    AssemblyContext<LIBMESH_TAG> &ctx) -> VectorFunctionType {
            VectorFunctionType ret = fun(right, ctx);

            const std::size_t n_functions = ret.size();
            const std::size_t n_quad_points = ret[0].size();

            for (std::size_t qp = 0; qp < n_quad_points; ++qp) {
                for (std::size_t i = 0; i < n_functions; ++i) {
                    VectorFunctionType::value_type::value_type v = ret[i][qp];
                    ret[i][qp] = left * v;
                }
            }

            return ret;
        }

        // template<class Op>
        inline static auto apply_binary(
            const Interpolate<UVector, TrialFunction<ProductFunctionSpace<LibMeshFunctionSpace>>> &left,
            const GradInterpolate<UVector, TrialFunction<LibMeshFunctionSpace>> &right,
            // const Op &op,
            const Minus &op,
            AssemblyContext<LIBMESH_TAG> &ctx) -> decltype(grad(right.expr().coefficient(), right.expr().fun(), ctx)) {
            auto f = fun(left.coefficient(), left.fun(), ctx);
            auto g = grad(right.expr().coefficient(), right.expr().fun(), ctx);
            auto ret = g;

            std::size_t n_quad_points = g.size();

            for (std::size_t qp = 0; qp < n_quad_points; ++qp) {
                subtract(f[qp], g[qp], ret[qp]);
            }

            return std::move(ret);
        }

        template <class Vector>
        inline static void subtract(const Vector &left, const Vector &right, Vector &result) {
            result = left;
            result -= right;
        }

        template <class Op>
        inline static FunctionType apply_binary(const Number<double> &left,
                                                const TestFunction<LibMeshFunctionSpace> &right,
                                                const Op &op,
                                                AssemblyContext<LIBMESH_TAG> &ctx) {
            FunctionType ret = ctx.test()[right.space_ptr()->subspace_id()]->get_phi();

            const double num = left;

            for (auto &v : ret) {
                for (auto &s : v) {
                    s = op.apply(num, s);
                }
            }

            return ret;
        }

        static auto inner(const std::vector<LMDenseVector> &left,
                          const std::vector<LMDenseVector> &right,
                          const AssemblyContext<LIBMESH_TAG> &ctx) -> std::vector<double> {
            assert(left.size() == right.size());
            const std::size_t n = left.size();

            std::vector<double> ret(n);
            for (std::size_t i = 0; i < n; ++i) {
                ret[i] = utopia::inner(left[i], right[i]);
            }

            return ret;
        }

        template <class Left, class Right>
        static auto inner(const std::vector<Left> &left,
                          const std::vector<Right> &right,
                          const AssemblyContext<LIBMESH_TAG> &ctx) -> std::vector<double> {
            assert(left.size() == right.size());
            const std::size_t n = left.size();

            std::vector<double> ret(n);
            for (std::size_t i = 0; i < n; ++i) {
                ret[i] = utopia::inner(left[i], right[i]);
            }

            return ret;
        }

        template <class Left, class Right>
        static auto inner(const FQValues<Left> &left,
                          const QValues<Right> &right,
                          const AssemblyContext<LIBMESH_TAG> &ctx) -> FQValues<double> {
            assert(left[0].size() == right.size());
            const std::size_t n = left.size();

            FQValues<double> ret(n);
            for (std::size_t i = 0; i < n; ++i) {
                auto n_quad_points = left[i].size();
                ret[i].resize(n_quad_points);
                for (std::size_t qp = 0; qp < n_quad_points; ++qp) {
                    ret[i][qp] = utopia::inner(left[i][qp], right[qp]);
                }
            }

            return ret;
        }

        template <class Left, class Right>
        static auto inner(const QValues<Left> &left,
                          const FQValues<Right> &right,
                          const AssemblyContext<LIBMESH_TAG> &ctx) -> FQValues<double> {
            assert(left.size() == right[0].size());

            const std::size_t n_quad_points = left.size();
            const std::size_t n_funs = right.size();

            std::vector<std::vector<double>> ret(n_funs);

            for (std::size_t i = 0; i < n_funs; ++i) {
                ret[i].resize(n_quad_points);
                for (std::size_t qp = 0; qp < n_quad_points; ++qp) {
                    ret[i][qp] = utopia::inner(left[qp], right[i][qp]);
                }
            }

            return ret;
        }

        template <class Left, int Order, class Right>
        static auto inner(const Tensor<Left, Order> &left,
                          const QValues<Right> &right,
                          const AssemblyContext<LIBMESH_TAG> &ctx) -> QValues<double> {
            auto n_quad_points = right.size();
            QValues<double> ret(n_quad_points);

            for (std::size_t qp = 0; qp < n_quad_points; ++qp) {
                ret[qp] = utopia::inner(left, right[qp]);
            }

            return ret;
        }

        template <class Left, int Order, class Right>
        static auto inner(const Tensor<Left, Order> &left,
                          const FQValues<Right> &right,
                          const AssemblyContext<LIBMESH_TAG> &ctx) -> FQValues<double> {
            auto n_functions = right.size();
            auto n_quad_points = right[0].size();

            FQValues<double> ret(n_functions);

            for (std::size_t i = 0; i < n_functions; ++i) {
                ret[i].resize(right[i].size());
                for (std::size_t qp = 0; qp < n_quad_points; ++qp) {
                    ret[qp] = utopia::inner(left, right[i][qp]);
                }
            }

            return ret;
        }

        template <class Derived, int Order>
        static Derived integrate(Tensor<Derived, Order> &&mat, AssemblyContext<LIBMESH_TAG> &ctx) {
            return std::move(mat.derived());
        }

        template <class Derived, int Order>
        static const Derived &integrate(const Tensor<Derived, Order> &mat, AssemblyContext<LIBMESH_TAG> &ctx) {
            return mat.derived();
        }

        static double integrate(const double val, AssemblyContext<LIBMESH_TAG> &ctx) {
            assert(false);
            return val;
        }

        static Scalar integrate(const QValues<double> &vals, AssemblyContext<LIBMESH_TAG> &ctx) {
            auto &&dx = ctx.dx();
            uint n_quad_points = dx.size();

            assert(vals.size() == dx.size());

            Scalar ret = 0.;
            for (uint qp = 0; qp < n_quad_points; qp++) {
                ret += vals[qp] * dx[qp];
            }

            return ret;
        }

        static Scalar integrate(const ConstantCoefficient<double, 0> &val, AssemblyContext<LIBMESH_TAG> &ctx) {
            auto &&dx = ctx.dx();
            uint n_quad_points = dx.size();

            Scalar ret = 0.;
            for (uint qp = 0; qp < n_quad_points; qp++) {
                ret += dx[qp];
            }

            return ret * val.expr();
        }

        template <class Trial, class Test>
        static Matrix bilinear_form(const Range &trial_range,
                                    const Range &test_range,
                                    Trial &&trial,
                                    Test &&test,
                                    AssemblyContext<LIBMESH_TAG> &ctx) {
            Matrix result;
            init_tensor(result, ctx);

            Write<Matrix> wt(result);
            auto &&dx = ctx.dx();

            uint n_quad_points = dx.size();

            auto s = size(result);

            if (s.n_dims() == 1) {
                s.set_dims(2);
                s.set(1, 1);
            }

            const unsigned int n_test = test_range.extent();
            const unsigned int n_trial = trial_range.extent();

            for (uint i = 0; i < n_trial; i++) {
                for (uint j = 0; j < n_test; j++) {
                    for (uint qp = 0; qp < n_quad_points; qp++) {
                        add(result,
                            test_range.begin() + j,
                            trial_range.begin() + i,
                            utopia::inner(get(trial, qp, i), get(test, qp, j)) * dx[qp]);
                    }
                }
            }

            return result;
        }

        template <class Fun, class Test>
        static Vector linear_form(const Range &test_range, Fun &&fun, Test &&test, AssemblyContext<LIBMESH_TAG> &ctx) {
            Vector result;
            init_tensor(result, ctx);
            Write<Vector> wt(result);

            auto &&dx = ctx.dx();
            uint n_quad_points = dx.size();

            // disp(fun);

            auto n_test = test_range.extent();
            for (uint j = 0; j < n_test; j++) {
                for (uint qp = 0; qp < n_quad_points; qp++) {
                    add(result, test_range.begin() + j, 0, utopia::inner(get(fun, qp, 0), get(test, qp, j)) * dx[qp]);
                }
            }

            return result;
        }
    };
}  // namespace utopia

#endif  // UTOPIA_LIBMESH_FE_BACKEND_HPP
