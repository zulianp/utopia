#ifndef UTOPIA_MOREVB_FEM_HPP
#define UTOPIA_MOREVB_FEM_HPP

#include "utopia_DeviceTensorContraction.hpp"
#include "utopia_DeviceTensorProduct.hpp"
#include "utopia_DiffController.hpp"
#include "utopia_ExtendedFunction.hpp"
#include "utopia_FEFunction.hpp"
#include "utopia_GradInterpolate.hpp"
#include "utopia_LinearElasticityView.hpp"
#include "utopia_PrincipalShapeStressView.hpp"
#include "utopia_PrincipalStrainsView.hpp"
#include "utopia_StrainView.hpp"
#include "utopia_TensorView4.hpp"
#include "utopia_TestFunctions.hpp"
#include "utopia_Tracer.hpp"
#include "utopia_Views.hpp"

namespace utopia {

    template <class FunctionSpace, int Dim = FunctionSpace::Dim>
    class MOREVBFem final : virtual public UnconstrainedExtendedTestFunction<typename FunctionSpace::Matrix,
                                                                             typename FunctionSpace::Vector>,
                            virtual public ConstrainedExtendedTestFunction<typename FunctionSpace::Matrix,
                                                                           typename FunctionSpace::Vector> {
    public:
        using Scalar = typename FunctionSpace::Scalar;
        using SizeType = typename FunctionSpace::SizeType;
        using Vector = typename FunctionSpace::Vector;
        using Matrix = typename FunctionSpace::Matrix;
        using Device = typename FunctionSpace::Device;

        using CSpace = typename FunctionSpace::template Subspace<1>;
        using CElem = typename CSpace::ViewDevice::Elem;

        // FIXME
        using Shape = typename FunctionSpace::Shape;
        using Quadrature = utopia::Quadrature<Shape, 2 * (Shape::Order)>;

        static const int C_NDofs = CSpace::NDofs;
        static const int NQuadPoints = Quadrature::NPoints;

        void read(Input & /*in*/) override {}

        MOREVBFem(FunctionSpace &space) : space_(space) {
            // needed for ML setup
            space_.create_vector(this->_x_eq_values);
            space_.create_vector(this->_eq_constrains_flg);
            space_.create_vector(x_exact_);
            x_exact_.set(0.0);

            this->local_x_ = std::make_shared<Vector>();
            space_.create_local_vector(*this->local_x_);

            this->local_coord_function_ = std::make_shared<Vector>();
            space_.create_local_vector(*this->local_coord_function_);

            this->init_constraints();
            this->init_coord_prod_function(coord_function_);
        }

        inline bool initialize_hessian(Matrix &H, Matrix & /*H_pre*/) const override {
            space_.create_matrix(H);
            return true;
        }

        inline bool update(const Vector & /*x*/) override { return true; }

        bool value(const Vector &x_const, Scalar &val) const override {
            UTOPIA_TRACE_REGION_BEGIN("MembraneFEM::value");

            CSpace C = space_.subspace(0);

            // update local vector x
            space_.global_to_local(x_const, *local_x_);
            auto c_coeff = std::make_shared<Coefficient<CSpace>>(C, local_x_);

            space_.global_to_local(coord_function_, *local_coord_function_);
            auto coord_coeff = std::make_shared<Coefficient<CSpace>>(C, local_coord_function_);

            FEFunction<CSpace> c_fun(c_coeff);
            FEFunction<CSpace> coord_fun(coord_coeff);
            ////////////////////////////////////////////////////////////////////////////
            Quadrature q;

            auto c_val = c_fun.value(q);
            auto coords_val = coord_fun.value(q);

            auto c_grad = c_fun.gradient(q);
            auto differential = C.differential(q);
            // auto c_shape = C.shape(q);

            val = 0.0;

            {
                auto C_view = C.view_device();
                auto c_view = c_val.view_device();
                auto coords_view = coords_val.view_device();

                auto c_grad_view = c_grad.view_device();
                auto differential_view = differential.view_device();

                Device::parallel_reduce(
                    space_.element_range(),
                    UTOPIA_LAMBDA(const SizeType &i) {
                        CElem c_e;
                        C_view.elem(i, c_e);

                        StaticVector<Scalar, NQuadPoints> c;
                        StaticVector<Scalar, NQuadPoints> coords;
                        c_view.get(c_e, c);
                        coords_view.get(c_e, coords);

                        auto c_grad_el = c_grad_view.make(c_e);
                        auto dx = differential_view.make(c_e);

                        Scalar el_energy = 0.0;
                        for (SizeType qp = 0; qp < NQuadPoints; ++qp) {
                            el_energy += 0.5 * inner(c_grad_el[qp], c_grad_el[qp]) * dx(qp);
                            el_energy -= (1. / 8. * std::pow((c[qp] + coords[qp] + 1.), 4.)) * dx(qp);
                        }

                        assert(el_energy == el_energy);
                        return el_energy;
                    },
                    val);
            }

            val = x_const.comm().sum(val);

            UTOPIA_TRACE_REGION_END("MembraneFEM::value");
            return true;
        }

        bool gradient(const Vector &x_const, Vector &g) const override {
            UTOPIA_TRACE_REGION_BEGIN("MembraneFEM::gradient");

            if (empty(g)) {
                space_.create_vector(g);
            } else {
                g.set(0.0);
            }

            CSpace C = space_;

            ///////////////////////////////////////////////////////////////////////////

            // update local vector x
            space_.global_to_local(x_const, *local_x_);
            auto c_coeff = std::make_shared<Coefficient<CSpace>>(C, local_x_);

            space_.global_to_local(coord_function_, *local_coord_function_);
            auto coord_coeff = std::make_shared<Coefficient<CSpace>>(C, local_coord_function_);

            FEFunction<CSpace> c_fun(c_coeff);
            FEFunction<CSpace> coord_fun(coord_coeff);

            Quadrature q;
            auto c_val = c_fun.value(q);
            auto coords_val = coord_fun.value(q);

            auto c_grad = c_fun.gradient(q);
            auto differential = C.differential(q);

            auto c_shape = C.shape(q);
            auto c_grad_shape = C.shape_grad(q);

            {
                auto C_view = C.view_device();
                auto c_view = c_val.view_device();
                auto coords_view = coords_val.view_device();

                auto c_grad_view = c_grad.view_device();

                auto differential_view = differential.view_device();

                auto c_shape_view = c_shape.view_device();
                auto c_grad_shape_view = c_grad_shape.view_device();

                auto g_view = space_.assembly_view_device(g);

                Device::parallel_for(
                    space_.element_range(), UTOPIA_LAMBDA(const SizeType &i) {
                        StaticVector<Scalar, C_NDofs> c_el_vec;
                        c_el_vec.set(0.0);
                        ////////////////////////////////////////////

                        CElem c_e;
                        C_view.elem(i, c_e);
                        StaticVector<Scalar, NQuadPoints> c;
                        StaticVector<Scalar, NQuadPoints> coords;
                        c_view.get(c_e, c);
                        coords_view.get(c_e, coords);

                        auto c_grad_el = c_grad_view.make(c_e);
                        auto dx = differential_view.make(c_e);
                        auto c_grad_shape_el = c_grad_shape_view.make(c_e);
                        auto c_shape_fun_el = c_shape_view.make(c_e);

                        ////////////////////////////////////////////
                        for (SizeType qp = 0; qp < NQuadPoints; ++qp) {
                            for (SizeType j = 0; j < C_NDofs; ++j) {
                                c_el_vec(j) += inner(c_grad_el[qp], c_grad_shape_el(j, qp)) * dx(qp);

                                const Scalar shape_test = c_shape_fun_el(j, qp);
                                c_el_vec(j) -= 0.5 * std::pow((c[qp] + coords[qp] + 1.0), 3) * shape_test * dx(qp);
                            }
                        }

                        C_view.add_vector(c_e, c_el_vec, g_view);
                    });
            }

            // check before boundary conditions
            // // if (this->check_derivatives_) {
            // DiffController<Matrix, Vector> diff_ctrl_;
            // diff_ctrl_.check_grad(*this, x_const, g);
            // // }

            space_.apply_zero_constraints(g);

            UTOPIA_TRACE_REGION_END("MembraneFEM::gradient");
            return true;
        }

        bool hessian(const Vector &x_const, Matrix &H) const override {
            UTOPIA_TRACE_REGION_BEGIN("MembraneFEM::hessian");

            if (empty(H)) {
                // if(use_dense_hessian_) {
                //     H = local_zeros({space_.n_dofs(), space_.n_dofs()}); //FIXME
                // } else {
                space_.create_matrix(H);
                // }
            } else {
                H *= 0.0;
            }

            CSpace C = space_;

            ////////////////////////////////////////////////////////////////////////////
            // update local vector x
            space_.global_to_local(x_const, *local_x_);
            auto c_coeff = std::make_shared<Coefficient<CSpace>>(C, local_x_);

            space_.global_to_local(coord_function_, *local_coord_function_);
            auto coord_coeff = std::make_shared<Coefficient<CSpace>>(C, local_coord_function_);

            FEFunction<CSpace> c_fun(c_coeff);
            FEFunction<CSpace> coord_fun(coord_coeff);

            ////////////////////////////////////////////////////////////////////////////

            Quadrature q;

            auto c_val = c_fun.value(q);
            auto coords_val = coord_fun.value(q);

            auto c_grad = c_fun.gradient(q);
            auto differential = C.differential(q);

            auto c_shape = C.shape(q);
            auto c_grad_shape = C.shape_grad(q);

            {
                auto C_view = C.view_device();

                auto space_view = space_.view_device();

                auto c_view = c_val.view_device();
                auto coords_view = coords_val.view_device();

                auto c_grad_view = c_grad.view_device();
                auto differential_view = differential.view_device();

                auto c_shape_view = c_shape.view_device();
                auto c_grad_shape_view = c_grad_shape.view_device();

                auto H_view = space_.assembly_view_device(H);

                Device::parallel_for(
                    space_.element_range(), UTOPIA_LAMBDA(const SizeType &i) {
                        StaticMatrix<Scalar, C_NDofs, C_NDofs> el_mat;
                        el_mat.set(0.0);

                        ////////////////////////////////////////////
                        CElem c_e;
                        C_view.elem(i, c_e);
                        StaticVector<Scalar, NQuadPoints> c;
                        StaticVector<Scalar, NQuadPoints> coords;
                        c_view.get(c_e, c);
                        coords_view.get(c_e, coords);

                        auto dx = differential_view.make(c_e);
                        auto c_grad_shape_el = c_grad_shape_view.make(c_e);
                        auto c_shape_fun_el = c_shape_view.make(c_e);

                        ////////////////////////////////////////////
                        for (SizeType qp = 0; qp < NQuadPoints; ++qp) {
                            for (SizeType l = 0; l < C_NDofs; ++l) {
                                auto &&c_grad_l = c_grad_shape_el(l, qp);
                                const Scalar c_shape_l = c_shape_fun_el(l, qp);

                                Scalar val1 = 0.5 * 3.0 * std::pow((c[qp] + 1.0), 2) * c_shape_l * dx(qp);
                                el_mat(l, l) -= val1;

                                for (SizeType j = l; j < C_NDofs; ++j) {
                                    Scalar val = inner(c_grad_shape_el(j, qp), c_grad_l) * dx(qp);

                                    val = (l == j) ? (0.5 * val) : val;
                                    el_mat(l, j) += val;
                                    el_mat(j, l) += val;
                                }
                            }
                        }

                        space_view.add_matrix(c_e, el_mat, H_view);
                    });
            }

            space_.apply_constraints(H);

            UTOPIA_TRACE_REGION_END("MembraneFEM::hessian");
            return true;
        }

    private:
        void init_constraints() {
            Vector lb;
            space_.create_vector(lb);

            using Point = typename FunctionSpace::Point;
            using Dev = typename FunctionSpace::Device;
            using Mesh = typename FunctionSpace::Mesh;
            using Elem = typename FunctionSpace::Shape;
            using ElemViewScalar = typename utopia::FunctionSpace<Mesh, 1, Elem>::ViewDevice::Elem;
            static const int NNodes = Elem::NNodes;

            auto C = this->space_;

            auto sampler = utopia::sampler(
                C, UTOPIA_LAMBDA(const Point &coords)->Scalar {
                    auto x = coords[0];
                    auto y = coords[1];

                    const auto PI = 3.14159265358979323846;

                    auto result = std::sin(5. * PI * x) * std::sin(PI * y);
                    result *= std::sin(PI * (1. - x)) * std::sin(PI * (1. - y));

                    return result;
                });

            {
                auto C_view = C.view_device();
                auto sampler_view = sampler.view_device();
                auto x_view = this->space_.assembly_view_device(lb);

                Dev::parallel_for(
                    this->space_.element_range(), UTOPIA_LAMBDA(const SizeType &i) {
                        ElemViewScalar e;
                        C_view.elem(i, e);

                        StaticVector<Scalar, NNodes> s;
                        sampler_view.assemble(e, s);
                        C_view.set_vector(e, s, x_view);
                    });
            }

            this->constraints_ = make_lower_bound_constraints(std::make_shared<Vector>(lb));
        }

    public:
        Vector initial_guess() const override {
            Vector solution;
            space_.create_vector(solution);
            rename("X", solution);

            // TODO:: add initial condition
            solution.set(0.0);
            space_.apply_constraints(solution);

            return solution;
        }

        // not known...
        const Vector &exact_sol() const override {
            std::cout << "MembraneFEM:: exact Solution not known, terminate... \n";
            return x_exact_;
        }

        Scalar min_function_value() const override {
            std::cout << "MembraneFEM:: min_function_value not known, terminate... \n";
            return 0;
        }

        std::string name() const override { return "MembraneFEM"; }

        SizeType dim() const override {
            Vector help;
            space_.create_vector(help);
            return size(help);
        }

        bool exact_sol_known() const override { return false; }

        bool parallel() const override { return true; }

        void init_coord_prod_function(Vector &rhs) {
            using Elem = typename FunctionSpace::Elem;
            static const int NNodes = Elem::NNodes;
            using ElementVector = utopia::StaticVector<Scalar, NNodes>;
            using Mesh = typename FunctionSpace::Mesh;
            using Point = typename Mesh::Point;

            auto coord_function = UTOPIA_LAMBDA(const Point &coords)->Scalar {
                Scalar x = coords[0];
                Scalar y = coords[1];

                return x + y;
            };

            space_.create_vector(rhs);
            Quadrature q;

            Projection<FunctionSpace, Quadrature, decltype(coord_function)> proj(space_, q, coord_function);

            auto proj_view = proj.view_device();

            {
                auto space_view = space_.view_device();
                auto rhs_view = space_.assembly_view_device(rhs);

                Device::parallel_for(
                    space_.element_range(), UTOPIA_LAMBDA(const SizeType &i) {
                        Elem e;
                        space_view.elem(i, e);

                        ElementVector el_vec;
                        el_vec.set(0.0);

                        proj_view.assemble(i, e, el_vec);
                        space_view.add_vector(e, el_vec, rhs_view);
                    });
            }
        }

    private:
        Vector coord_function_;

        FunctionSpace &space_;
        Vector x_exact_;
        std::shared_ptr<Vector> local_x_;
        std::shared_ptr<Vector> local_coord_function_;
    };

}  // namespace utopia

// clean-up macros
#undef UNROLL_FACTOR
#undef U_MIN
#endif  // UTOPIA_MOREVB_FEM_HPP
