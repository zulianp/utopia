#ifndef UTOPIA_VOLDEV_GENERIC_PHASE_FIELD_HPP
#define UTOPIA_VOLDEV_GENERIC_PHASE_FIELD_HPP

#include "utopia_CoefStrainView.hpp"
#include "utopia_DeviceTensorContraction.hpp"
#include "utopia_DeviceTensorProduct.hpp"
#include "utopia_DiffController.hpp"
#include "utopia_ExtendedFunction.hpp"
#include "utopia_FEFunction.hpp"
#include "utopia_GenericPhaseFieldFormulation.hpp"
#include "utopia_GradInterpolate.hpp"
#include "utopia_LinearElasticityView.hpp"
#include "utopia_StrainView.hpp"
#include "utopia_TensorView4.hpp"
#include "utopia_Tracer.hpp"
#include "utopia_Views.hpp"
#include "utopia_petsc_NeumannBoundaryConditions.hpp"

#include "utopia_AppBase.hpp"

#define UNROLL_FACTOR 4
#define U_MIN(a, b) ((a) < (b) ? (a) : (b))

namespace utopia {

    template <class FunctionSpace, int Dim = FunctionSpace::Dim, class PFFormulation = AT1 >
    class VolDevGenericPhaseField final : public GenericPhaseFieldFormulation<FunctionSpace, Dim, PFFormulation> {
    public:
        using Scalar = typename FunctionSpace::Scalar;
        using Point = typename FunctionSpace::Point;
        using SizeType = typename FunctionSpace::SizeType;
        using Vector = typename FunctionSpace::Vector;
        using Matrix = typename FunctionSpace::Matrix;
        using Device = typename FunctionSpace::Device;

        using USpace = typename FunctionSpace::template Subspace<Dim>;
        using CSpace = typename FunctionSpace::template Subspace<1>;

        using UElem = typename USpace::ViewDevice::Elem;
        using CElem = typename CSpace::ViewDevice::Elem;
        using MixedElem = typename FunctionSpace::ViewDevice::Elem;
        using Parameters = utopia::PFFracParameters<FunctionSpace>;

        // FIXME
        using Shape = typename FunctionSpace::Shape;
        using Quadrature = utopia::Quadrature<Shape, 2 * (Shape::Order)>;
        static const int NQuadPoints = Quadrature::NPoints;

        static const int C_NDofs = CSpace::NDofs;
        static const int U_NDofs = USpace::NDofs;

        VolDevGenericPhaseField(FunctionSpace &space)
            : GenericPhaseFieldFormulation<FunctionSpace, Dim, PFFormulation>(space) {
            this->params_.fill_in_isotropic_elast_tensor();
        }

        VolDevGenericPhaseField(FunctionSpace &space, const Parameters &params)
            : GenericPhaseFieldFormulation<FunctionSpace, Dim, PFFormulation>(space, params) {
            this->params_.fill_in_isotropic_elast_tensor();
        }



        bool value(const Vector &x_const, Scalar &val) const override {
            UTOPIA_TRACE_REGION_BEGIN("VolDevGenericPhaseField::value");

            USpace U;
            this->space_.subspace(1, U);
            CSpace C = this->space_.subspace(0);

            ///////////////////////////////////////////////////////////////////////////

            // update local vector x
            this->space_.global_to_local(x_const, *this->local_x_);
            auto u_coeff = std::make_shared<Coefficient<USpace>>(U, this->local_x_);
            auto c_coeff = std::make_shared<Coefficient<CSpace>>(C, this->local_x_);

            // udpate local pressure field
            this->space_.global_to_local(this->pressure_field_, *this->local_pressure_field_);
            auto p_coeff = std::make_shared<Coefficient<CSpace>>(C, this->local_pressure_field_);

            // update c_old
            this->space_.global_to_local(this->x_old_, *this->local_c_old_);
            auto c_old_coeff = std::make_shared<Coefficient<CSpace>>(C, this->local_c_old_);

            FEFunction<CSpace> c_old_fun(c_old_coeff);
            FEFunction<CSpace> press_fun(p_coeff);
            FEFunction<CSpace> c_fun(c_coeff);
            FEFunction<USpace> u_fun(u_coeff);
            ////////////////////////////////////////////////////////////////////////////

            Quadrature q;

            auto c_val = c_fun.value(q);
            auto c_old = c_old_fun.value(q);
            auto p_val = press_fun.value(q);

            auto c_grad = c_fun.gradient(q);
            auto u_val = u_fun.value(q);
            auto differential = C.differential(q);

            val = 0.0;

            CoefStrain<USpace, Quadrature> strain(u_fun.coefficient(), q);

            {
                auto U_view = U.view_device();
                auto C_view = C.view_device();

                auto c_view = c_val.view_device();
                auto c_old_view = c_old.view_device();
                auto p_view = p_val.view_device();

                auto c_grad_view = c_grad.view_device();
                auto u_view = u_val.view_device();

                auto strain_view = strain.view_device();
                auto differential_view = differential.view_device();

                Device::parallel_reduce(
                    this->space_.element_range(),
                    UTOPIA_LAMBDA(const SizeType &i) {
                        CElem c_e;
                        C_view.elem(i, c_e);

                        StaticVector<Scalar, NQuadPoints> c;
                        StaticVector<Scalar, NQuadPoints> c_old;
                        StaticVector<Scalar, NQuadPoints> p;
                        c_view.get(c_e, c);
                        c_old_view.get(c_e, c_old);
                        p_view.get(c_e, p);

                        UElem u_e;
                        U_view.elem(i, u_e);
                        auto el_strain = strain_view.make(u_e);
                        auto c_grad_el = c_grad_view.make(c_e);

                        auto dx = differential_view.make(c_e);

                        ////////////////////////////////////////////
                        bool update_elast_tensor = true;
                        Point centroid;
                        c_e.centroid(centroid);
                        this->non_const_params().update(centroid, update_elast_tensor);
                        ////////////////////////////////////////////

                        Scalar el_energy = 0.0;
                        Scalar tr = 0.0;
                        for (SizeType qp = 0; qp < NQuadPoints; ++qp) {
                            el_energy += energy(this->params_, c[qp], c_grad_el[qp], tr, el_strain.strain[qp]) * dx(qp);

                            if (this->params_.use_pressure) {
                                el_energy += PFFormulation::degradation(c[qp], this->params_) *
                                             p[qp] * tr * dx(qp);
                            }

                            if (this->params_.use_penalty_irreversibility) {
                                auto c_cold = c[qp] - c_old[qp];
                                auto c_cold_bracket = c_cold < 0.0 ? c_cold : 0.0;
                                el_energy +=
                                    this->params_.penalty_param_irreversible / 2.0 * c_cold_bracket * c_cold_bracket * dx(qp);

                                if ( PFFormulation::penalise_negative_phase_field_values ){
                                    //std::cout << "Should not be here" << std::endl;
                                    auto c_at_qp = c[qp];
                                    auto c_neg_bracket = c_at_qp < 0.0 ? -c_at_qp : 0.0;
                                    el_energy +=
                                        this->params_.penalty_param_non_neg / 2.0 * c_neg_bracket * c_neg_bracket * dx(qp);
                                }

                            }
                        }

                        assert(el_energy == el_energy);
                        return el_energy;
                    },
                    val);
            }

            val = x_const.comm().sum(val);

            assert(val == val);

            if (!empty(this->force_field_)) {
                // MAYBE -= dot(x_const, this->force_field_);
                val += dot(x_const, this->force_field_);
            }

            this->add_pf_constraints(x_const);   //TAKEN AWAY BUT STILL ACTIVE IN VOL DEV SPLIT -- CHECK

            UTOPIA_TRACE_REGION_END("VolDevGenericPhaseField::value");
            return true;
        }

        bool elastic_energy(const Vector &x_const, Scalar &val) const override {
            UTOPIA_TRACE_REGION_BEGIN("IsotropicGenericPhaseField::elastic_energy");

            USpace U;
            this->space_.subspace(1, U);
            CSpace C = this->space_.subspace(0);

            ///////////////////////////////////////////////////////////////////////////

            // update local vector x
            this->space_.global_to_local(x_const, *this->local_x_);
            auto u_coeff = std::make_shared<Coefficient<USpace>>(U, this->local_x_);
            auto c_coeff = std::make_shared<Coefficient<CSpace>>(C, this->local_x_);

            // udpate local pressure field
            this->space_.global_to_local(this->pressure_field_, *this->local_pressure_field_);
            auto p_coeff = std::make_shared<Coefficient<CSpace>>(C, this->local_pressure_field_);

            // update c_old
            this->space_.global_to_local(this->x_old_, *this->local_c_old_);
            auto c_old_coeff = std::make_shared<Coefficient<CSpace>>(C, this->local_c_old_);

            FEFunction<CSpace> c_old_fun(c_old_coeff);
            FEFunction<CSpace> press_fun(p_coeff);
            FEFunction<CSpace> c_fun(c_coeff);
            FEFunction<USpace> u_fun(u_coeff);
            ////////////////////////////////////////////////////////////////////////////

            Quadrature q;

            auto c_val = c_fun.value(q);
            auto c_old = c_old_fun.value(q);
            auto p_val = press_fun.value(q);

            auto c_grad = c_fun.gradient(q);
            auto u_val = u_fun.value(q);
            auto differential = C.differential(q);

            val = 0.0;

            CoefStrain<USpace, Quadrature> strain(u_coeff, q);

            {
                auto U_view = U.view_device();
                auto C_view = C.view_device();

                auto c_view = c_val.view_device();
                auto c_old_view = c_old.view_device();
                auto p_view = p_val.view_device();

                auto c_grad_view = c_grad.view_device();
                auto u_view = u_val.view_device();

                auto strain_view = strain.view_device();
                auto differential_view = differential.view_device();

                Device::parallel_reduce(
                    this->space_.element_range(),
                    UTOPIA_LAMBDA(const SizeType &i) {
                        CElem c_e;
                        C_view.elem(i, c_e);

                        StaticVector<Scalar, NQuadPoints> c;
                        StaticVector<Scalar, NQuadPoints> c_old;
                        c_view.get(c_e, c);
                        c_old_view.get(c_e, c_old);

                        UElem u_e;
                        U_view.elem(i, u_e);
                        auto el_strain = strain_view.make(u_e);

                        auto dx = differential_view.make(c_e);

                        ////////////////////////////////////////////
                        bool update_elast_tensor = true;
                        Point centroid;
                        c_e.centroid(centroid);
                        this->non_const_params().update(centroid, update_elast_tensor);
                        ////////////////////////////////////////////

                        Scalar el_energy = 0.0;
                        Scalar tr = 0.0;

                        for (SizeType qp = 0; qp < NQuadPoints; ++qp) {

                            el_energy += elastic_energy(this->params_, c[qp], tr, el_strain.strain[qp]) * dx(qp);
                        }

                        assert(el_energy == el_energy);
                        return el_energy;
                    },
                    val);
            }

            val = x_const.comm().sum(val);

            assert(val == val);

            UTOPIA_TRACE_REGION_END("VolDevGenericPhaseField::elastic_energy");
            return true;
        }

        bool fracture_energy(const Vector &x_const, Scalar &val) const override {
            UTOPIA_TRACE_REGION_BEGIN("VolDevGenericPhaseField::fracture_energy");

            USpace U;
            this->space_.subspace(1, U);
            CSpace C = this->space_.subspace(0);

            ///////////////////////////////////////////////////////////////////////////

            // update local vector x
            this->space_.global_to_local(x_const, *this->local_x_);
            auto u_coeff = std::make_shared<Coefficient<USpace>>(U, this->local_x_);
            auto c_coeff = std::make_shared<Coefficient<CSpace>>(C, this->local_x_);

            // udpate local pressure field
            this->space_.global_to_local(this->pressure_field_, *this->local_pressure_field_);
            auto p_coeff = std::make_shared<Coefficient<CSpace>>(C, this->local_pressure_field_);

            // update c_old
            this->space_.global_to_local(this->x_old_, *this->local_c_old_);
            auto c_old_coeff = std::make_shared<Coefficient<CSpace>>(C, this->local_c_old_);

            FEFunction<CSpace> c_old_fun(c_old_coeff);
            FEFunction<CSpace> press_fun(p_coeff);
            FEFunction<CSpace> c_fun(c_coeff);
            FEFunction<USpace> u_fun(u_coeff);
            ////////////////////////////////////////////////////////////////////////////

            Quadrature q;

            auto c_val = c_fun.value(q);
            auto c_old = c_old_fun.value(q);
            auto p_val = press_fun.value(q);

            auto c_grad = c_fun.gradient(q);
            auto u_val = u_fun.value(q);
            auto differential = C.differential(q);

            val = 0.0;

            CoefStrain<USpace, Quadrature> strain(u_coeff, q);

            {
                auto U_view = U.view_device();
                auto C_view = C.view_device();

                auto c_view = c_val.view_device();
                auto c_old_view = c_old.view_device();
                auto p_view = p_val.view_device();

                auto c_grad_view = c_grad.view_device();

                auto differential_view = differential.view_device();

                Device::parallel_reduce(
                    this->space_.element_range(),
                    UTOPIA_LAMBDA(const SizeType &i) {
                        CElem c_e;
                        C_view.elem(i, c_e);

                        StaticVector<Scalar, NQuadPoints> c;
                        StaticVector<Scalar, NQuadPoints> c_old;
                        c_view.get(c_e, c);
                        c_old_view.get(c_e, c_old);

                        auto c_grad_el = c_grad_view.make(c_e);

                        auto dx = differential_view.make(c_e);

                        ////////////////////////////////////////////
                        bool update_elast_tensor = true;
                        Point centroid;
                        c_e.centroid(centroid);
                        this->non_const_params().update(centroid, update_elast_tensor);
                        ////////////////////////////////////////////

                        Scalar el_energy = 0.0;

                        for (SizeType qp = 0; qp < NQuadPoints; ++qp) {
                            el_energy += 
                                GenericPhaseFieldFormulation<FunctionSpace, Dim,PFFormulation>::fracture_energy(
                                    this->params_, c[qp], c_grad_el[qp]) *
                                    dx(qp);
                        }

                        assert(el_energy == el_energy);
                        return el_energy;
                    },
                    val);
            }

            val = x_const.comm().sum(val);

            UTOPIA_TRACE_REGION_END("VolDevGenericPhaseField::fracture_energy");
            return true;
        }

        bool gradient(const Vector &x_const, Vector &g) const override {
            UTOPIA_TRACE_REGION_BEGIN("VolDevGenericPhaseField::gradient");

            if (empty(g)) {
                this->space_.create_vector(g);
            } else {
                g.set(0.0);
            }

            USpace U;
            this->space_.subspace(1, U);
            CSpace C = this->space_.subspace(0);

            ///////////////////////////////////////////////////////////////////////////

            // update local vector x
            this->space_.global_to_local(x_const, *this->local_x_);
            auto u_coeff = std::make_shared<Coefficient<USpace>>(U, this->local_x_);
            auto c_coeff = std::make_shared<Coefficient<CSpace>>(C, this->local_x_);

            // udpate local pressure field
            this->space_.global_to_local(this->pressure_field_, *this->local_pressure_field_);
            auto p_coeff = std::make_shared<Coefficient<CSpace>>(C, this->local_pressure_field_);

            // update c_old
            this->space_.global_to_local(this->x_old_, *this->local_c_old_);
            auto c_old_coeff = std::make_shared<Coefficient<CSpace>>(C, this->local_c_old_);

            FEFunction<CSpace> c_old_fun(c_old_coeff);
            FEFunction<CSpace> press_fun(p_coeff);
            FEFunction<CSpace> c_fun(c_coeff);
            FEFunction<USpace> u_fun(u_coeff);

            ////////////////////////////////////////////////////////////////////////////

            Quadrature q;

            auto c_val = c_fun.value(q);
            auto c_old = c_old_fun.value(q);
            auto press_val = press_fun.value(q);

            auto c_grad = c_fun.gradient(q);
            auto u_val = u_fun.value(q);
            auto differential = C.differential(q);

            // auto v_grad_shape = U.shape_grad(q);
            auto c_shape = C.shape(q);
            auto c_grad_shape = C.shape_grad(q);

            CoefStrain<USpace, Quadrature> strain(u_coeff, q);
            Strain<USpace, Quadrature> ref_strain_u(U, q);

            {
                auto U_view = U.view_device();
                auto C_view = C.view_device();

                auto c_view = c_val.view_device();
                auto c_old_view = c_old.view_device();
                auto p_view = press_val.view_device();

                auto c_grad_view = c_grad.view_device();
                auto u_view = u_val.view_device();

                auto strain_view = strain.view_device();
                auto differential_view = differential.view_device();

                // auto v_grad_shape_view = v_grad_shape.view_device();
                auto c_shape_view = c_shape.view_device();
                auto c_grad_shape_view = c_grad_shape.view_device();

                auto g_view = this->space_.assembly_view_device(g);
                auto ref_strain_u_view = ref_strain_u.view_device();

                Device::parallel_for(
                    this->space_.element_range(), UTOPIA_LAMBDA(const SizeType &i) {
                        StaticMatrix<Scalar, Dim, Dim> stress;  //, strain_p;
                        StaticVector<Scalar, U_NDofs> u_el_vec;
                        StaticVector<Scalar, C_NDofs> c_el_vec;

                        u_el_vec.set(0.0);
                        c_el_vec.set(0.0);

                        ////////////////////////////////////////////

                        UElem u_e;
                        U_view.elem(i, u_e);
                        auto el_strain = strain_view.make(u_e);
                        // auto u_grad_shape_el = v_grad_shape_view.make(u_e);
                        auto &&u_strain_shape_el = ref_strain_u_view.make(u_e);

                        ////////////////////////////////////////////

                        CElem c_e;
                        C_view.elem(i, c_e);
                        StaticVector<Scalar, NQuadPoints> c;
                        c_view.get(c_e, c);

                        StaticVector<Scalar, NQuadPoints> c_old;
                        c_old_view.get(c_e, c_old);

                        StaticVector<Scalar, NQuadPoints> p;
                        p_view.get(c_e, p);

                        auto c_grad_el = c_grad_view.make(c_e);
                        auto dx = differential_view.make(c_e);
                        auto c_grad_shape_el = c_grad_shape_view.make(c_e);
                        auto c_shape_fun_el = c_shape_view.make(c_e);

                        ////////////////////////////////////////////
                        bool update_elast_tensor = true;
                        Point centroid;
                        c_e.centroid(centroid);
                        this->non_const_params().update(centroid, update_elast_tensor);
                        ////////////////////////////////////////////

                        Scalar tr_strain_u, gc, elast;

                        for (SizeType qp = 0; qp < NQuadPoints; ++qp) {

                            compute_stress(this->params_, c[qp], el_strain.strain[qp], stress, tr_strain_u, gc, elast); //also computes gc, tr_strain, and gradient of elastic energy wrt c

                            for (SizeType j = 0; j < U_NDofs; ++j) {
                                auto &&strain_test = u_strain_shape_el(j, qp);
                                u_el_vec(j) += inner(stress, strain_test) * dx(qp);

                                if (this->params_.use_pressure) {
                                    u_el_vec(j) += gc * p[qp] * sum(diag(strain_test)) * dx(qp);
                                }
                            }

                            for (SizeType j = 0; j < C_NDofs; ++j) {
                                const Scalar shape_test = c_shape_fun_el(j, qp);
                                const Scalar frac = GenericPhaseFieldFormulation<FunctionSpace, Dim,PFFormulation>::
                                    grad_fracture_energy_wrt_c(
                                        this->params_, c[qp], c_grad_el[qp], shape_test, c_grad_shape_el(j, qp));

                                c_el_vec(j) += (elast * shape_test + frac) * dx(qp);

                                if (this->params_.use_pressure) {
                                    const Scalar der_c_pres = 
                                        PFFormulation::degradation_deriv(c[qp], this->params_) *p[qp] * tr_strain_u * shape_test;
                                    c_el_vec(j) += der_c_pres * dx(qp);
                                }


                                if (this->params_.use_penalty_irreversibility) {
                                    //c_old[qp] < 0.0 ? 0 : c_old[qp];
                                    auto c_cold = c[qp] - c_old[qp];
                                    auto c_cold_bracket = c_cold < 0.0 ? c_cold : 0.0;
                                    c_el_vec(j) += this->params_.penalty_param_irreversible * c_cold_bracket * shape_test * dx(qp);
                                            //std::cout << "pen: " << this->params_.penalty_param_irreversible << std::endl;

                                    if (PFFormulation::penalise_negative_phase_field_values){
                                        //std::cout << "Should not be here" << std::endl;
                                        auto c_at_qp = c[qp];
                                        auto c_neg_bracket = c_at_qp < 0.0 ? -c_at_qp : 0.0;
                                        c_el_vec(j) += this->params_.penalty_param_non_neg * c_neg_bracket * shape_test * dx(qp);
                                                //std::cout << "pen: " << this->params_.penalty_param_non_neg << std::endl;
                                    }
                                }
                            }
                        }

                        U_view.add_vector(u_e, u_el_vec, g_view);
                        C_view.add_vector(c_e, c_el_vec, g_view);
                    });
            }

            // check before boundary conditions
            if (this->check_derivatives_) {
                this->diff_ctrl_.check_grad(*this, x_const, g);
            }

            if (!empty(this->force_field_)) {
                // MAYBE g -= this->force_field_;
                g += this->force_field_;
            }

            this->space_.apply_zero_constraints(g);

//             fully broken case is treated as Dirichlet BC
            if (this->params_.use_crack_set_irreversibiblity) {
                this->apply_zero_constraints_irreversibiblity(g, x_const);
            }

            UTOPIA_TRACE_REGION_END("VolDevGenericPhaseField::gradient");
            return true;
        }

        bool hessian(const Vector &x_const, Matrix &H) const override {
            UTOPIA_TRACE_REGION_BEGIN("VolDevGenericPhaseField::hessian");

            if (empty(H)) {
                this->space_.create_matrix(H);
            } else {
                H *= 0.0;
            }

            USpace U;
            this->space_.subspace(1, U);
            CSpace C = this->space_.subspace(0);


            ////////////////////////////////////////////////////////////////////////////

            // update local vector x
            this->space_.global_to_local(x_const, *this->local_x_);
            auto u_coeff = std::make_shared<Coefficient<USpace>>(U, this->local_x_);
            auto c_coeff = std::make_shared<Coefficient<CSpace>>(C, this->local_x_);

            // udpate local pressure field
            this->space_.global_to_local(this->pressure_field_, *this->local_pressure_field_);
            auto p_coeff = std::make_shared<Coefficient<CSpace>>(C, this->local_pressure_field_);

            // update c_old
            this->space_.global_to_local(this->x_old_, *this->local_c_old_);
            auto c_old_coeff = std::make_shared<Coefficient<CSpace>>(C, this->local_c_old_);

            FEFunction<CSpace> c_old_fun(c_old_coeff);
            FEFunction<CSpace> press_fun(p_coeff);
            FEFunction<CSpace> c_fun(c_coeff);
            FEFunction<USpace> u_fun(u_coeff);

            ////////////////////////////////////////////////////////////////////////////

            Quadrature q;

            auto c_val = c_fun.value(q);
            auto c_old = c_old_fun.value(q);  //E.P Added old value of damage at quadrature
            auto p_val = press_fun.value(q);

            auto c_grad = c_fun.gradient(q);
            auto u_val = u_fun.value(q);
            auto differential = C.differential(q);

            // auto v_grad_shape = U.shape_grad(q);
            auto c_shape = C.shape(q);
            auto c_grad_shape = C.shape_grad(q);

            // value based
            CoefStrain<USpace, Quadrature> strain(u_coeff, q);

            // reference based
            ShapeStress<USpace, Quadrature> p_stress(U, q, this->params_.mu, this->params_.lambda);
            Strain<USpace, Quadrature> ref_strain_u(U, q);

            {
                UTOPIA_TRACE_REGION_BEGIN("VolDevGenericPhaseField::hessian_local_assembly");

                auto U_view = U.view_device();
                auto C_view = C.view_device();
                auto space_view = this->space_.view_device();

                auto c_view = c_val.view_device();
                auto c_old_view = c_old.view_device();      //E.P getting view device of old damage within scope
                auto p_view = p_val.view_device();

                auto c_grad_view = c_grad.view_device();
                auto u_view = u_val.view_device();

                auto strain_view = strain.view_device();
                auto differential_view = differential.view_device();

                // auto v_grad_shape_view = v_grad_shape.view_device();
                auto c_shape_view = c_shape.view_device();
                auto c_grad_shape_view = c_grad_shape.view_device();

                // FIXME
                auto p_stress_view = p_stress.view_device();

                auto H_view = this->space_.assembly_view_device(H);
                auto ref_strain_u_view = ref_strain_u.view_device();

                Device::parallel_for(
                    this->space_.element_range(), UTOPIA_LAMBDA(const SizeType &i) {
                        // StaticMatrix<Scalar, Dim, Dim> strain_n, strain_p;
                        StaticMatrix<Scalar, U_NDofs + C_NDofs, U_NDofs + C_NDofs> el_mat;
                        StaticMatrix<Scalar, Dim, Dim> stress_positive;

                        MixedElem e;
                        space_view.elem(i, e);
                        el_mat.set(0.0);

                        ////////////////////////////////////////////
                        UElem u_e;
                        U_view.elem(i, u_e);
                        auto el_strain = strain_view.make(u_e);
                        auto &&u_strain_shape_el = ref_strain_u_view.make(u_e);

                        ////////////////////////////////////////////
                        CElem c_e;
                        C_view.elem(i, c_e);
                        StaticVector<Scalar, NQuadPoints> c;
                        StaticVector<Scalar, NQuadPoints> p;
                        c_view.get(c_e, c);
                        p_view.get(c_e, p);

                        StaticVector<Scalar, NQuadPoints> c_old;    //E.P Added c_old vector for penalty irreversability
                        c_old_view.get(c_e, c_old);

                        auto dx = differential_view.make(c_e);
                        auto c_grad_shape_el = c_grad_shape_view.make(c_e);
                        auto c_shape_fun_el = c_shape_view.make(c_e);

                        ////////////////////////////////////////////
                        bool update_elast_tensor = true;
                        Point centroid;
                        c_e.centroid(centroid);
                        this->non_const_params().update(centroid, update_elast_tensor);
                        //Getting new material parameter values

                        ////////////////////////////////////////////

                        Scalar tr_strain_u, eep;

                        for (SizeType qp = 0; qp < NQuadPoints; ++qp) {
                            compute_positive_quantities(
                                this->params_, el_strain.strain[qp], stress_positive, tr_strain_u, eep); //Positive Strain energy  calculated here

                            //const Scalar eep = elastic_energy(this->params_, c[qp], tr_strain_u, el_strain.strain[qp]);
                            //const Scalar eep_fix = strain_energy(this->params_, tr_strain_u, el_strain.strain[qp]);

                            // pragma GCCunroll(C_NDofs)
                            for (SizeType l = 0; l < C_NDofs; ++l) {
                                const Scalar c_shape_l = c_shape_fun_el(l, qp);
                                auto &&c_grad_l = c_grad_shape_el(l, qp);

                                // SYMMETRIC VERSION
                                for (SizeType j = l; j < C_NDofs; ++j) {
                                    const Scalar c_shape_j_l_prod = c_shape_fun_el(j, qp) * c_shape_l;

                                    Scalar val = bilinear_cc(this->params_,
                                                             c[qp],
                                                             eep,
                                                             // c_shape_j,
                                                             // c_shape_l,
                                                             c_shape_j_l_prod,
                                                             c_grad_shape_el(j, qp),
                                                             c_grad_l) *
                                                 dx(qp);

                                    if (this->params_.use_pressure) {
                                        val += PFFormulation::degradation_deriv2(c[qp], this->params_) * p[qp] * tr_strain_u *
                                               c_shape_j_l_prod * dx(qp);
                                    }


                                    if (this->params_.use_penalty_irreversibility) {
                                        auto c_cold = c[qp] - c_old[qp];
                                        auto c_heaviside = c_cold <= 0.0 ? 1.0 : 0.0;
                                        val += c_heaviside* this->params_.penalty_param_irreversible *
                                               c_shape_j_l_prod * dx(qp);


                                        if (PFFormulation::penalise_negative_phase_field_values ){
                                            //std::cout << "Should not be here" << std::endl;
                                            auto c_at_qp = c[qp];
                                            auto c_neg_bracket = c_at_qp <= 0.0 ? 1.0 : 0.0;
                                            val += this->params_.penalty_param_non_neg * c_neg_bracket * 
                                                    c_shape_j_l_prod * dx(qp);
                                        }
                                    }

                                    val = (l == j) ? (0.5 * val) : val;

                                    el_mat(l, j) += val;
                                    el_mat(j, l) += val;
                                }
                            }

                            for (SizeType l = 0; l < U_NDofs; ++l) {
                                auto &&u_strain_shape_l = u_strain_shape_el(l, qp);

                                for (SizeType j = l; j < U_NDofs; ++j) {
                                // Varying stress tensor
                                     Scalar val = bilinear_uu(this->params_,
                                                             c[qp],
                                                             el_strain.strain[qp],
                                                             u_strain_shape_el(j, qp),
                                                             u_strain_shape_l) *
                                                   dx(qp);


                                    val = (l == j) ? (0.5 * val) : val;
                                    el_mat(C_NDofs + l, C_NDofs + j) += val;
                                    el_mat(C_NDofs + j, C_NDofs + l) += val;
                                }
                            }

                            //////////////////////////////////////////////////////////////////////////////////////////////////////

                            if (this->params_.turn_off_uc_coupling == false ||
                                this->params_.turn_off_cu_coupling == false) {
                                // #pragma clang loop unroll_count(U_MIN(C_NDofs, UNROLL_FACTOR))
                                // #pragma GCC unroll U_MIN(C_NDofs, UNROLL_FACTOR)

                                //compute_stress(this->params_, tr_strain_u, el_strain.strain[qp], stress);
                                for (SizeType c_i = 0; c_i < C_NDofs; ++c_i) {
                                    // CHANGE (pre-compute/store shape fun)
                                    const Scalar c_shape_i = c_shape_fun_el(c_i, qp);
                                    for (SizeType u_i = 0; u_i < U_NDofs; ++u_i) {
                                        auto &&strain_shape = u_strain_shape_el(u_i, qp);

                                        Scalar val =
                                            bilinear_uc(this->params_, c[qp], stress_positive, strain_shape, c_shape_i) * dx(qp);

                                        if (this->params_.use_pressure) {
                                            const Scalar tr_strain_shape = sum(diag(strain_shape));
                                            val += PFFormulation::degradation_deriv(            //SHOULD THIS BE deriv2 !!?? BUG? Check
                                                       c[qp], this->params_) *
                                                   p[qp] * tr_strain_shape * c_shape_i * dx(qp);
                                        }

                                        // not symetric, but more numerically stable
                                        if (this->params_.turn_off_cu_coupling == false) {
                                            el_mat(c_i, C_NDofs + u_i) += val;
                                        }

                                        if (this->params_.turn_off_uc_coupling == false) {
                                            el_mat(C_NDofs + u_i, c_i) += val;
                                        }
                                    }
                                }
                            }
                        }

                        space_view.add_matrix(e, el_mat, H_view);
                    });

                UTOPIA_TRACE_REGION_END("VolDevGenericPhaseField::hessian_local_assembly");
            }

            // check before boundary conditions
            if (this->check_derivatives_) {
                this->diff_ctrl_.check_hessian(*this, x_const, H);
            }

            this->space_.apply_constraints(H);

            if (this->params_.use_crack_set_irreversibiblity) {
                this->apply_zero_constraints_irreversibiblity(H, x_const);
            }

            UTOPIA_TRACE_REGION_END("VolDevGenericPhaseField::hessian");
            return true;
        }

        //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

        template <class Grad>
        UTOPIA_INLINE_FUNCTION static Scalar bilinear_uu(const Parameters &params,
                                                         const Scalar &phase_field_value,
                                                         const Grad &strain,
                                                         const Grad &strain_trial,
                                                         const Grad &strain_test) {
            const Scalar strain0tr = trace(strain);

            Tensor4th<Scalar, Dim, Dim, Dim, Dim> Jacobian_neg, Jacobian_mult;

            if (strain0tr < 0) {
                Jacobian_neg = params.kappa * params.I4sym;
            }

            Scalar gc = PFFormulation::degradation(phase_field_value, params );

            // would be nicer, if this works without 4th order tensor...
            Jacobian_mult = params.elast_tensor - Jacobian_neg;
            Jacobian_mult = (gc * Jacobian_mult) + Jacobian_neg;

            Scalar val = inner(strain_trial, contraction(Jacobian_mult, strain_test));

            return val;
        }

        template <class GradShape>
        UTOPIA_INLINE_FUNCTION static Scalar bilinear_cc(const Parameters &params,
                                                         const Scalar &phase_field_value,
                                                         const Scalar &strain_energy,
                                                         // const Scalar &shape_trial,
                                                         // const Scalar &shape_test,
                                                         const Scalar &shape_prod,
                                                         const GradShape &grad_trial,
                                                         const GradShape &grad_test) {
            return GenericPhaseFieldFormulation<FunctionSpace, Dim,PFFormulation>::diffusion_c(
                        params, grad_trial, grad_test) +
                   GenericPhaseFieldFormulation<FunctionSpace,Dim, PFFormulation>::reaction_c(
                        params, phase_field_value, shape_prod) +
                   elastic_deriv_cc(params, phase_field_value, strain_energy, shape_prod);
        }


        template <class Stress, class FullStrain>
        UTOPIA_INLINE_FUNCTION static Scalar bilinear_uc(const Parameters &params,
                                                         const Scalar &phase_field_value,
                                                         const Stress &stress,
                                                         const FullStrain &full_strain,
                                                         const Scalar &c_trial_fun) {
            return PFFormulation::degradation_deriv(phase_field_value, params) * c_trial_fun * inner(stress, full_strain);
        }


        UTOPIA_INLINE_FUNCTION static Scalar elastic_deriv_cc(const Parameters &params,
                                                              const Scalar &phase_field_value,
                                                              const Scalar &strain_energy,
                                                              // const Scalar &trial,
                                                              // const Scalar &test
                                                              const Scalar &shape_prod) {
            const Scalar dcc = (1.0 - params.regularization) * PFFormulation::degradation_deriv2(phase_field_value, params);

            return dcc * shape_prod * strain_energy;
        }


        template <class Strain, class Stress>
        UTOPIA_INLINE_FUNCTION static void compute_stress(const Parameters &params,
                                                          const Scalar &phase_field_value,
                                                          const Strain &strain,
                                                          Stress &stress,
                                                          Scalar &tr,
                                                          Scalar &gc,
                                                          Scalar &elast_energy) {
            tr = trace(strain);
            const Strain strain_dev = strain - (((1. / Dim) * tr) * device::identity<Scalar>());

            const Scalar tr_negative = device::min(tr, 0.0);
            const Scalar tr_positive = tr - tr_negative;

            stress = ((params.kappa * tr_positive) * device::identity<Scalar>());
            stress += ((2.0 * params.mu) * strain_dev);

            gc = PFFormulation::degradation(phase_field_value, params);
            stress = (gc * stress) + ((params.kappa * tr_negative) * device::identity<Scalar>());

            // energy positive
            const Scalar energy_positive =
                (0.5 * params.kappa * tr_positive * tr_positive) + (params.mu * inner(strain_dev, strain_dev));

            elast_energy =
                PFFormulation::degradation_deriv(phase_field_value, params) *
                energy_positive;
        }


        template <class Strain, class Stress>
        UTOPIA_INLINE_FUNCTION static void compute_positive_quantities(const Parameters &params,
                                                                       const Strain &strain,
                                                                       Stress &stress_positive,
                                                                       Scalar &tr,
                                                                       Scalar &energy_positive) {
            tr = trace(strain);
            const Strain strain_dev = strain - (((1. / Dim) * tr) * device::identity<Scalar>());

            const Scalar tr_negative = device::min(tr, 0.0);
            const Scalar tr_positive = tr - tr_negative;

            stress_positive = ((params.kappa * tr_positive) * device::identity<Scalar>());
            stress_positive += ((2.0 * params.mu) * strain_dev);

            energy_positive =
                (0.5 * params.kappa * tr_positive * tr_positive) + (params.mu * inner(strain_dev, strain_dev));
        }


        template <class Grad, class Strain>
        UTOPIA_INLINE_FUNCTION static Scalar energy(const Parameters &params,
                                                    // c
                                                    const Scalar &phase_field_value,
                                                    const Grad &phase_field_grad,
                                                    // u
                                                    Scalar &trace,
                                                    const Strain &strain) {
            return GenericPhaseFieldFormulation<FunctionSpace, Dim, PFFormulation>::fracture_energy(
                       params, phase_field_value, phase_field_grad) +
                   elastic_energy(params, phase_field_value, trace, strain);
        }


        template <class Strain>
        UTOPIA_INLINE_FUNCTION static Scalar elastic_energy(const Parameters &params,
                                                            const Scalar &phase_field_value,
                                                            Scalar &tr,
                                                            const Strain &strain) {
            tr = trace(strain);
            Strain strain_dev = strain - ((1. / Dim) * tr * device::identity<Scalar>());

            const Scalar tr_negative = device::min(tr, 0.0);
            const Scalar tr_positive = tr - tr_negative;

            const Strain strain_dev2 = transpose(strain_dev) * strain_dev;
            Scalar tr2dv = trace(strain_dev2);

            // energy positive
            const Scalar energy_positive = (0.5 * params.kappa * tr_positive * tr_positive) + (params.mu * tr2dv);
            const Scalar energy_negative = (0.5 * params.kappa * tr_negative * tr_negative);

            const Scalar energy =
                (PFFormulation::degradation(phase_field_value, params) *
                 energy_positive) +
                energy_negative;

            return energy;
        }

        void write_to_file(const std::string &output_path, const Vector &x, const Scalar time) override {
            PhaseFieldFracBase<FunctionSpace, Dim>::write_to_file(output_path, x, time);

            // Post-processing functions
            // And write outputs
            this->export_strain_and_stress(output_path, x, time);
            if (mpi_world_rank() == 0 ) std::cout << "Saving file: " << output_path << std::endl;
        }
    };

}  // namespace utopia

// clean-up macros
#undef UNROLL_FACTOR
#undef U_MIN
#endif
