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


        bool elastic_energy_in_middle_layer(const Vector &x_const, Scalar &val) const override {
            UTOPIA_TRACE_REGION_BEGIN("IsotropicGenericPhaseField::elastic_energy_in_middle_layer");

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

                        if (centroid[1] > this->non_const_params().bottom_layer_height &&
                            centroid[1] < this->non_const_params().top_layer_height ){
                            //integrate element contribution in middle layer
                            for (SizeType qp = 0; qp < NQuadPoints; ++qp) {
                                Scalar tr = trace(el_strain.strain[qp]);
                                el_energy += elastic_energy(this->params_, c[qp], tr, el_strain.strain[qp]) * dx(qp);
                            }
                         }

                        assert(el_energy == el_energy);
                        return el_energy;
                    },
                    val);
            }

            val = x_const.comm().sum(val);

            assert(val == val);

            UTOPIA_TRACE_REGION_END("IsotropicGenericPhaseField::elastic_energy_in_middle_layer");
            return true;
        }

        bool export_elast_frac_energies(std::string output_path, const Vector &x_const, const Scalar time) const {
            UTOPIA_TRACE_REGION_BEGIN("VolDevGeneric::export_elast_frac_energies");

            static const int total_components = 3; //elastic energy, fracture energy

            using WSpace = typename FunctionSpace::template Subspace<1>;
            using SSpace = typename FunctionSpace::template Subspace<total_components>;
            using SElem = typename SSpace::ViewDevice::Elem;


            Vector w;
            Vector g;
            // Getting displacement subspace
            USpace U;
            this->space_.subspace(1, U);

            WSpace W(this->space_.mesh().clone(1));
            CSpace CC = this->space_.subspace(0);

            /// Creating strain subspace

            // cloning mesh
            auto strain_mesh = this->space_.mesh().clone(total_components);
            assert(strain_mesh->n_components() == total_components);
            // Creating Subspace with cloned mesh

            SSpace S(std::move(strain_mesh));
            assert(S.n_dofs() == CC.n_dofs() * total_components);

            S.create_vector(g);
            W.create_vector(w);

            assert(g.size() == w.size() * total_components);

            ///////////////////////////////////////////////////////////////////////////

            // update local vector x
            this->space_.global_to_local(x_const, *this->local_x_);  // Gets the vector local to the MPI processor
            auto u_coeff = std::make_shared<Coefficient<USpace>>(U, this->local_x_);  // Sets stage for accessing the element node variables
            auto c_coeff = std::make_shared<Coefficient<CSpace>>(CC, this->local_x_);

            // getting FEFunction Space which contains objects for shape function manipulation
            FEFunction<USpace> u_fun(u_coeff);
            FEFunction<CSpace> c_fun(c_coeff);

            {
                ////////////////////////////////////////////////////////////////////////////

                // Quadrature for shape function integration
                Quadrature q;

                // Creating objects for Nodal and Gradient interpolation
                auto u_val = u_fun.value(q);
                auto u_grad = u_fun.gradient(q);
                auto c_val = c_fun.value(q);
                auto c_grad = c_fun.gradient(q);

                // What is thAis for ???
                auto differential = W.differential(q);

                // auto v_grad_shape = U.shape_grad(q);
                auto c_shape = W.shape(q);            // Getting shape functions from FunctionSpace
                auto c_grad_shape = W.shape_grad(q);  // Getting derivative of shape functions from FunctionSpace

                CoefStrain<USpace, Quadrature> strain(u_coeff, q);  // displacement coefficients
                // Strain<USpace, Quadrature> ref_strain_u(U, q); //Test strains (just shape functions gradients for
                // strain)

                auto U_view = U.view_device();
                auto W_view = W.view_device();
                auto S_view = S.view_device();

                auto CC_view = CC.view_device(); //EP

                auto c_view = c_val.view_device();//EP
                auto u_view = u_val.view_device();

                auto strain_view = strain.view_device();
                auto differential_view = differential.view_device();

                // auto v_grad_shape_view = v_grad_shape.view_device();
                auto c_shape_view = c_shape.view_device();  // scalar shape functions
                auto c_grad_view = c_grad.view_device();


                // Preparing the vector for which the Strain function space nows the dimensions (nodes*components), so
                // that we can write on this later
                auto g_view = S.assembly_view_device(g);
                auto w_view = W.assembly_view_device(w);

                // auto ref_strain_u_view = ref_strain_u.view_device();

                Device::parallel_for(
                    this->space_.element_range(), UTOPIA_LAMBDA(const SizeType &i) {
                        Scalar elastic_value, fracture_value, elastic_pos, elastic_neg;
                        StaticVector<Scalar, total_components * C_NDofs> energies_vec;
                        StaticVector<Scalar, C_NDofs> weight_el_vec;
                        StaticMatrix<Scalar, Dim, Dim> stress;  //, strain_p;

                        energies_vec.set(0.0);
                        weight_el_vec.set(0.0);

                        ////////////////////////////////////////////

                        UElem u_e;
                        U_view.elem(i, u_e);
                        auto el_strain =
                            strain_view.make(u_e);  // el_strain.strain[qp] gives matrix of strain at int point

                        SElem s_e;
                        S_view.elem(i, s_e);  // just needed for add_vector into g

                        CElem w_e;
                        W_view.elem(i, w_e);  // getting element for storing wieghts in WSpace

                        CElem cc_e;
                        CC_view.elem(i, cc_e);  // getting element for storing wieghts in CSpace

                        StaticVector<Scalar, NQuadPoints> c;
                        c_view.get(cc_e, c);

                        ////////////////////////////////////////////


                        auto c_grad_el = c_grad_view.make(cc_e);

                        auto dx = differential_view.make(w_e);
                        auto c_shape_fun_el = c_shape_view.make(w_e);  // shape functions (scalar)

                        ////////////////////////////////////////////
                        bool update_elast_tensor = true;
                        Point centroid;
                        w_e.centroid(centroid);
                        this->non_const_params().update(centroid, update_elast_tensor);
                        ////////////////////////////////////////////

                        // loop over all nodes, and for each node, we integrate the strain at the int point weightwd by
                        // the distance to the node (shape function)
                        for (SizeType n = 0; n < C_NDofs; n++) {
                            elastic_value = 0.0;
                            elastic_pos   = 0.0;
                            fracture_value = 0.0;
                            for (SizeType qp = 0; qp < NQuadPoints; ++qp) {
                                auto shape = c_shape_fun_el(n, qp);  // shape function at N and Quadrature point
                                auto weight = dx(qp);                // no need for weights! we want length instead

                                // Calculate strain at quadrature point
                                auto &epsi = el_strain.strain[qp];
                                Scalar tr = trace(el_strain.strain[qp]), elastic_pos_qp;
                                compute_positive_quantities(this->params_, epsi, stress, tr, elastic_pos_qp );

                                 fracture_value += shape * weight *
                                       GenericPhaseFieldFormulation<FunctionSpace, Dim,PFFormulation>::fracture_energy(
                                                this->params_, c[qp], c_grad_el[qp]); //sum fracture energy at integration point

                                elastic_value += shape * weight *
                                        elastic_energy(this->params_, c[qp], tr, el_strain.strain[qp]) ;


                                elastic_pos += shape * weight * PFFormulation::degradation(c[qp], this->params_)*elastic_pos_qp;

                                // getting nodal weight for normalisation
                                weight_el_vec[n] += shape * weight;
                            }

                            // now we need to accumulate the matrix strain into engineering strain vector
                            int offset = C_NDofs;
                            energies_vec[           n ] = elastic_value;
                            energies_vec[offset   + n ] = fracture_value;
                            energies_vec[2*offset + n ] = elastic_pos;
                        }

                        // now adding element contribution to global strain and weight vector
                        S_view.add_vector(s_e, energies_vec, g_view);
                        W_view.add_vector(w_e, weight_el_vec, w_view);
                    });  // end of parallel for

            }  // destruction of view activates MPI Synchronisation

            //            int weight_index = (i - (i % strain_components) ) / strain_components;

            {
                // disp(g.size());
                // disp(w.size());

                // viewing strain vector we just created
                auto strain_and_stress_view = local_view_device(g);
                auto weight_view = local_view_device(w);
                auto r = local_range_device(w);  // range of vector w (using primitivo di utopio)
                parallel_for(
                    r, UTOPIA_LAMBDA(int i) {
                        auto wi = weight_view.get(i);  // extracts vector component
                        for (int k = 0; k < total_components; k++) {
                            int nodal_offset =
                                i * total_components;  // vector g is 2*straincomponents bigger than vector of weights w
                            auto si = strain_and_stress_view.get(
                                nodal_offset + k);  // get k'th strain corresponding to node i with weight i
                            strain_and_stress_view.set(nodal_offset + k,
                                                       si / wi);  // normalise the strain value by the weight wi
                        }
                    });
            }  // incase backed PETSC needs synchronisation (create view in scopes and destroy them when not needed)

            rename("elastic fracture energy", g);
            output_path += "_Energ_" + std::to_string(time) + ".vtr";
            S.write(output_path, g);  // Function space knows how to write

            UTOPIA_TRACE_REGION_END("VolDevGeneric::export_elast_frac_energies");
            return true;
        }

        bool export_strain_and_stress(std::string output_path, const Vector &x_const, const Scalar time) const override{
            UTOPIA_TRACE_REGION_BEGIN("PhaseFieldFracBase::strain");

            static const int strain_components = (Dim - 1) * 3;
            static const int Total_components = strain_components * 2;

            using WSpace = typename FunctionSpace::template Subspace<1>;
            using SSpace = typename FunctionSpace::template Subspace<Total_components>;
            using SElem = typename SSpace::ViewDevice::Elem;

            Vector w;
            Vector g;
            // Getting displacement subspace
            USpace U;
            this->space_.subspace(1, U);

            WSpace C(this->space_.mesh().clone(1));
            CSpace CC = this->space_.subspace(0);

            /// Creating strain subspace

            // cloning mesh
            auto strain_mesh = this->space_.mesh().clone(Total_components);
            assert(strain_mesh->n_components() == Total_components);
            // Creating Subspace with cloned mesh

            SSpace S(std::move(strain_mesh));

            assert(S.n_dofs() == C.n_dofs() * Total_components);

            S.create_vector(g);
            C.create_vector(w);

            assert(g.size() == w.size() * Total_components);

            ///////////////////////////////////////////////////////////////////////////

            // update local vector x
            this->space_.global_to_local(x_const, *this->local_x_);  // Gets the vector local to the MPI processor
            auto u_coeff = std::make_shared<Coefficient<USpace>>(
                U, this->local_x_);  // Sets stage for getting accessing the element node variables
            auto c_coeff = std::make_shared<Coefficient<CSpace>>(CC, this->local_x_);

            // getting FEFunction Space which contains objects for shape function manipulation
            FEFunction<USpace> u_fun(u_coeff);
            FEFunction<CSpace> c_fun(c_coeff);

            {
                ////////////////////////////////////////////////////////////////////////////

                // Quadrature for shape function integration
                Quadrature q;

                // Creating objects for Nodal and Gradient interpolation
                auto u_val = u_fun.value(q);
                auto u_grad = u_fun.gradient(q);
                auto c_val = c_fun.value(q);

                // What is thAis for ???
                auto differential = C.differential(q);

                // auto v_grad_shape = U.shape_grad(q);
                auto c_shape = C.shape(q);            // Getting shape functions from FunctionSpace
                auto c_grad_shape = C.shape_grad(q);  // Getting derivative of shape functions from FunctionSpace

                CoefStrain<USpace, Quadrature> strain(u_coeff, q);  // displacement coefficients
                // Strain<USpace, Quadrature> ref_strain_u(U, q); //Test strains (just shape functions gradients for
                // strain)

                auto U_view = U.view_device();
                auto C_view = C.view_device();
                auto S_view = S.view_device();
                auto CC_view = CC.view_device(); //EP

                auto c_view = c_val.view_device();//EP
                auto u_view = u_val.view_device();

                auto strain_view = strain.view_device();
                auto differential_view = differential.view_device();

                // auto v_grad_shape_view = v_grad_shape.view_device();
                auto c_shape_view = c_shape.view_device();  // scalar shape functions
                // auto c_grad_shape_view = c_grad_shape.view_device();

                // Preparing the vector for which the Strain function space nows the dimensions (nodes*components), so
                // that we can write on this later
                auto g_view = S.assembly_view_device(g);
                auto w_view = C.assembly_view_device(w);

                // auto ref_strain_u_view = ref_strain_u.view_device();

                Device::parallel_for(
                    this->space_.element_range(), UTOPIA_LAMBDA(const SizeType &i) {
                        StaticMatrix<Scalar, Dim, Dim> strain_value, stress_value;
                        StaticVector<Scalar, Total_components * C_NDofs> strain_and_stress_el_vec;
                        StaticVector<Scalar, C_NDofs> weight_el_vec;
                        StaticMatrix<Scalar, Dim, Dim> stress;  //, strain_p;

                        strain_and_stress_el_vec.set(0.0);
                        weight_el_vec.set(0.0);

                        ////////////////////////////////////////////

                        UElem u_e;
                        U_view.elem(i, u_e);
                        auto el_strain =
                            strain_view.make(u_e);  // el_strain.strain[qp] gives matrix of strain at int point

                        SElem s_e;
                        S_view.elem(i, s_e);  // just needed for add_vector into g

                        // auto u_grad_shape_el = v_grad_shape_view.make(u_e);
                        // auto &&u_strain_shape_el = ref_strain_u_view.make(u_e);

                        ////////////////////////////////////////////

                        CElem c_e;
                        C_view.elem(i, c_e);  // getting element for storing wieghts in CSpace

                        CElem cc_e;
                        CC_view.elem(i, cc_e);  // getting element for storing wieghts in CSpace

                        StaticVector<Scalar, NQuadPoints> c;
                        c_view.get(cc_e, c);

                        auto dx = differential_view.make(c_e);
                        auto c_shape_fun_el = c_shape_view.make(c_e);  // shape functions (scalar)

                        ////////////////////////////////////////////
                        bool update_elast_tensor = true;
                        Point centroid;
                        c_e.centroid(centroid);
                        this->non_const_params().update(centroid, update_elast_tensor);
                        ////////////////////////////////////////////

                        // loop over all nodes, and for each node, we integrate the strain at the int point weightwd by
                        // the distance to the node (shape function)
                        for (SizeType n = 0; n < C_NDofs; n++) {
                            strain_value.set(0.0);
                            stress_value.set(0.0);
                            for (SizeType qp = 0; qp < NQuadPoints; ++qp) {
                                auto shape = c_shape_fun_el(n, qp);  // shape function at N and Quadrature point
                                auto weight = dx(qp);                // no need for weights! we want length instead

                                // Calculate strain at quadrature point
                                const auto &epsi = el_strain.strain[qp];

                                // calculate stress at quadrature
                                Scalar tr, gc, elast_energy;
                                compute_stress(this->params_, //need kappa updated
                                               c[qp],
                                               epsi,
                                               stress,  // gets stress (already degraded) at quadrature point
                                               tr,
                                               gc,
                                               elast_energy);

                                strain_value +=
                                    epsi * shape * weight;  // matrix of strains added to existing nodal strain (

                                stress_value += stress * shape * weight;  // Sum stress at integration point

                                // getting nodal weight for normalisation
                                weight_el_vec[n] += shape * weight;
                            }

                            // now we need to accumulate the matrix strain into engineering strain vector
                            int offset = C_NDofs, idx{0};
                            for (int r = 0; r < Dim; ++r) {
                                for (int c = r; c < Dim; c++) {
                                    strain_and_stress_el_vec[idx * offset + n] = stress_value(r, c);
                                    if (strain_components<Total_components)
                                        strain_and_stress_el_vec[(strain_components + idx)*offset + n ] = strain_value(r,c);
                                    idx++;
                                }
                            }
                        }

                        // now adding element contribution to global strain and weight vector
                        S_view.add_vector(s_e, strain_and_stress_el_vec, g_view);
                        C_view.add_vector(c_e, weight_el_vec, w_view);
                    });  // end of parallel for

            }  // destruction of view activates MPI Synchronisation

            //            int weight_index = (i - (i % strain_components) ) / strain_components;

            {
                // disp(g.size());
                // disp(w.size());

                // viewing strain vector we just created
                auto strain_and_stress_view = local_view_device(g);
                auto weight_view = local_view_device(w);
                auto r = local_range_device(w);  // range of vector w (using primitivo di utopio)
                parallel_for(
                    r, UTOPIA_LAMBDA(int i) {
                        auto wi = weight_view.get(i);  // extracts vector component
                        for (int k = 0; k < Total_components; k++) {
                            int nodal_offset =
                                i * Total_components;  // vector g is 2*straincomponents bigger than vector of weights w
                            auto si = strain_and_stress_view.get(
                                nodal_offset + k);  // get k'th strain corresponding to node i with weight i
                            strain_and_stress_view.set(nodal_offset + k,
                                                       si / wi);  // normalise the strain value by the weight wi

//                            if (strain_components != total_components) {
//                                auto sig_i =
//                                    strain_and_stress_view.get(nodal_offset + strain_components +
//                                                               k);  // get stress component which is offset additionally
//                                                                    // in the g vector by the strain components
//                                strain_and_stress_view.set(nodal_offset + strain_components + k, sig_i / wi);
//                            }
                            // assert( std::signbit(si) == std::signbit(sig_i));
                        }
                    });
            }  // incase backed PETSC needs synchronisation (create view in scopes and destroy them when not needed)

            rename("stress and strain", g);
            output_path += "_strainstress_" + std::to_string(time) + ".vtr";
            S.write(output_path, g);  // Function space knows how to write

            UTOPIA_TRACE_REGION_END("PhaseFieldFracBase::strain");
            return true;
        }


        void write_to_file(const std::string &output_path, const Vector &x, const Scalar time) override {
            PhaseFieldFracBase<FunctionSpace, Dim>::write_to_file(output_path, x, time);

            // Post-processing functions
            // And write outputs
            export_strain_and_stress(output_path, x, time); //different for VolDev formulation
            export_elast_frac_energies(output_path,x,time); //diff for VolDev Dormulation
            if (mpi_world_rank() == 0 ) std::cout << "Saving file: " << output_path << std::endl;
        }
    };

}  // namespace utopia

// clean-up macros
#undef UNROLL_FACTOR
#undef U_MIN
#endif
