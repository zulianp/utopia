#ifndef UTOPIA_PHASE_FIELD_HPP
#define UTOPIA_PHASE_FIELD_HPP

#include "utopia_CoefStrainView.hpp"
#include "utopia_DeviceTensorContraction.hpp"
#include "utopia_DiffController.hpp"
#include "utopia_ExtendedFunction.hpp"
#include "utopia_FEFunction.hpp"
#include "utopia_GradInterpolate.hpp"
#include "utopia_LinearElasticityView.hpp"
#include "utopia_PhaseFieldBase.hpp"
#include "utopia_StrainView.hpp"
#include "utopia_TensorView4.hpp"
#include "utopia_Tracer.hpp"
#include "utopia_Views.hpp"
#include "utopia_petsc_NeumannBoundaryConditions.hpp"

#include "utopia_AppBase.hpp"

#define UNROLL_FACTOR 4
#define U_MIN(a, b) ((a) < (b) ? (a) : (b))

namespace utopia {

    template <class FunctionSpace, int Dim = FunctionSpace::Dim>
    class IsotropicPhaseFieldForBrittleFractures final : public PhaseFieldFracBase<FunctionSpace, Dim> {
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

        IsotropicPhaseFieldForBrittleFractures(FunctionSpace &space) : PhaseFieldFracBase<FunctionSpace, Dim>(space) {}

        IsotropicPhaseFieldForBrittleFractures(FunctionSpace &space, const Parameters &params)
            : PhaseFieldFracBase<FunctionSpace, Dim>(space, params) {}

        /*
        bool export_strain(std::string output_path , const Vector &x_const, const Scalar time) const {
            UTOPIA_TRACE_REGION_BEGIN("IsotropicPhaseFieldForBrittleFractures::strain");

            static const int strain_components = (Dim - 1) * 3;

            using WSpace = typename FunctionSpace::template Subspace<1>;
            using SSpace = typename FunctionSpace::template Subspace<strain_components>;
            using SElem = typename SSpace::ViewDevice::Elem;

            Vector w;
            Vector g;
            // Getting displacement subspace
            USpace U;
            this->space_.subspace(1, U);

            WSpace C(this->space_.mesh().clone(1));

            /// Creating strain subspace

            // cloning mesh
            auto strain_mesh = this->space_.mesh().clone(strain_components);
            assert(strain_mesh->n_components() == strain_components);
            // Creating Subspace with cloned mesh

            SSpace S(std::move(strain_mesh));

            assert(S.n_dofs() == C.n_dofs() * strain_components);

            S.create_vector(g);
            C.create_vector(w);

            assert(g.size() == w.size() * strain_components);

            ///////////////////////////////////////////////////////////////////////////

            // update local vector x
            this->space_.global_to_local(x_const, *this->local_x_);  // Gets the vector local to the MPI processor
            auto u_coeff = std::make_shared<Coefficient<USpace>>(
                U, this->local_x_);  // Sets stage for getting accessing the element node variables

            // getting FEFunction Space which contains objects for shape function manipulation
            FEFunction<USpace> u_fun(u_coeff);

            {
                ////////////////////////////////////////////////////////////////////////////

                // Quadrature for shape function integration
                Quadrature q;

                // Creating objects for Nodal and Gradient interpolation
                auto u_val = u_fun.value(q);
                auto u_grad = u_fun.gradient(q);

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
                        StaticMatrix<Scalar, Dim, Dim> strain_value;  //, strain_p;
                        StaticVector<Scalar, strain_components * C_NDofs> strain_el_vec;
                        StaticVector<Scalar, C_NDofs> weight_el_vec;

                        strain_el_vec.set(0.0);
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

                        auto dx = differential_view.make(c_e);
                        auto c_shape_fun_el = c_shape_view.make(c_e);  // shape functions (scalar)

                        // loop over all nodes, and for each node, we integrate the strain at the int point weightwd by
                        // the distance to the node (shape function)
                        for (SizeType n = 0; n < C_NDofs; n++) {
                            strain_value.set(0.0);
                            for (SizeType qp = 0; qp < NQuadPoints; ++qp) {
                                // const Scalar tr_strain_u = trace(el_strain.strain[qp]);

                                auto shape = c_shape_fun_el(n, qp);  // shape function at N and Quadrature point
                                auto weight = dx(qp);                // no need for weights! we want length instead

                                auto &epsi = el_strain.strain[qp];  // strain at integration point
                                strain_value +=
                                    epsi * shape * weight;  // matrix of strains added to existing nodal strain (

                                // getting nodal weight for normalisation
                                weight_el_vec[n] += shape * weight;
                            }

                            // now we need to accumulate the matrix strain into engineering strain vector
                            int offset = n * strain_components;
                            for (int r = 0; r < Dim; ++r) {
                                for (int c = r; c < Dim; c++) {
                                    strain_el_vec[offset++] = strain_value(r, c);
                                }
                            }
                        }

                        // now adding element contribution to global strain and weight vector
                        S_view.add_vector(s_e, strain_el_vec, g_view);
                        C_view.add_vector(c_e, weight_el_vec, w_view);
                    });  // end of parallel for

            }  // destruction of view activates MPI Synchronisation

            //            int weight_index = (i - (i % strain_components) ) / strain_components;

            {
                // disp(g.size());
                // disp(w.size());

                // viewing strain vector we just created
                auto strain_view = local_view_device(g);
                auto weight_view = local_view_device(w);
                auto r = local_range_device(w);  // range of vector w (using primitivo di utopio)
                parallel_for(
                    r, UTOPIA_LAMBDA(int i) {
                        auto wi = weight_view.get(i);  // extracts vector component
                        for (int k = 0; k < strain_components; k++) {
                            int idx = i * strain_components + k;
                            auto si = strain_view.get(idx);
                            strain_view.set(idx, si / wi);
                        }
                    });
            }  // incase backed PETSC needs synchronisation (create view in scopes and destroy them when not needed)

            rename("strain", g);
            output_path += "_strain_" + std::to_string(time) + ".vtr";
            S.write(output_path, g);  // Function space knows how to write

            UTOPIA_TRACE_REGION_END("IsotropicPhaseFieldForBrittleFractures::strain");
            return true;
        }
        */

        bool value(const Vector &x_const, Scalar &val) const override {
            UTOPIA_TRACE_REGION_BEGIN("IsotropicPhaseFieldForBrittleFractures::value");

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
                        bool update_elast_tensor = false;
                        Point centroid;
                        c_e.centroid(centroid);
                        this->non_const_params().update(centroid, update_elast_tensor);
                        ////////////////////////////////////////////

                        Scalar el_energy = 0.0;
                        for (SizeType qp = 0; qp < NQuadPoints; ++qp) {
                            Scalar tr = trace(el_strain.strain[qp]);
                            if (this->params_.use_pressure) {
                                el_energy += PhaseFieldFracBase<FunctionSpace, Dim>::quadratic_degradation(
                                                 this->params_, c[qp]) *
                                             p[qp] * tr * dx(qp);
                            }

                            el_energy += energy(this->params_, c[qp], c_grad_el[qp], tr, el_strain.strain[qp]) * dx(qp);

                            if (this->params_.use_penalty_irreversibility) {
                                // if (mpi_world_rank()==0) utopia::out() << "Penalty: " << this->params_.penalty_param
                                // << std::endl;

                                auto c_cold = c[qp] - c_old[qp];
                                auto c_cold_bracket = c_cold < 0.0 ? c_cold : 0.0;
                                el_energy += this->params_.penalty_param_irreversible / 2.0 * c_cold_bracket *
                                             c_cold_bracket * dx(qp);
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
                val -= dot(x_const, this->force_field_);
            }

            if (this->params_.use_crack_set_irreversibiblity) this->add_pf_constraints(x_const);

            UTOPIA_TRACE_REGION_END("IsotropicPhaseFieldForBrittleFractures::value");
            return true;
        }

        bool elastic_energy(const Vector &x_const, Scalar &val) const override {
            UTOPIA_TRACE_REGION_BEGIN("IsotropicPhaseFieldForBrittleFractures::elastic_energy");

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
                        bool update_elast_tensor = false;
                        Point centroid;
                        c_e.centroid(centroid);
                        this->non_const_params().update(centroid, false);
                        ////////////////////////////////////////////

                        Scalar el_energy = 0.0;

                        for (SizeType qp = 0; qp < NQuadPoints; ++qp) {
                            Scalar tr = trace(el_strain.strain[qp]);

                            el_energy += elastic_energy(this->params_, c[qp], tr, el_strain.strain[qp]) * dx(qp);
                        }

                        assert(el_energy == el_energy);
                        return el_energy;
                    },
                    val);
            }

            val = x_const.comm().sum(val);

            assert(val == val);

            UTOPIA_TRACE_REGION_END("IsotropicPhaseFieldForBrittleFractures::elastic_energy");
            return true;
        }

        bool fracture_energy(const Vector &x_const, Scalar &val) const override {
            UTOPIA_TRACE_REGION_BEGIN("IsotropicPhaseFieldForBrittleFractures::fracture_energy");

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
                        bool update_elast_tensor = false;
                        Point centroid;
                        c_e.centroid(centroid);
                        this->non_const_params().update(centroid, update_elast_tensor);
                        ////////////////////////////////////////////

                        Scalar el_energy = 0.0;

                        for (SizeType qp = 0; qp < NQuadPoints; ++qp) {
                            el_energy += PhaseFieldFracBase<FunctionSpace, Dim>::fracture_energy(
                                             this->params_, c[qp], c_grad_el[qp]) *
                                         dx(qp);
                        }

                        assert(el_energy == el_energy);
                        return el_energy;
                    },
                    val);
            }

            val = x_const.comm().sum(val);

            UTOPIA_TRACE_REGION_END("IsotropicPhaseFieldForBrittleFractures::fracture_energy");
            return true;
        }

        bool gradient(const Vector &x_const, Vector &g) const override {
            UTOPIA_TRACE_REGION_BEGIN("IsotropicPhaseFieldForBrittleFractures::gradient");

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
                        bool update_elast_tensor = false;
                        Point centroid;
                        c_e.centroid(centroid);
                        this->non_const_params().update(centroid, update_elast_tensor);
                        ////////////////////////////////////////////

                        for (SizeType qp = 0; qp < NQuadPoints; ++qp) {
                            const Scalar tr_strain_u = trace(el_strain.strain[qp]);

                            compute_stress(this->params_, tr_strain_u, el_strain.strain[qp], stress);
                            Scalar gc_qp =
                                PhaseFieldFracBase<FunctionSpace, Dim>::quadratic_degradation(this->params_, c[qp]);
                            stress =
                                (gc_qp * (1.0 - this->params_.regularization) + this->params_.regularization) * stress;

                            for (SizeType j = 0; j < U_NDofs; ++j) {
                                auto &&strain_test = u_strain_shape_el(j, qp);
                                u_el_vec(j) += inner(stress, strain_test) * dx(qp);

                                if (this->params_.use_pressure) {
                                    u_el_vec(j) += gc_qp * p[qp] * sum(diag(strain_test)) * dx(qp);
                                }
                            }

                            const Scalar elast =
                                grad_elastic_energy_wrt_c(this->params_, c[qp], tr_strain_u, el_strain.strain[qp]);

                            for (SizeType j = 0; j < C_NDofs; ++j) {
                                const Scalar shape_test = c_shape_fun_el(j, qp);
                                const Scalar frac = PhaseFieldFracBase<FunctionSpace, Dim>::grad_fracture_energy_wrt_c(
                                    this->params_, c[qp], c_grad_el[qp], shape_test, c_grad_shape_el(j, qp));

                                if (this->params_.use_pressure) {
                                    const Scalar der_c_pres =
                                        PhaseFieldFracBase<FunctionSpace, Dim>::quadratic_degradation_deriv(
                                            this->params_, c[qp]) *
                                        p[qp] * tr_strain_u * shape_test;
                                    c_el_vec(j) += der_c_pres * dx(qp);
                                }

                                c_el_vec(j) += (elast * shape_test + frac) * dx(qp);

                                if (this->params_.use_penalty_irreversibility) {
                                    auto c_cold = c[qp] - c_old[qp];
                                    auto c_cold_bracket = c_cold < 0.0 ? c_cold : 0.0;
                                    c_el_vec(j) +=
                                        this->params_.penalty_param_irreversible * c_cold_bracket * shape_test * dx(qp);
                                }
                            }
                        }

                        //                         printf("\n------ u gradient ------\n");
                        //                         for (int i = 0; i < U_NDofs; i++) {
                        //                            double val = u_el_vec(i);
                        //                            if (std::abs(val) < 1e-4) {
                        //                                val = 0;
                        //                            }
                        //                            printf("%.4g, ", val);
                        //                         }
                        //                         printf("\n------- c gradient ------\n");
                        //                         for (int j = 0; j < C_NDofs; j++) {
                        //                                 double val = c_el_vec(j);
                        //                                 if (std::abs(val) < 1e-4) {
                        //                                     val = 0;
                        //                                 }
                        //                                 printf("%.4g, ", val);
                        //                             }
                        //                         printf("\n");

                        U_view.add_vector(u_e, u_el_vec, g_view);
                        C_view.add_vector(c_e, c_el_vec, g_view);
                    });
            }

            if (!empty(this->force_field_)) {
                g -= this->force_field_;
                //                g += this->force_field_;
            }

            // check before boundary conditions
            if (this->check_derivatives_) {
                this->diff_ctrl_.check_grad(*this, x_const, g);
            }

            this->space_.apply_zero_constraints(g);

            // fully broken case is treated as Dirichlet BC
            if (this->params_.use_crack_set_irreversibiblity) {
                this->apply_zero_constraints_irreversibiblity(g, x_const);
            }

            UTOPIA_TRACE_REGION_END("IsotropicPhaseFieldForBrittleFractures::gradient");
            return true;
        }

        bool hessian(const Vector &x_const, Matrix &H) const override {
            UTOPIA_TRACE_REGION_BEGIN("IsotropicPhaseFieldForBrittleFractures::hessian");

            if (empty(H)) {
                // if(use_dense_hessian_) {
                //     H = local_zeros({this->space_.n_dofs(), this->space_.n_dofs()}); //FIXME
                // } else {
                this->space_.create_matrix(H);
                // }
            } else {
                H *= 0.0;
            }

            USpace U;
            this->space_.subspace(1, U);
            CSpace C = this->space_.subspace(0);

            // auto &x     = const_cast<Vector &>(x_const);t
            // auto &press = const_cast<Vector &>(pressure_field_);

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
            auto c_old = c_old_fun.value(q);  // E.P Added old value of damage at quadrature
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
            Strain<USpace, Quadrature> ref_strain_u(U, q);

            {
                UTOPIA_TRACE_REGION_BEGIN("IsotropicPhaseFieldForBrittleFractures::hessian_local_assembly");

                auto U_view = U.view_device();
                auto C_view = C.view_device();
                auto space_view = this->space_.view_device();

                auto c_view = c_val.view_device();
                auto c_old_view = c_old.view_device();  // E.P getting view device of old damage within scope
                auto p_view = p_val.view_device();

                auto c_grad_view = c_grad.view_device();
                auto u_view = u_val.view_device();

                auto strain_view = strain.view_device();
                auto differential_view = differential.view_device();

                // auto v_grad_shape_view = v_grad_shape.view_device();
                auto c_shape_view = c_shape.view_device();
                auto c_grad_shape_view = c_grad_shape.view_device();

                auto H_view = this->space_.assembly_view_device(H);
                auto ref_strain_u_view = ref_strain_u.view_device();

                Device::parallel_for(
                    this->space_.element_range(), UTOPIA_LAMBDA(const SizeType &i) {
                        // StaticMatrix<Scalar, Dim, Dim> strain_n, strain_p;
                        StaticMatrix<Scalar, U_NDofs + C_NDofs, U_NDofs + C_NDofs> el_mat;
                        StaticMatrix<Scalar, Dim, Dim> stress;

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

                        StaticVector<Scalar, NQuadPoints> c_old;  // E.P Added c_old vector for penalty irreversability
                        c_old_view.get(c_e, c_old);

                        auto dx = differential_view.make(c_e);
                        auto c_grad_shape_el = c_grad_shape_view.make(c_e);
                        auto c_shape_fun_el = c_shape_view.make(c_e);

                        ////////////////////////////////////////////
                        bool update_elast_tensor = false;
                        Point centroid;
                        c_e.centroid(centroid);

                        this->non_const_params().update(centroid, update_elast_tensor);

                        // Getting new material parameter values
                        const Scalar mu = this->params_.mu;
                        const Scalar lambda = this->params_.lambda;

                        ////////////////////////////////////////////

                        for (SizeType qp = 0; qp < NQuadPoints; ++qp) {
                            const Scalar tr_strain_u = trace(el_strain.strain[qp]);

                            // const Scalar eep = elastic_energy(this->params_, c[qp], tr_strain_u,
                            // el_strain.strain[qp]);
                            const Scalar eep_fix = strain_energy(this->params_, tr_strain_u, el_strain.strain[qp]);

                            // pragma GCCunroll(C_NDofs)
                            for (SizeType l = 0; l < C_NDofs; ++l) {
                                const Scalar c_shape_l = c_shape_fun_el(l, qp);
                                auto &&c_grad_l = c_grad_shape_el(l, qp);

                                // SYMMETRIC VERSION
                                for (SizeType j = l; j < C_NDofs; ++j) {
                                    const Scalar c_shape_j_l_prod = c_shape_fun_el(j, qp) * c_shape_l;

                                    Scalar val = bilinear_cc(this->params_,
                                                             c[qp],
                                                             eep_fix,
                                                             c_shape_j_l_prod,
                                                             c_grad_shape_el(j, qp),
                                                             c_grad_l) *
                                                 dx(qp);

                                    if (this->params_.use_pressure) {
                                        std::cout << "USE PRESSURE IS ON" << std::endl;
                                        val += PhaseFieldFracBase<FunctionSpace, Dim>::quadratic_degradation_deriv2(
                                                   this->params_, c[qp]) *
                                               p[qp] * tr_strain_u * c_shape_j_l_prod * dx(qp);
                                    }

                                    if (this->params_.use_penalty_irreversibility) {
                                        auto c_cold = c[qp] - c_old[qp];
                                        auto c_heaviside = c_cold <= 0.0 ? 1.0 : 0.0;
                                        val += c_heaviside * this->params_.penalty_param_irreversible *
                                               c_shape_j_l_prod * dx(qp);
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
                                    auto element_stress =
                                        2.0 * mu * u_strain_shape_el(j, qp) +
                                        lambda * trace(u_strain_shape_el(j, qp)) * (device::identity<Scalar>());
                                    Scalar val = PhaseFieldFracBase<FunctionSpace, Dim>::bilinear_uu(
                                                     this->params_,
                                                     c[qp],
                                                     // p_stress_view.stress(j, qp),       //constant material props
                                                     element_stress,  // hetero material props
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

                                compute_stress(this->params_, tr_strain_u, el_strain.strain[qp], stress);
                                for (SizeType c_i = 0; c_i < C_NDofs; ++c_i) {
                                    // CHANGE (pre-compute/store shape fun)
                                    const Scalar c_shape_i = c_shape_fun_el(c_i, qp);

                                    // #pragma clang loop unroll_count(U_MIN(U_NDofs,
                                    // UNROLL_FACTOR)) #pragma GCC unroll U_MIN(U_NDofs,
                                    // UNROLL_FACTOR)
                                    for (SizeType u_i = 0; u_i < U_NDofs; ++u_i) {
                                        auto &&strain_shape = u_strain_shape_el(
                                            u_i, qp);  // E.P ?? should this be strain shape functions or just strain??

                                        Scalar val =
                                            bilinear_uc(this->params_, c[qp], stress, strain_shape, c_shape_i) * dx(qp);

                                        if (this->params_.use_pressure) {
                                            std::cout << "USE PRESSURE IS ON" << std::endl;
                                            const Scalar tr_strain_shape = sum(diag(strain_shape));
                                            val += PhaseFieldFracBase<FunctionSpace, Dim>::quadratic_degradation_deriv(
                                                       this->params_, c[qp]) *
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
                        }  // end of quadrature point llop

                        //                         printf("\n------- Hess ------------\n");
                        //                         for (int i = 0; i < U_NDofs + C_NDofs; i++) {
                        //                             for (int j = 0; j < U_NDofs + C_NDofs; j++) {
                        //                                 double val = el_mat(i, j);
                        //                                 if (std::abs(val) < 1e-4) {
                        //                                     val = 0;
                        //                                 }
                        //                                 printf("%.4g, ", val);
                        //                             }
                        //                             printf("\n");
                        //                         }
                        //                         printf("---------------------\n");

                        space_view.add_matrix(e, el_mat, H_view);
                    });

                UTOPIA_TRACE_REGION_END("IsotropicPhaseFieldForBrittleFractures::hessian_local_assembly");
            }

            // check before boundary conditions
            if (this->check_derivatives_) {
                this->diff_ctrl_.check_hessian(*this, x_const, H);
            }

            this->space_.apply_constraints(H);

            if (this->params_.use_crack_set_irreversibiblity) {
                this->apply_zero_constraints_irreversibiblity(H, x_const);
            }

            UTOPIA_TRACE_REGION_END("IsotropicPhaseFieldForBrittleFractures::hessian");
            return true;
        }

        //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

        template <class GradShape>
        UTOPIA_INLINE_FUNCTION static Scalar bilinear_cc(const Parameters &params,
                                                         const Scalar &phase_field_value,
                                                         const Scalar &elastic_energy,
                                                         // const Scalar &shape_trial,
                                                         // const Scalar &shape_test,
                                                         const Scalar &shape_prod,
                                                         const GradShape &grad_trial,
                                                         const GradShape &grad_test) {
            return PhaseFieldFracBase<FunctionSpace, Dim>::diffusion_c(params, grad_trial, grad_test) +
                   reaction_c(params, shape_prod) +
                   elastic_deriv_cc(params, phase_field_value, elastic_energy, shape_prod);
        }

        // (sigma+(phi_u), epsilon(u)) * g'_c * phi_c
        template <class Stress, class FullStrain>
        UTOPIA_INLINE_FUNCTION static Scalar bilinear_cu(const Parameters &params,
                                                         const Scalar &phase_field_value,
                                                         const Stress &stress_p,
                                                         const FullStrain &full_strain,
                                                         const Scalar &c_trial_fun) {
            return ((1.0 - params.regularization) *
                    PhaseFieldFracBase<FunctionSpace, Dim>::quadratic_degradation_deriv(params, phase_field_value)) *
                   c_trial_fun * inner(stress_p, full_strain);
        }

        template <class Stress, class FullStrain>
        UTOPIA_INLINE_FUNCTION static Scalar bilinear_uc(const Parameters &params,
                                                         const Scalar &phase_field_value,
                                                         const Stress &stress,
                                                         const FullStrain &full_strain,
                                                         const Scalar &c_trial_fun) {
            return PhaseFieldFracBase<FunctionSpace, Dim>::quadratic_degradation_deriv(params, phase_field_value) *
                   c_trial_fun * inner(stress, full_strain);
        }

        UTOPIA_INLINE_FUNCTION static Scalar reaction_c(const Parameters &params,
                                                        // const Scalar &trial,
                                                        // const Scalar &test,
                                                        const Scalar &shape_prod) {
            return (params.fracture_toughness / params.length_scale) * shape_prod;
        }

        UTOPIA_INLINE_FUNCTION static Scalar elastic_deriv_cc(const Parameters &params,
                                                              const Scalar &phase_field_value,
                                                              const Scalar &elastic_energy,
                                                              // const Scalar &trial,
                                                              // const Scalar &test
                                                              const Scalar &shape_prod) {
            const Scalar dcc =
                (1.0 - params.regularization) *
                PhaseFieldFracBase<FunctionSpace, Dim>::quadratic_degradation_deriv2(params, phase_field_value);
            return dcc * shape_prod * elastic_energy;
        }

        template <class Strain>
        UTOPIA_INLINE_FUNCTION static Scalar grad_elastic_energy_wrt_c(const Parameters &params,
                                                                       const Scalar &phase_field_value,
                                                                       const Scalar &trace,
                                                                       const Strain &strain) {
            return (PhaseFieldFracBase<FunctionSpace, Dim>::quadratic_degradation_deriv(params, phase_field_value) *
                    (1.0 - params.regularization)) *
                   strain_energy(params, trace, strain);
        }

        template <class Strain, class Stress>
        UTOPIA_INLINE_FUNCTION static void compute_stress(const Parameters &params,
                                                          const Scalar &tr,
                                                          const Strain &strain,
                                                          Stress &stress) {
            stress = (2.0 * params.mu * strain) + (params.lambda * tr * (device::identity<Scalar>()));
        }

        template <class Grad, class Strain>
        UTOPIA_INLINE_FUNCTION static Scalar energy(const Parameters &params,
                                                    // c
                                                    const Scalar &phase_field_value,
                                                    const Grad &phase_field_grad,
                                                    // u
                                                    const Scalar &trace,
                                                    const Strain &strain) {
            //            std::cout << "frac en: " << PhaseFieldFracBase<FunctionSpace, Dim>::fracture_energy(
            //                        params, phase_field_value, phase_field_grad) << "   elas en: " <<
            //                        elastic_energy(params, phase_field_value, trace, strain) << std::endl;

            return PhaseFieldFracBase<FunctionSpace, Dim>::fracture_energy(
                       params, phase_field_value, phase_field_grad) +
                   elastic_energy(params, phase_field_value, trace, strain);
        }

        template <class Strain>
        UTOPIA_INLINE_FUNCTION static Scalar strain_energy(const Parameters &params,
                                                           const Scalar trace,
                                                           const Strain &strain) {
            return 0.5 * params.lambda * trace * trace + params.mu * inner(strain, strain);
        }

        template <class Strain>
        UTOPIA_INLINE_FUNCTION static Scalar elastic_energy(const Parameters &params,
                                                            const Scalar &phase_field_value,
                                                            const Scalar &trace,
                                                            const Strain &strain) {
            return (PhaseFieldFracBase<FunctionSpace, Dim>::quadratic_degradation(params, phase_field_value) *
                        (1.0 - params.regularization) +
                    params.regularization) *
                   strain_energy(params, trace, strain);
        }

        void write_to_file(const std::string &output_path, const Vector &x, const Scalar time) override {
            PhaseFieldFracBase<FunctionSpace, Dim>::write_to_file(output_path, x, time);

            // Post-processing functions
            // And write outputs
            this->export_strain_and_stress(output_path, x, time);
            if (mpi_world_rank() == 0) std::cout << "Saving file: " << output_path << std::endl;
        }
    };

}  // namespace utopia

// clean-up macros
#undef UNROLL_FACTOR
#undef U_MIN
#endif
