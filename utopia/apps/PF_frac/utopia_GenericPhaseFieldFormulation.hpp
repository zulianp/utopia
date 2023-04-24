#ifndef UTOPIA_GENERIC_PHASE_FIELD_FORMULATION_HPP
#define UTOPIA_GENERIC_PHASE_FIELD_FORMULATION_HPP

#include "utopia_DeviceTensorContraction.hpp"
#include "utopia_DeviceTensorProduct.hpp"
#include "utopia_DiffController.hpp"
#include "utopia_ExtendedFunction.hpp"
#include "utopia_FEFunction.hpp"
#include "utopia_GradInterpolate.hpp"
#include "utopia_LinearElasticityView.hpp"
#include "utopia_PhaseFieldBase.hpp"
#include "utopia_PrincipalShapeStressView.hpp"
#include "utopia_PrincipalStrainsView.hpp"
#include "utopia_StrainView.hpp"
#include "utopia_TensorView4.hpp"
#include "utopia_Tracer.hpp"
#include "utopia_petsc_NeumannBoundaryConditions.hpp"

#define UNROLL_FACTOR 4
#define U_MIN(a, b) ((a) < (b) ? (a) : (b))

namespace utopia {

    // E.P See linear softening law in Wu & Nguyen 2017 - Length-scale insensitive phase field damage model for brittle fracture
    struct CHZ_Linear {

          static constexpr double a2_ = -0.5;

           template<class FunctionSpace>
           UTOPIA_INLINE_FUNCTION static double damage_normalisation(const PFFracParameters<FunctionSpace> & p) {
               // c0 = pi
               return  p.fracture_toughness/(p.length_scale*M_PI);
           }

           template <typename C>
           UTOPIA_INLINE_FUNCTION static C local_dissipation(const C &c) {
               // alpha * [( epsi + (1-epsi)* alpha]   where epsi = 2 for linear

               return 2.0*c - c*c;
           }

           template <typename C>
           UTOPIA_INLINE_FUNCTION static C local_dissipation_deriv(const C &c) {
               return 2.0 - 2.0*c;
           }

           template <typename C>
           UTOPIA_INLINE_FUNCTION static C local_dissipation_deriv2(const C &) {
               return -2.0;
           }

           // E.P: Penalty for AT1 Model is used herein
           // this computation follows eq. 60 from "On penalization in variational
           // phase-field models of britlle fracture, Gerasimov, Lorenzis"
           template<class FunctionSpace>
           static void configure_penalty_irreversibility( PFFracParameters<FunctionSpace> & p){
               assert(p.use_penalty_irreversibility);
               typename FunctionSpace::Scalar tol2 = p.penalty_tol * p.penalty_tol;
               p.penalty_param_irreversible = p.fracture_toughness / p.length_scale * (27.0 /(64.0 * tol2) );
           }

           template<class FunctionSpace>
           static void configure_penalty_non_negative( PFFracParameters<FunctionSpace> & p){
               assert(p.use_penalty_irreversibility);
               typename FunctionSpace::Scalar L = (p.Length_x + p.Length_y + p.Length_z)/FunctionSpace::Dim;
               p.penalty_param_non_neg = p.fracture_toughness / p.length_scale * 9.0/64.0 * (L/p.length_scale - 2.0 )/(p.penalty_tol_non_neg);
               if (mpi_world_rank()==0)
                 utopia::out() << "Lengthscale: " << p.length_scale << "  Penalty CHZ_linear n_neg: " << p.penalty_param_non_neg << std::endl;
           }


           template<typename C, class FunctionSpace>
           UTOPIA_INLINE_FUNCTION static C degradation(const C &c, const PFFracParameters<FunctionSpace> & p) {
               const double a1_ = 4.0/(M_PI*p.length_scale) * p.E * p.fracture_toughness / (p.tensile_strength*p.tensile_strength) ;
               C imc  = 1.0 - c;
               C imc2 = imc*imc;
               return imc2/(imc2 + a1_*c*(1. + a2_*c));
           }

           template <typename C, class FunctionSpace>
           UTOPIA_INLINE_FUNCTION static C degradation_deriv( const C &c, const PFFracParameters<FunctionSpace> & p) {
               const double a1_ = 4.0/(M_PI*p.length_scale) * p.E * p.fracture_toughness / (p.tensile_strength*p.tensile_strength) ;
               C imc  = 1.0 - c;
               C imc2 = imc*imc;
               return - a1_*imc*(2.0*a2_*c + c + 1.0)/std::pow(a1_*c*(a2_*c+1.0) + imc2  ,2.0);
           }

           template <typename C, class FunctionSpace>
           UTOPIA_INLINE_FUNCTION static C degradation_deriv2( const C &c, const PFFracParameters<FunctionSpace> & p) {
               const double a1_ = 4.0/(M_PI*p.length_scale) * p.E * p.fracture_toughness / (p.tensile_strength*p.tensile_strength) ;
               C imc  = 1.0 - c;
               C imc2 = imc*imc;
               return 2.0*a1_*(2.0*a2_*c + c - a2_)/std::pow(a1_*c*(a2_*c+1.0) + imc2  ,2.0)
                       + 2.0*a1_*imc*(2.0*a2_*c + c + 1.0)*(2.0*a1_*a2_*c + a1_ - 2.0*imc )/std::pow(a1_*c*(a2_*c+1.0) + imc2  , 3.0);
           }

           static const bool penalise_negative_phase_field_values = false; //NOT WORKING, but AT1 models need to penalise negative phase field values

           static const bool enforce_min_crack_driving_force = false; //NOT WORKING, keep false.

           template<class FunctionSpace>
           UTOPIA_INLINE_FUNCTION static double min_crack_driving_force(const PFFracParameters<FunctionSpace> & p) {
                      return 0.5*p.tensile_strength*p.tensile_strength/p.E;
                  }


           //Value of phase field when damage localises in 1D bar
           static double CriticalDamage() {return 0.0 ; } //elastic phase


       };


    struct AT1 {
    public:
        template <class FunctionSpace>
        UTOPIA_INLINE_FUNCTION static double damage_normalisation(const PFFracParameters<FunctionSpace> &p) {
            return 3.0 / 8.0 * p.fracture_toughness / p.length_scale;
        }

        template <typename C>
        UTOPIA_INLINE_FUNCTION static C local_dissipation(const C &c) {
            return c;
        }

        template <typename C>
        UTOPIA_INLINE_FUNCTION static C local_dissipation_deriv(const C & c) {
            return 1.0;
        }

        template <typename C>
        UTOPIA_INLINE_FUNCTION static C local_dissipation_deriv2(const C &) {
            return 0.;
        }

        // E.P: Penalty for AT1 Model
        // this computation follows eq. 60 from "On penalization in variational
        // phase-field models of britlle fracture, Gerasimov, Lorenzis"
        template <class FunctionSpace>
        static void configure_penalty_irreversibility(PFFracParameters<FunctionSpace> &p) {
            assert(p.use_penalty_irreversibility);
            typename FunctionSpace::Scalar tol2 = p.penalty_tol * p.penalty_tol;
            p.penalty_param_irreversible = p.fracture_toughness / p.length_scale * (27.0 / (64.0 * tol2));
            if (mpi_world_rank() == 0)
                utopia::out() << "Lengthscale: " << p.length_scale
                              << "  Penalty AT1 Irrev: " << p.penalty_param_irreversible << std::endl;
        }

        template <class FunctionSpace>
        static void configure_penalty_non_negative(PFFracParameters<FunctionSpace> &p) {
            assert(p.use_penalty_irreversibility);
            typename FunctionSpace::Scalar L = (p.Length_x + p.Length_y + p.Length_z) / FunctionSpace::Dim;
            p.penalty_param_non_neg = p.fracture_toughness / p.length_scale * 9.0 / 64.0 * (L / p.length_scale - 2.0) /
                                      (p.penalty_tol_non_neg);
            if (mpi_world_rank() == 0)
                utopia::out() << "Lengthscale: " << p.length_scale << "  Penalty AT1 n_neg: " << p.penalty_param_non_neg
                              << std::endl;
        }

       template <typename C, class FunctionSpace>
       UTOPIA_INLINE_FUNCTION static C degradation(const C &c, const PFFracParameters<FunctionSpace> & ) {
           C imc = 1.0 - c;
           return imc * imc;
       }

       template <typename C, class FunctionSpace>
       UTOPIA_INLINE_FUNCTION static C degradation_deriv( const C &c, const PFFracParameters<FunctionSpace> & ) {
           C imc = 1.0 - c;
           return -2.0 * imc;
       }

       template <typename C, class FunctionSpace>
       UTOPIA_INLINE_FUNCTION static C degradation_deriv2( const C &, const PFFracParameters<FunctionSpace> & ) {
           return 2.0;
       }

        static const bool penalise_negative_phase_field_values =
            false;  // NOT WORKING, but AT1 models need to penalise negative phase field values
        
        static const bool enforce_min_crack_driving_force = false;

        template <class FunctionSpace>
        UTOPIA_INLINE_FUNCTION static double min_crack_driving_force(const PFFracParameters<FunctionSpace> &p) {
            return 3.0 / 16.0 * p.fracture_toughness / p.length_scale;
        }

       //Value of phase field when damage localises in 1D bar
       static double CriticalDamage() {return 0.0 ; } //elastic phase

    };

    struct AT2 {
    public:
        template <class FunctionSpace>
        UTOPIA_INLINE_FUNCTION static double damage_normalisation(const PFFracParameters<FunctionSpace> &p) {
            return 1.0 / 2.0 * p.fracture_toughness / p.length_scale;
        }

        template <typename C>
        UTOPIA_INLINE_FUNCTION static C local_dissipation(const C &c) {
            return c * c;
        }

        template <typename C>
        UTOPIA_INLINE_FUNCTION static C local_dissipation_deriv(const C &c) {
            return 2.0 * c;
        }

        template <typename C>
        UTOPIA_INLINE_FUNCTION static C local_dissipation_deriv2(const C &c) {
            return 2.0;
        }

        // E.P: AT2 Model
        // this computation follows eq. 50 from "On penalization in variational
        // phase-field models of britlle fracture, Gerasimov, Lorenzis"
        template <class FunctionSpace>
        static void configure_penalty_irreversibility(PFFracParameters<FunctionSpace> &p) {
            assert(p.use_penalty_irreversibility);
            typename FunctionSpace::Scalar tol2 = p.penalty_tol * p.penalty_tol;
            p.penalty_param_irreversible = p.fracture_toughness / p.length_scale * (1.0 / tol2 - 1.0);
            if (mpi_world_rank() == 0)
                utopia::out() << "Lengthscale: " << p.length_scale
                              << "  Penalty AT2 Irrev: " << p.penalty_param_irreversible << std::endl;
        }

        template <class FunctionSpace>
        static void configure_penalty_non_negative(PFFracParameters<FunctionSpace> &p) {
            p.penalty_param_non_neg = 0.0;  // not needed in AT2
            if (mpi_world_rank() == 0)
                utopia::out() << "Lengthscale: " << p.length_scale << "  Penalty AT2 n_neg: " << p.penalty_param_non_neg
                              << std::endl;
        }

        

       template <typename C, class FunctionSpace>
       UTOPIA_INLINE_FUNCTION static C degradation(const C &c, const PFFracParameters<FunctionSpace> &) {
           C imc = 1.0 - c;
           return imc * imc;
       }

       template <typename C, class FunctionSpace>
       UTOPIA_INLINE_FUNCTION static C degradation_deriv( const C &c, const PFFracParameters<FunctionSpace> &) {
           C imc = 1.0 - c;
           return -2.0 * imc;
       }

       template <typename C, class FunctionSpace>
       UTOPIA_INLINE_FUNCTION static C degradation_deriv2( const C &, const PFFracParameters<FunctionSpace> &) {
           return 2.0;
       }

        static const bool penalise_negative_phase_field_values =
            false;  // AT2 models do not need to penalise negative phase field values

        static const bool enforce_min_crack_driving_force = false;

        template <class FunctionSpace>
        UTOPIA_INLINE_FUNCTION static double min_crack_driving_force(const PFFracParameters<FunctionSpace> &p) {
            return 3.0 / 16.0 * p.fracture_toughness / p.length_scale;
        }
      
        //Value of phase field when damage localises in 1D bar (disregarding the stability of the bar)
        static double CriticalDamage() {return 0.25 ; } //elastic phase
  
    };
   


    template <class FunctionSpace, int Dim = FunctionSpace::Dim, class PFFormulation = AT1>
    class GenericPhaseFieldFormulation : public PhaseFieldFracBase<FunctionSpace, Dim> {
    public:
        using Scalar = typename FunctionSpace::Scalar;
        using SizeType = typename FunctionSpace::SizeType;
        using Vector = typename FunctionSpace::Vector;
        using Matrix = typename FunctionSpace::Matrix;
        using Device = typename FunctionSpace::Device;

        using USpace = typename FunctionSpace::template Subspace<Dim>;
        using CSpace = typename FunctionSpace::template Subspace<1>;

        using UElem = typename USpace::ViewDevice::Elem;
        using CElem = typename CSpace::ViewDevice::Elem;
        using MixedElem = typename FunctionSpace::ViewDevice::Elem;

        using Parameters = typename PhaseFieldFracBase<FunctionSpace, Dim>::PFFracParameters;
        using HeteroParamsFunction = typename Parameters::HeteroParamsFunction;

        using Point = typename FunctionSpace::Point;

        // FIXME
        using Shape = typename FunctionSpace::Shape;
        using Quadrature = utopia::Quadrature<Shape, 2 * (Shape::Order)>;

        static const int C_NDofs = CSpace::NDofs;
        static const int U_NDofs = USpace::NDofs;

        static const int NQuadPoints = Quadrature::NPoints;

        GenericPhaseFieldFormulation(FunctionSpace &space) : PhaseFieldFracBase<FunctionSpace, Dim>(space) {}

        GenericPhaseFieldFormulation(FunctionSpace &space, const Parameters &params)
            : PhaseFieldFracBase<FunctionSpace, Dim>(space, params) {}

        //Tensile strength (for models with elastic phase AT1 and CHZ )
        double CriticalStress(){
            // sigma_c = Eo * Uc / L
            return this->params_.E * CriticalDisplacement() / this->params_.Length_x;
        }

        //For a 1D homogeneous Bar
        double CriticalDisplacement(){
            const double c = PFFormulation::CriticalDamage();
            // Uc = L * sqrt( - 2 w'/E' )
            return this->params_.Length_x * 
                std::sqrt(- 2.0 * PFFormulation::damage_normalisation(this->params_) * PFFormulation::local_dissipation_deriv(c)
                            / (PFFormulation::degradation_deriv(c, this->params_) * this->params_.E)  ); }

        void read(Input &in) {
            PhaseFieldFracBase<FunctionSpace, Dim>::read(in);

            if (this->params_.use_penalty_irreversibility) {
                PFFormulation::configure_penalty_irreversibility(this->params_);
                PFFormulation::configure_penalty_non_negative(this->params_);
            }

            if (mpi_world_rank() == 0) {
                double disp_x;
                in.get("disp_x", disp_x);
                std::cout << "Tensile Strength:  " << CriticalStress()
                          << "\nCrit Displacement: " << CriticalDisplacement()
                          << "\nTime at failure: " << CriticalDisplacement() / disp_x << std::endl;
            }
        }
        ////////////////////////////////////////////////////////////////////////////////////
        virtual bool fracture_energy(const Vector & /*x_const*/, Scalar & /*val*/) const = 0;
        virtual bool elastic_energy(const Vector & /*x_const*/, Scalar & /*val*/) const = 0;

        ////////////////////////////////////////////////////////////////////////////////////

        /// VALUE TERMS /////
        // E.P Generic implementation of   w1(w + l^2 |nabla alpha|^2 )
        template <typename PhaseFieldValue, class Grad>
        UTOPIA_INLINE_FUNCTION static PhaseFieldValue fracture_energy(const Parameters &params,
                                                                      const PhaseFieldValue &phase_field_value,
                                                                      const Grad &phase_field_grad) {
            return PFFormulation::damage_normalisation(params) *
                   (PFFormulation::local_dissipation(phase_field_value) +
                    params.length_scale * params.length_scale * inner(phase_field_grad, phase_field_grad));
        }

        /// Gradient Terms ////
        template <typename PhaseFieldValue, class Grad, typename TestFunction, class GradTest>
        UTOPIA_INLINE_FUNCTION static PhaseFieldValue grad_fracture_energy_wrt_c(
            const Parameters &params,
            const PhaseFieldValue &phase_field_value,
            const Grad &phase_field_grad,
            const TestFunction &test_function,
            const GradTest &grad_test_function) {
            return PFFormulation::damage_normalisation(params) *
                   (PFFormulation::local_dissipation_deriv(phase_field_value) * test_function +
                    2.0 * params.length_scale * params.length_scale * inner(phase_field_grad, grad_test_function));
        }

        /// Hessian Terms /////
        template <class Grad>
        UTOPIA_INLINE_FUNCTION static auto diffusion_c(const Parameters &params,
                                                       const Grad &g_trial,
                                                       const Grad &g_test) {
            return PFFormulation::damage_normalisation(params) * 2 * params.length_scale * params.length_scale *
                   inner(g_trial, g_test);
        }

        template <typename PhaseFieldValue>
        UTOPIA_INLINE_FUNCTION static Scalar reaction_c(const Parameters &params,
                                                        const PhaseFieldValue &phase_field_value,
                                                        // const Scalar &trial,
                                                        // const Scalar &test,
                                                        const Scalar &shape_prod) {
            // std::cout << PFFormulation::damage_normalisation(params) *
            // PFFormulation::local_dissipation_deriv2(phase_field_value)  * shape_prod << std::endl;
            return PFFormulation::damage_normalisation(params) *
                   PFFormulation::local_dissipation_deriv2(phase_field_value) * shape_prod;
        }

        template <typename PhaseFieldValue, class StressShape, class Grad>
        UTOPIA_INLINE_FUNCTION static PhaseFieldValue bilinear_uu(const Parameters &params,
                                                                  const PhaseFieldValue &phase_field_value,
                                                                  const StressShape &stress,
                                                                  const Grad &strain_test) {
            const auto gc =
                ((1.0 - params.regularization) * PFFormulation::degradation(phase_field_value, params) + params.regularization);
            return inner(gc * stress, strain_test);
        }

        bool fracture_energy_in_middle_layer(const Vector &x_const, Scalar &val) const {
            UTOPIA_TRACE_REGION_BEGIN("GenericPhaseFieldFormulation::fracture_energy_in_middle_layer");

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

                        if (centroid[1] > this->non_const_params().bottom_layer_height &&
                            centroid[1] < this->non_const_params().top_layer_height ){
                            //only for middle layer
                            for (SizeType qp = 0; qp < NQuadPoints; ++qp) {
                                el_energy +=
                                    GenericPhaseFieldFormulation<FunctionSpace, Dim,PFFormulation>::fracture_energy(
                                        this->params_, c[qp], c_grad_el[qp]) *
                                        dx(qp);
                            }
                        }
                        assert(el_energy == el_energy);
                        return el_energy;
                    },
                    val);
            }

            val = x_const.comm().sum(val);

            UTOPIA_TRACE_REGION_END("GenericPhaseFieldFormulation::fracture_energy_in_middle_layer");
            return true;
        }


        virtual bool export_strain_and_stress(std::string output_path, const Vector &x_const, const Scalar time) const override{
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
                                auto &epsi = el_strain.strain[qp];

                                // calculate stress at quadrature
                                const Scalar tr_strain_u = trace(el_strain.strain[qp]);
                                this->compute_stress(this->params_,
                                               tr_strain_u,
                                               el_strain.strain[qp],
                                               stress);  // gets stress at quadrature point

                                strain_value +=
                                    epsi * shape * weight;  // matrix of strains added to existing nodal strain (

                                stress_value += PFFormulation::degradation(c[qp], this->params_) * stress * shape * weight;  // Sum stress at integration point

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

        bool export_material_params(std::string output_path){
            UTOPIA_TRACE_REGION_BEGIN("GenericPhaseFieldFormulation::export_mechanical_params");

            static const int total_components = 7.0; //E, nu, Gc, lambda, mu, tensile_strength, crit_disp

            using PSpace = typename FunctionSpace::template Subspace<total_components>;
            using SElem  = typename PSpace::ViewDevice::Elem;
            using WSpace = typename FunctionSpace::template Subspace<1>;

            Vector w;
            Vector g;

            /// Creating strain subspace
            // cloning mesh
            auto param_mesh = this->space_.mesh().clone(total_components);
            assert(param_mesh->n_components() == total_components);
            // Creating Subspace with cloned mesh

            PSpace S(std::move(param_mesh));
            WSpace C(this->space_.mesh().clone(1));


            S.create_vector(g);
            C.create_vector(w);
            ///////////////////////////////////////////////////////////////////////////

            {
                ////////////////////////////////////////////////////////////////////////////



                auto S_view = S.view_device();
                auto C_view = C.view_device();

                // Preparing the vector for which the parameter function space knows the dimensions (nodes*components), so
                // that we can write on this later
                auto g_view = S.assembly_view_device(g);
                auto w_view = C.assembly_view_device(w);



            Device::parallel_for(
                this->space_.element_range(), UTOPIA_LAMBDA(const SizeType &i) {
                    StaticVector<Scalar, total_components * C_NDofs> material_params;
                    StaticVector<Scalar, C_NDofs> node_count;
                    material_params.set(.0);
                    node_count.set(1.0);

                    SElem s_e;
                    S_view.elem(i, s_e);  // just needed for add_vector into g (and node coords)

                    CElem c_e;
                    C_view.elem(i, c_e);  // getting element for storing wieghts in CSpace

                    for (SizeType n = 0; n < C_NDofs; n++) {

                        ////////////////////////////////////////////
                        bool update_elast_tensor = false;
                        Point coord;
                        s_e.node(n, coord);
                        this->non_const_params().update(coord, update_elast_tensor);
                        ////////////////////////////////////////////

                        Scalar Gc = this->params_.fracture_toughness;
                        Scalar mu = this->params_.mu ;
                        Scalar l  = this->params_.lambda;
                        Scalar E  = mu*(3.*l + 2.*mu)/(l+mu);
                        Scalar nu = E/(2.0*mu) - 1.;
                        Scalar tens_strength = CriticalStress();  //std::pow(3./8. * Gc*E/this->params_.length_scale, 0.5);  //AT1 hardcoded
                        Scalar crit_disp = CriticalDisplacement();


                        material_params[n] = E;
                        material_params[C_NDofs   + n] = nu;
                        material_params[C_NDofs*2 + n] = Gc;
                        material_params[C_NDofs*3 + n] = l;
                        material_params[C_NDofs*4 + n] = mu;
                        material_params[C_NDofs*5 + n] = tens_strength;
                        material_params[C_NDofs*6 + n] = crit_disp;

                    }

                    C_view.add_vector(c_e, node_count, w_view);
                    S_view.add_vector(s_e, material_params, g_view);
                        }); //end of parallel loop
            }

            {

                auto mat_params = local_view_device(g);
                auto node_count = local_view_device(w);
                auto r = local_range_device(w);  // range of vector w (using primitivo di utopio)
                parallel_for(
                    r, UTOPIA_LAMBDA(int i) {
                        auto wi = node_count.get(i);  // extracts vector component
                        for (int k = 0; k < total_components; k++) {

                            int nodal_offset = i * total_components;                //vector g is 2*straincomponents bigger than vector of weights w
                            auto si = mat_params.get(nodal_offset + k);      //get k'th strain corresponding to node i with weight i
                            mat_params.set(nodal_offset + k, si / wi);       //normalise the strain value by the weight wi

                        }
                    });
            }  // incase backed PETSC needs synchronisation (create view in scopes and destroy them when not needed)


            rename("E nu Gc lambda mu ft Uc", g);
            output_path += "_params.vtr";
            S.write(output_path, g);  // Function space knows how to write

            UTOPIA_TRACE_REGION_END("GenericPhaseFieldFormulation::export_mechanical_params");
            return true;

        }

    protected:
    };

}  // namespace utopia

// clean-up macros
#undef UNROLL_FACTOR
#undef U_MIN
#endif
