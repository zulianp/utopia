#ifndef UTOPIA_GENERIC_PHASE_FIELD_FORMULATION_HPP
#define UTOPIA_GENERIC_PHASE_FIELD_FORMULATION_HPP

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
#include "utopia_Tracer.hpp"
#include "utopia_Views.hpp"
#include "utopia_petsc_NeumannBoundaryConditions.hpp"
#include "utopia_PhaseFieldBase.hpp"

#define UNROLL_FACTOR 4
#define U_MIN(a, b) ((a) < (b) ? (a) : (b))

namespace utopia {


struct AT1 {

public:
       template<class FunctionSpace>
       UTOPIA_INLINE_FUNCTION static double damage_normalisation(const PFFracParameters<FunctionSpace> & p) {
           return 3.0/8.0 * p.fracture_toughness/p.length_scale ;
       }

       template <typename C>
       UTOPIA_INLINE_FUNCTION static C local_dissipation(const C &c) {
           return c;
       }

       template <typename C>
       UTOPIA_INLINE_FUNCTION static C local_dissipation_deriv(const C &c) {
           return 1.0;
       }

       template <typename C>
       UTOPIA_INLINE_FUNCTION static C local_dissipation_deriv2(const C &c) {
           return 0.;
       }


       template <typename C>
       UTOPIA_INLINE_FUNCTION static C degradation(const C &c) {
           C imc = 1.0 - c;
           return imc * imc;
       }

       template <typename C>
       UTOPIA_INLINE_FUNCTION static C degradation_deriv( const C &c) {
           C imc = 1.0 - c;
           return -2.0 * imc;
       }

       template <typename C>
       UTOPIA_INLINE_FUNCTION static C degradation_deriv2( const C &) {
           return -2.0; //THIS WAS CHANGED TO NEGATIVE
       }

   };


struct AT2 {

public:
       template<class FunctionSpace>
       UTOPIA_INLINE_FUNCTION static double damage_normalisation(const PFFracParameters<FunctionSpace> & p) {
           return 1.0/2.0 * p.fracture_toughness/p.length_scale ;
       }

       template <typename C>
       UTOPIA_INLINE_FUNCTION static C local_dissipation(const C &c) {
           return c*c;
       }

       template <typename C>
       UTOPIA_INLINE_FUNCTION static C local_dissipation_deriv(const C &c) {
           return 2.0*c;
       }

       template <typename C>
       UTOPIA_INLINE_FUNCTION static C local_dissipation_deriv2(const C &c) {
           return 2.0;
       }


       template <typename C>
       UTOPIA_INLINE_FUNCTION static C degradation(const C &c) {
           C imc = 1.0 - c;
           return imc * imc;
       }

       template <typename C>
       UTOPIA_INLINE_FUNCTION static C degradation_deriv( const C &c) {
           C imc = 1.0 - c;
           return -2.0 * imc;
       }

       template <typename C>
       UTOPIA_INLINE_FUNCTION static C degradation_deriv2( const C &) {
           return -2.0; //THIS WAS CHANGED TO NEGATIVE
       }

   };








    template <class FunctionSpace, int Dim = FunctionSpace::Dim, class PFFormulation=AT1>
    class GenericPhaseFieldFormulation : public PhaseFieldFracBase<FunctionSpace,Dim> {
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

        using Parameters = typename PhaseFieldFracBase<FunctionSpace,Dim>::PFFracParameters;
        using HeteroParamsFunction = typename Parameters::HeteroParamsFunction;

        // FIXME
        using Shape = typename FunctionSpace::Shape;
        using Quadrature = utopia::Quadrature<Shape, 2 * (Shape::Order)>;

        static const int C_NDofs = CSpace::NDofs;
        static const int U_NDofs = USpace::NDofs;

        static const int NQuadPoints = Quadrature::NPoints;

        GenericPhaseFieldFormulation(FunctionSpace &space) :
            PhaseFieldFracBase<FunctionSpace, Dim>(space)
        { //No extra construction
        }


        GenericPhaseFieldFormulation(FunctionSpace &space, const Parameters &params) : PhaseFieldFracBase<FunctionSpace,Dim>(space, params)
        {

        }

        ////////////////////////////////////////////////////////////////////////////////////
        virtual bool fracture_energy(const Vector & /*x_const*/, Scalar & /*val*/) const = 0;
        virtual bool elastic_energy(const Vector & /*x_const*/, Scalar & /*val*/) const = 0;

        ////////////////////////////////////////////////////////////////////////////////////



        /// VALUE TERMS /////
        //E.P Generic implementation of   w1(w + l^2 |nabla alpha|^2 )
        template <typename PhaseFieldValue, class Grad>
        UTOPIA_INLINE_FUNCTION static PhaseFieldValue fracture_energy(const Parameters &params,
                                                                      const PhaseFieldValue &phase_field_value,
                                                                      const Grad &phase_field_grad) {
           return  pf_formulation_.damage_normalisation(params) *
                    ( pf_formulation_.local_dissipation(phase_field_value) + params.length_scale*params.length_scale * inner(phase_field_grad, phase_field_grad) );
        }


        /// Gradient Terms ////
        template <typename PhaseFieldValue, class Grad, typename TestFunction, class GradTest>
        UTOPIA_INLINE_FUNCTION static PhaseFieldValue grad_fracture_energy_wrt_c(
            const Parameters &params,
            const PhaseFieldValue &phase_field_value,
            const Grad &phase_field_grad,
            const TestFunction &test_function,
            const GradTest &grad_test_function) {
            return pf_formulation_.damage_normalisation(params) *
                    ( pf_formulation_.local_dissipation_deriv(phase_field_value)*test_function  +
                      2.0*params.length_scale*params.length_scale * inner(phase_field_grad, grad_test_function) );
            }


        /// Hessian Terms /////
        template <class Grad>
        UTOPIA_INLINE_FUNCTION static auto diffusion_c(const Parameters &params,
                                                       const Grad &g_trial,
                                                       const Grad &g_test) {
            return pf_formulation_.damage_normalisation(params) * 2 * params.length_scale * params.length_scale * inner(g_trial, g_test);
        }

        template <typename PhaseFieldValue>
        UTOPIA_INLINE_FUNCTION static Scalar reaction_c(const Parameters &params,
                                                        const PhaseFieldValue & phase_field_value,
                                                        // const Scalar &trial,
                                                        // const Scalar &test,
                                                        const Scalar &shape_prod) {
            return pf_formulation_.damage_normalisation(params) * pf_formulation_.local_dissipation_deriv2(phase_field_value)  * shape_prod;
        }



        template <typename PhaseFieldValue, class StressShape, class Grad>
        UTOPIA_INLINE_FUNCTION static PhaseFieldValue bilinear_uu(const Parameters &params,
                                                                  const PhaseFieldValue &phase_field_value,
                                                                  const StressShape &stress,
                                                                  const Grad &strain_test) {
            const auto gc = ((1.0 - params.regularization) * pf_formulation_.degradation(phase_field_value) +
                             params.regularization);
            return inner(gc * stress, strain_test);
        }



    protected:

        static PFFormulation pf_formulation_;      //Phase field formulation

    };

}  // namespace utopia

// clean-up macros
#undef UNROLL_FACTOR
#undef U_MIN
#endif
