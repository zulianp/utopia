#ifndef UTOPIA_BC_SETUP_HPP
#define UTOPIA_BC_SETUP_HPP


#include "utopia_Base.hpp"
#include "utopia_RangeDevice.hpp"

//include edsl components
#include "utopia_AppRunner.hpp"
#include "utopia_Core.hpp"
#include "utopia_PetscDM.hpp"
#include "utopia_petsc_Matrix.hpp"
#include "utopia_AssemblyView.hpp"
#include "utopia_DeviceView.hpp"
#include "utopia_petsc.hpp"
#include "utopia_ConjugateGradient.hpp"
#include "utopia_TrivialPreconditioners.hpp"
#include "utopia_LaplacianView.hpp"
#include "utopia_MPITimeStatistics.hpp"
#include "utopia_BratuFE.hpp"
#include "utopia_PoissonFE.hpp"
#include "utopia_MassMatrixView.hpp"
#include "utopia_petsc_dma_FunctionSpace.hpp"
#include "utopia_petsc_DirichletBoundaryConditions.hpp"
#include "utopia_LinearElasticityView.hpp"
#include "utopia_GradInterpolate.hpp"
#include "utopia_PrincipalStrainsView.hpp"
#include "utopia_IsotropicPhaseField.hpp"
#include "utopia_FEFunction.hpp"
#include "utopia_SampleView.hpp"
#include "utopia_MPRGP.hpp"
#include "utopia_TrustRegionVariableBound.hpp"
#include "utopia_QuasiNewtonBound.hpp"
#include "utopia_Backtracking.hpp"
#include "utopia_LBFGS.hpp"
#include "utopia_QuasiTrustRegionVariableBound.hpp"
#include "utopia_InitialCondition.hpp"


#include <random>
#include <cmath>
#include <chrono>

namespace utopia {

    template<class FunctionSpace>
    class BCSetup : public Configurable
    {
        public:
            using Scalar    = typename FunctionSpace::Scalar;

            BCSetup(FunctionSpace & space):
            space_(space)
            {

            }

            virtual ~BCSetup()
            {

            }

            void read(Input &in) override
            {

            }            

            virtual void emplace_time_dependent_BC(const Scalar & time) = 0; 

        protected:
            FunctionSpace & space_; 
    };


    template<class FunctionSpace>
    class PFFracFixAllDisp : public BCSetup<FunctionSpace>
    {
        public:
            using Scalar    = typename FunctionSpace::Scalar;
            using Vector   = typename FunctionSpace::Vector;

            PFFracFixAllDisp(FunctionSpace & space): BCSetup<FunctionSpace>(space)
            {

            }

            void emplace_time_dependent_BC(const Scalar & /*time*/) override
            {
                static const int Dim   = FunctionSpace::Dim;

                using Point          = typename FunctionSpace::Point;
                this->space_.reset_bc(); 

                for(int d = 1; d < Dim + 1; ++d) {
                    this->space_.emplace_dirichlet_condition(
                        SideSet::left(),
                        UTOPIA_LAMBDA(const Point &p) -> Scalar {
                        return 0.0;
                        },
                        d
                    );

                    this->space_.emplace_dirichlet_condition(
                        SideSet::right(),
                        UTOPIA_LAMBDA(const Point &p) -> Scalar {
                            return 0.0;
                        },
                        d
                        );

                    this->space_.emplace_dirichlet_condition(
                        SideSet::top(),
                        UTOPIA_LAMBDA(const Point &p) -> Scalar {
                            return 0.0;
                        },
                        d
                        ); 

                    this->space_.emplace_dirichlet_condition(
                        SideSet::bottom(),
                        UTOPIA_LAMBDA(const Point &p) -> Scalar {
                            return 0.0;
                        },
                        d
                        );


                    if(Dim ==3){
                        this->space_.emplace_dirichlet_condition(
                            SideSet::front(),
                            UTOPIA_LAMBDA(const Point &p) -> Scalar {
                                return 0.0;
                            },
                            d
                            );

                        this->space_.emplace_dirichlet_condition(
                            SideSet::back(),
                            UTOPIA_LAMBDA(const Point &p) -> Scalar {
                                return 0.0;
                            },
                            d
                            );                                                    
                    }

                }
            }
    };   




    template<class FunctionSpace>
    class PFFracTension2D : public BCSetup<FunctionSpace>
    {
        public:
            using Scalar    = typename FunctionSpace::Scalar;
            using Vector    = typename FunctionSpace::Vector;

            PFFracTension2D(FunctionSpace & space, const Scalar & disp_y=1.0): BCSetup<FunctionSpace>(space), disp_y_(disp_y)
            {

            }

            void read(Input &in) override
            {
                in.get("disp_y", disp_y_);   
            }               

            void emplace_time_dependent_BC(const Scalar & time) override
            {
                static const int Dim   = FunctionSpace::Dim;

                using Point          = typename FunctionSpace::Point;
                this->space_.reset_bc(); 

                this->space_.emplace_dirichlet_condition(
                    SideSet::left(),
                    UTOPIA_LAMBDA(const Point &p) -> Scalar {
                        return 0.0;
                    },
                    1
                    );

                this->space_.emplace_dirichlet_condition(
                    SideSet::right(),
                    UTOPIA_LAMBDA(const Point &p) -> Scalar {
                        return disp_y_ * time;
                    },
                    1
                    );

                for(int d = 2; d < Dim + 1; ++d) {
                    this->space_.emplace_dirichlet_condition(
                        SideSet::left(),
                        UTOPIA_LAMBDA(const Point &p) -> Scalar {
                        return 0.0;
                        },
                        d
                    );

                    this->space_.emplace_dirichlet_condition(
                        SideSet::right(),
                        UTOPIA_LAMBDA(const Point &p) -> Scalar {
                            return 0.0;
                        },
                        d
                        );
                    }
            }

            private:
                Scalar disp_y_; 
    };        
}


#endif // UTOPIA_BC_SETUP_HPP