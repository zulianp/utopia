#ifndef UTOPIA_PF2D_FRAC_NET_SAMPLER_HPP
#define UTOPIA_PF2D_FRAC_NET_SAMPLER_HPP


#include "utopia_Base.hpp"

//FIXME: this file causes nvcc to fail
#ifndef KOKKOS_ENABLE_CUDA

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
#include "utopia_MPITimeStatistics.hpp"
#include "utopia_MassMatrixView.hpp"
#include "utopia_petsc_dma_FunctionSpace.hpp"
#include "utopia_petsc_DirichletBoundaryConditions.hpp"
#include "utopia_LinearElasticityView.hpp"
#include "utopia_GradInterpolate.hpp"
#include "utopia_PrincipalStrainsView.hpp"
#include "utopia_IsotropicPhaseField.hpp"
#include "utopia_FEFunction.hpp"
#include "utopia_SampleView.hpp"
#include "utopia_InitialConditionBase.hpp"

#include <random>
#include <cmath>
#include <chrono>

namespace utopia {

    template<class T>
    struct Point2D{

        Point2D()
        {

        }

        Point2D(const T & xx, const T & yy):
        x(xx), y(yy)
        {

        }

        T x;
        T y;

        void describe()
        {
            std::cout<<"(" << x << " , "<< y << " ) \n";
        }
    };


    template<class T>
    class FracNetSamplerParams2D : public Configurable
    {

        public:
            FracNetSamplerParams2D(): 
            x_min(0.0), 
            x_max(1.0), 
            y_min(0.0), 
            y_max(1.0),
            alpha_min(0.0), 
            alpha_max(180), 
            pow_dist_coef_length(2.8), 
            min_length(0.0), 
            seed(3)
            {

            }                  

            void read(Input &in) override
            {
                in.get("x_min", x_min);
                in.get("x_max", x_max);

                in.get("y_min", y_min);
                in.get("y_max", y_max);

                in.get("alpha_min", alpha_min);
                in.get("alpha_max", alpha_max);

                in.get("pow_dist_coef_length", pow_dist_coef_length);
                in.get("min_length", min_length);

                in.get("seed", seed); 
            }

        public:
            T x_min; 
            T x_max; 

            T y_min; 
            T y_max;                 

            T alpha_min; 
            T alpha_max; 

            T pow_dist_coef_length;     
            T min_length;    

            T seed;                          
    };    

    template<class T>
    class Rectangle
    {
        public:
            Rectangle(const Point2D<T> & A, const Point2D<T> & B, const Point2D<T> & C, const Point2D<T> & D):
            A_(A), B_(B), C_(C), D_(D)
            {

            }

            Rectangle(const T & width, const FracNetSamplerParams2D<T> & params)
            {
                randomly_generate(width, params);
            }

            Rectangle(const Point2D<T> & A, const T & length, const T & width, const T & theta ):
            A_(A)
            {
                this->generate_rectangle(length, width, theta);
            }

            // generate from line 
            Rectangle(const Point2D<T> & p1, const Point2D<T> & p2, const T & width)
            {
                Point2D<T>  V, P, N; 
                V.x = p2.x - p1.x;
                V.y = p2.y - p1.y;

                P.x = V.y; 
                P.y = -V.x; 

                auto length = std::sqrt(P.x * P.x + P.y * P.y); 
                N.x = P.x /length;
                N.y = P.y /length;                 

                A_.x = p1.x + N.x * width / 2.0;
                A_.y = p1.y + N.y * width / 2.0;
                B_.x = p1.x - N.x * width / 2.0;
                B_.y = p1.y - N.y * width / 2.0;
                C_.x = p2.x + N.x * width / 2.0;
                C_.y = p2.y + N.y * width / 2.0;
                D_.x = p2.x - N.x * width / 2.0;
                D_.y = p2.y - N.y * width / 2.0;                 
            }


            bool belongs_to_rectangle(const T & x_coord, const T & y_coord)
            {
                Point2D<T> M;
                M.x = x_coord;
                M.y = y_coord;

                return belongs_to_rectangle(M);
            }

            bool belongs_to_rectangle(Point2D<T> M)
            {
                Point2D<T> AB, AM, BD, BM;
                build_vector(A_, B_, AB);
                build_vector(A_, M, AM);
                build_vector(B_, D_, BD);
                build_vector(B_, M, BM);

                T dotABAM = vec_dot(AB, AM);
                T dotABAB = vec_dot(AB, AB);
                T dotBDBM = vec_dot(BD, BM);
                T dotBDBD = vec_dot(BD, BD);

                return ((0.0 <= dotABAM) && (dotABAM <= dotABAB) && (0.0 <= dotBDBM) && (dotBDBM <= dotBDBD));
            }


            void describe()
            {
                std::cout<<"A: "<< A_.x << " "<< A_.y << "  \n";
                std::cout<<"B: "<< B_.x << " "<< B_.y << "  \n";
                std::cout<<"C: "<< C_.x << " "<< C_.y << "  \n";
                std::cout<<"D: "<< D_.x << " "<< D_.y << "  \n";
                std::cout<<"------------------------------  \n";
            }


        private:
            void randomly_generate(const T & width, const FracNetSamplerParams2D<T> & params)
            {
                //unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
                // const unsigned seed = 3;
                static std::default_random_engine generator (params.seed);

                // this one needs to be replaced
                std::uniform_real_distribution<> distr_point_xdir(params.x_min, params.x_max);
                std::uniform_real_distribution<> distr_point_ydir(params.y_min, params.y_max);

                std::uniform_int_distribution<> distr_angle(params.alpha_min, params.alpha_max);

                A_.x = distr_point_xdir(generator);
                A_.y = distr_point_ydir(generator);

                T theta = distr_angle(generator);

                // length should be driven from power distribution
                // std::uniform_real_distribution<> distr_length(3.0*width, 0.15);

                // length should be driven from power distribution
                std::uniform_real_distribution<> distr_length(std::max(params.x_min, params.y_min), std::min(params.x_max, params.y_max));
                // const T x_min = 3.0*width > 0.04 ? 3.0*width : 0.04;
                

                const T x_min = (params.min_length==0)? 3.0*width : params.min_length; 
                const T r = distr_length(generator);
                T length = x_min * std::pow( (1.-r), (-1./(params.pow_dist_coef_length - 1.)));

                generate_rectangle(length, width, theta);
            }

            void generate_rectangle(const T & a, const T & b, const T & theta){
                const T pi = std::acos(-1.0);
                T theta_rad = theta * pi/180.0;

                B_.x = A_.x + (a * std::cos(theta_rad));
                B_.y = A_.y + (a * std::sin(theta_rad));

                C_.x = A_.x + (b * std::cos(theta_rad + pi/2.0));
                C_.y = A_.y + (b * std::sin(theta_rad + pi/2.0));

                D_.x = A_.x + ((a * std::cos(theta_rad)) + (b * std::cos(theta_rad + pi/2.0)));
                D_.y = A_.y + ((a * std::sin(theta_rad)) + (b * std::sin(theta_rad + pi/2.0)));
            }

            void build_vector(const Point2D<T> & A, const Point2D<T> & B, Point2D<T> & result)
            {
                result.x = B.x - A.x;
                result.y = B.y - A.y;
            }

            T vec_dot(const Point2D<T> & A, const Point2D<T> & B)
            {
                return (A.x * B.x) + (A.y * B.y);
            }

        private:
            Point2D<T> A_;
            Point2D<T> B_;
            Point2D<T> C_;
            Point2D<T> D_;
    };


    template<class FunctionSpace>
    class InitialCondidtionPFFracNet2D :public InitialCondition<FunctionSpace>
    {
        public:

            // using Comm           = typename FunctionSpace::Comm;
            using Mesh           = typename FunctionSpace::Mesh;
            using Elem           = typename FunctionSpace::Shape;
            using ElemView       = typename FunctionSpace::ViewDevice::Elem;
            using SizeType       = typename FunctionSpace::SizeType;
            using Scalar         = typename FunctionSpace::Scalar;
            using Dev            = typename FunctionSpace::Device;
            using Point          = typename FunctionSpace::Point;
            using ElemViewScalar = typename utopia::FunctionSpace<Mesh, 1, Elem>::ViewDevice::Elem;
            static const int NNodes = Elem::NNodes;

            InitialCondidtionPFFracNet2D(FunctionSpace &space, const SizeType & PF_component, const SizeType & num_fracs=10):
            InitialCondition<FunctionSpace>(space), PF_component_(PF_component), num_fracs_(num_fracs), pressure0_(1.0)
            {

            }

            void read(Input &in) override
            {
                in.get("num_fracs", num_fracs_);
                in.get("pressure0", pressure0_);

                sampler_params_.read(in); 
            }



            void init(PetscVector &x) override
            {
                using CoeffVector = utopia::StaticVector<Scalar, NNodes>;
                // un-hard-code
                auto C = this->space_.subspace(PF_component_);

                auto width =  3.0 * this->space_.mesh().min_spacing();

                // if(mpi_world_rank()==0){
                //     std::cout<<"width: "<< width << "  \n";
                // }

                std::vector<Rectangle<Scalar>> rectangles;

                for(auto r=0; r < num_fracs_; r++){
                    rectangles.push_back(Rectangle<Scalar>(width, sampler_params_));
                }

                auto sampler = utopia::sampler(C, [&rectangles](const Point &x) -> Scalar {

                    for(auto r=0; r < rectangles.size(); r++){
                        if(rectangles[r].belongs_to_rectangle(x[0], x[1]))
                            return 1.0;
                    }
                    return 0.0;
                });

                {
                    auto C_view       = C.view_device();
                    auto sampler_view = sampler.view_device();
                    auto x_view       = this->space_.assembly_view_device(x);

                    Dev::parallel_for(this->space_.element_range(), UTOPIA_LAMBDA(const SizeType &i) {
                        ElemViewScalar e;
                        C_view.elem(i, e);

                        CoeffVector s;
                        sampler_view.assemble(e, s);
                        C_view.set_vector(e, s, x_view);
                    });
                }
            }

            void init(PetscVector &sol_vec, PetscVector &press_vec) override
            {
                PetscVector sol_vec_copy = sol_vec;

                // un-hard-code
                auto C = this->space_.subspace(PF_component_);
                auto width =  5.0 * this->space_.mesh().min_spacing();

                // if(mpi_world_rank()==0){
                //     std::cout<<"width: "<< width << "  \n";
                // }

                const Point2D<Scalar> A(0.5, 0.5);
                Rectangle<Scalar> rectangle(A, width, width, 0.0);

                auto sampler = utopia::sampler(C, [&rectangle](const Point &x) -> Scalar {
                    if(rectangle.belongs_to_rectangle(x[0], x[1])){
                        return 1.0;
                    }
                    else{
                        return 0.0;
                    }
                });


                Scalar p = pressure0_;

                auto press_sampler = utopia::sampler(C, [&rectangle, p](const Point &x) -> Scalar {
                    if(rectangle.belongs_to_rectangle(x[0], x[1])){
                        return p;
                    }
                    else{
                        return p/10.0;
                    }
                });


                {
                    auto C_view             = C.view_device();
                    auto sampler_view       = sampler.view_device();
                    auto press_sampler_view = press_sampler.view_device();

                    auto sol_view       = this->space_.assembly_view_device(sol_vec);
                    auto press_view     = this->space_.assembly_view_device(press_vec);

                    Dev::parallel_for(this->space_.element_range(), UTOPIA_LAMBDA(const SizeType &i) {
                        ElemViewScalar e;
                        C_view.elem(i, e);

                        StaticVector<Scalar, NNodes> s;
                        sampler_view.assemble(e, s);
                        C_view.set_vector(e, s, sol_view);

                        press_sampler_view.assemble(e, s);
                        C_view.set_vector(e, s, press_view);

                    });
                }

                // add new fracture to existing ones
                sol_vec += sol_vec_copy;
            }

        private:
            SizeType PF_component_;
            SizeType num_fracs_;
            Scalar pressure0_;

            FracNetSamplerParams2D<Scalar>  sampler_params_; 

    };

}

#endif

#endif
