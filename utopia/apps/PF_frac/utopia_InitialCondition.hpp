#ifndef UTOPIA_INITIAL_CONDITION_PF_HPP
#define UTOPIA_INITIAL_CONDITION_PF_HPP

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

#include <random>
#include <cmath>
#include <chrono>

namespace utopia {

    template<class FunctionSpace>
    class InitialCondition : public Configurable
    {

        public:
            InitialCondition(FunctionSpace &space): space_(space){ 

            }

            virtual ~InitialCondition(){

            }

            void read(Input &in) override
            {

            }                  

        virtual void init(PetscVector &x) = 0; 


        protected:
            FunctionSpace & space_; 
    };


    template<class FunctionSpace>
    class InitialCondidtionPFTension : public InitialCondition<FunctionSpace>
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


            InitialCondidtionPFTension(FunctionSpace &space, const SizeType & PF_component):    InitialCondition<FunctionSpace>(space), 
                                                                                                PF_component_(PF_component) 
            {

            }


            void init(PetscVector &x) override
            {
                // un-hard-code
                auto C = this->space_.subspace(PF_component_);

                auto sampler = utopia::sampler(C, UTOPIA_LAMBDA(const Point &x) -> Scalar {
                    Scalar f = 0.0;
                        if(  x[0] > (0.5-this->space_.mesh().min_spacing()) && x[0] < (0.5 + this->space_.mesh().min_spacing())  && x[1]  < 0.5 ){
                            f = 1.0; 
                        }
                        else{
                            f = 0.0; 
                        }
                    return f;
                });

                {
                    auto C_view       = C.view_device();
                    auto sampler_view = sampler.view_device();
                    auto x_view       = this->space_.assembly_view_device(x);

                    Dev::parallel_for(this->space_.local_element_range(), UTOPIA_LAMBDA(const SizeType &i) {
                        ElemViewScalar e;
                        C_view.elem(i, e);

                        StaticVector<Scalar, NNodes> s;
                        sampler_view.assemble(e, s);
                        C_view.set_vector(e, s, x_view);
                    });
                }
            }

        private:
            SizeType PF_component_; 
    };


    template<class FunctionSpace>
    class InitialCondidtionPFTbar : public InitialCondition<FunctionSpace>
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


            InitialCondidtionPFTbar(FunctionSpace &space, const SizeType & PF_component):   InitialCondition<FunctionSpace>(space), 
                                                                                            PF_component_(PF_component) 
            {

            }


            void init(PetscVector &x) override
            {
                // un-hard-code
                auto C = this->space_.subspace(PF_component_);

                auto sampler = utopia::sampler(C, UTOPIA_LAMBDA(const Point &x) -> Scalar {
                    Scalar f = 0.0;
                        if(  x[0] > (0.5-this->space_.mesh().min_spacing()) && x[0] < (0.5 + this->space_.mesh().min_spacing())  && x[1]  < 0.5  && x[1]  > 0.3  ){
                            f = 1.0; 
                        }
                        else if((x[0] > 0.3) && (x[0] < 0.7) && (x[1]  > 0.7 - this->space_.mesh().min_spacing())  &&  (x[1]  < 0.7 + this->space_.mesh().min_spacing())){
                            f = 1.0; 
                        }
                        else{
                            f = 0.0; 
                        }

                    return f;
                });

                {
                    auto C_view       = C.view_device();
                    auto sampler_view = sampler.view_device();
                    auto x_view       = this->space_.assembly_view_device(x);

                    Dev::parallel_for(this->space_.local_element_range(), UTOPIA_LAMBDA(const SizeType &i) {
                        ElemViewScalar e;
                        C_view.elem(i, e);

                        StaticVector<Scalar, NNodes> s;
                        sampler_view.assemble(e, s);
                        C_view.set_vector(e, s, x_view);
                    });
                }
            }

        private:
            SizeType PF_component_; 
    }; 


    template<class T>
    struct Point2D{
        T x; 
        T y; 
        void describe()
        {
            std::cout<<"(" << x << " , "<< y << " ) \n"; 
        }
    }; 

    template<class T>
    class Rectangle
    {
        public:
            Rectangle(const Point2D<T> & A, const Point2D<T> & B, const Point2D<T> & C, const Point2D<T> & D): 
            A_(A), B_(B), C_(C), D_(D)
            {

            }

            Rectangle(const T & width)
            {
                randomly_generate(width); 
            }        
            
            Rectangle(const Point2D<T> & A, const T & length, const T & width, const T & theta ): 
            A_(A)
            {
                this->generate_rectangle(length, width, theta); 
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
            void randomly_generate(const T & width)
            {
                //unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
                const unsigned seed = 3; 
                static std::default_random_engine generator (seed);

                // this one needs to be replaced 
                std::uniform_real_distribution<> distr_point(0.0, 1.0);                 
                std::uniform_int_distribution<> distr_angle(0.0, 180); 

                // this one should be driven from power distribution 
                std::uniform_real_distribution<> distr_length(2.0*width, 0.4);                 

                A_.x = distr_point(generator); 
                A_.y = distr_point(generator); 

                T length = distr_length(generator); 
                T theta = distr_angle(generator); 

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
    class InitialCondidtionPFFracNet : public InitialCondition<FunctionSpace>
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


            InitialCondidtionPFFracNet(FunctionSpace &space, const SizeType & PF_component, const SizeType & num_fracs=10):    
            InitialCondition<FunctionSpace>(space), PF_component_(PF_component), num_fracs_(num_fracs)
            {

            }

            void read(Input &in) override
            {
                in.get("num_fracs", num_fracs_);
            }                


            void init(PetscVector &x) override
            {
                // un-hard-code
                auto C = this->space_.subspace(PF_component_);

                auto width =  3.0 * this->space_.mesh().min_spacing(); 

                if(mpi_world_rank()==0){
                    std::cout<<"width: "<< width << "  \n"; 
                }

                std::vector<Rectangle<Scalar>> rectangles; 

                for(auto r=0; r < num_fracs_; r++){
                    rectangles.push_back(Rectangle<Scalar>(width)); 
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

                    Dev::parallel_for(this->space_.local_element_range(), UTOPIA_LAMBDA(const SizeType &i) {
                        ElemViewScalar e;
                        C_view.elem(i, e);

                        StaticVector<Scalar, NNodes> s;
                        sampler_view.assemble(e, s);
                        C_view.set_vector(e, s, x_view);
                    });
                }
            }

        private:
            SizeType PF_component_; 
            SizeType num_fracs_; 
    };

}

#endif