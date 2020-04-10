#ifndef UTOPIA_INITIAL_CONDITION_PF_HPP
#define UTOPIA_INITIAL_CONDITION_PF_HPP


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
        virtual void init(PetscVector &sol_vec, PetscVector & /*press_vec*/ ){
            this->init(sol_vec);
        }


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

                    Dev::parallel_for(this->space_.element_range(), UTOPIA_LAMBDA(const SizeType &i) {
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

                    Dev::parallel_for(this->space_.element_range(), UTOPIA_LAMBDA(const SizeType &i) {
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
    struct Point3D{

        Point3D()
        {

        }

        Point3D(const T & xx, const T & yy, const T && zz):
        x(xx), y(yy), z(zz)
        {

        }

        T x;
        T y;
        T z;

        void describe()
        {
            std::cout<<"( " << x << " , "<< y << " , "<< z << " )" <<  " \n";
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
                std::uniform_real_distribution<> distr_point(0.1, 0.9);
                std::uniform_int_distribution<> distr_angle(0.0, 180);


                A_.x = distr_point(generator);
                A_.y = distr_point(generator);

                T theta = distr_angle(generator);

                // length should be driven from power distribution
                // std::uniform_real_distribution<> distr_length(3.0*width, 0.15);

                // length should be driven from power distribution
                std::uniform_real_distribution<> distr_length(0.0, 1.0);
                const T alpha = 2.8;
                const T x_min = 3.0*width > 0.04 ? 3.0*width : 0.04;
                const T r = distr_length(generator);
                T length = x_min * std::pow( (1.-r), (-1./(alpha-1.)));

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


    template<class T>
    class Paralleloid
    {
        public:
            Paralleloid(const T & width)
            {
                points_.resize(8);
                randomly_generate(width);
            }

            Paralleloid(const Point3D<T> & A, const T & length, const T & width, const T & theta, const T & gamma )
            {
                points_.resize(8);
                points_[0] = A;
                this->generate_paralleloid(length, width, theta, gamma);
            }

            bool belongs_to_paralleloid(const T & x_coord, const T & y_coord, const T & z_coord)
            {
                Point3D<T> M;
                M.x = x_coord;
                M.y = y_coord;
                M.z = z_coord;

                return belongs_to_paralleloid(M);
            }

            bool belongs_to_paralleloid(Point3D<T> M)
            {

                Point3D<T> u, v, w;
                build_vector(points_[0], points_[4], u);
                build_vector(points_[0], points_[1], v);
                build_vector(points_[0], points_[2], w);

                const T dotMu = vec_dot(M, u);
                const T dotMv = vec_dot(M, v);
                const T dotMw = vec_dot(M, w);


                bool flg1 =  (vec_dot(u, points_[4]) <= dotMu) && (dotMu <= vec_dot(u, points_[0]));
                bool flg2 =  (vec_dot(v, points_[1]) <= dotMv) && (dotMv <= vec_dot(v, points_[0]));
                bool flg3 =  (vec_dot(w, points_[2]) <= dotMw) && (dotMw <= vec_dot(w, points_[0]));

                return flg1 && flg2 && flg3;

            }


            void describe()
            {
                std::cout<<"A: "<< points_[0].x << " "<< points_[0].y << " " << points_[0].z  <<" \n";
                std::cout<<"B: "<< points_[1].x << " "<< points_[1].y << " " << points_[1].z  <<" \n";
                std::cout<<"C: "<< points_[2].x << " "<< points_[2].y << " " << points_[2].z  <<" \n";
                std::cout<<"D: "<< points_[3].x << " "<< points_[3].y << " " << points_[3].z  <<" \n \n";

                std::cout<<"E: "<< points_[4].x << " "<< points_[4].y <<  "  "<< points_[4].z << "  \n";
                std::cout<<"F: "<< points_[5].x << " "<< points_[5].y <<  "  "<< points_[5].z << "  \n";
                std::cout<<"G: "<< points_[6].x << " "<< points_[6].y <<  "  "<< points_[6].z << "  \n";
                std::cout<<"H: "<< points_[7].x << " "<< points_[7].y <<  "  "<< points_[7].z << "  \n";
                std::cout<<"------------------------------  \n";
            }

        private:
            void randomly_generate(const T & depth)
            {
                //unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
                const unsigned seed = 15;
                static std::default_random_engine generator (seed);

                // this one needs to be replaced
                std::uniform_real_distribution<> distr_point_x(4.0, 6.0);
                std::uniform_real_distribution<> distr_point_y(4.0, 6.0);
                std::uniform_real_distribution<> distr_point_z(4.0, 6.0);
                std::uniform_int_distribution<> distr_angle(0, 180);


                points_[0].x = distr_point_x(generator);
                points_[0].y = distr_point_y(generator);
                points_[0].z = distr_point_z(generator);


                const T theta = distr_angle(generator);
                const T gamma = distr_angle(generator);

                // // length should be driven from power distribution
                // std::uniform_real_distribution<> distr_length(0.0, 1.0);
                // const T alpha1 = 2.8;
                // const T x_min = 3.0*depth > 0.04 ? 3.0*depth : 0.04;
                // const T r = distr_length(generator);
                // const T length = x_min * std::pow( (1.-r), (-1./(alpha1-1.)));


                // const T alpha2 = 2.95;
                // const T r2 = distr_length(generator);
                // const T width = x_min * std::pow( (1.-r2), (-1./(alpha2-1.)));


                const T width = 0.15;
                const T length = 0.15;

                std::uniform_int_distribution<> distr_dir(0.0, 6);
                const int i = distr_dir(generator);

                if(i==0){
                    generate_paralleloid(length, depth, width, theta, gamma);
                }
                else if(i==1){
                    generate_paralleloid(length, width, depth, theta, gamma);
                }
                else if(i==3){
                    generate_paralleloid(depth, length, width, theta, gamma);
                }
                else if(i==4){
                    generate_paralleloid(depth, width, length, theta, gamma);
                }
                else if(i==5){
                    generate_paralleloid(length, depth, width, theta, gamma);
                }
                else{
                    generate_paralleloid(width, depth, length, theta, gamma);
                }


            }

            void generate_paralleloid(const T & a, const T & b, const T & c, const T & theta, const T & gamma){
                const T pi = std::acos(-1.0);
                const T theta_rad = theta * pi/180.0;
                const T gamma_rad = gamma * pi/180.0;

                const auto cos_theta = std::cos(theta_rad);
                const auto sin_theta = std::sin(theta_rad);

                const auto cos_theta_pi = std::cos(theta_rad + pi/2.0);
                const auto sin_theta_pi = std::sin(theta_rad + pi/2.0);


                points_[1].x = points_[0].x + (a * cos_theta);
                points_[1].y = points_[0].y + (a * sin_theta);
                points_[1].z = points_[0].z;

                points_[2].x = points_[0].x + (b * cos_theta_pi);
                points_[2].y = points_[0].y + (b * sin_theta_pi);
                points_[2].z = points_[0].z;

                points_[3].x = points_[0].x + ((a * cos_theta) + (b * cos_theta_pi));
                points_[3].y = points_[0].y + ((a * sin_theta) + (b * sin_theta_pi));
                points_[3].z = points_[0].z;

                points_[4].x = points_[0].x;
                points_[4].y = points_[0].y;
                points_[4].z = points_[0].z + c;

                points_[5].x = points_[1].x;
                points_[5].y = points_[1].y;
                points_[5].z = points_[1].z + c;

                points_[6].x = points_[2].x;
                points_[6].y = points_[2].y;
                points_[6].z = points_[2].z + c;

                points_[7].x = points_[3].x;
                points_[7].y = points_[3].y;
                points_[7].z = points_[3].z + c;


                // // rotate 
                // const auto cos_gamma = std::cos(gamma_rad); 
                // const auto sin_gamma = std::sin(gamma_rad); 

                // for(auto p=0; p < points_.size(); p++){
                //     auto ax = points_[p].x; 
                //     auto ay = points_[p].z; 
                //     points_[p].x = (ax* cos_gamma) - (ay * sin_gamma); 
                //     points_[p].z = (ay* cos_gamma) + (ax * sin_gamma); 
                // }

            }

            void build_vector(const Point3D<T> & A, const Point3D<T> & B, Point3D<T> & result)
            {
                result.x = A.x - B.x;
                result.y = A.y - B.y;
                result.z = A.z - B.z;
            }

            T vec_dot(const Point3D<T> & A, const Point3D<T> & B)
            {
                return (A.x * B.x) + (A.y * B.y) + (A.z * B.z);
            }


            void vec_cross(const Point3D<T> & A, const Point3D<T> & B, Point3D<T> & result)
            {
                result.x = A.y * B.z - A.z*B.y;
                result.y = A.z * B.x - A.x*B.z;
                result.z = A.x * B.y - A.y*B.x;
            }

        private:
            std::vector< Point3D<T> > points_;

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
            InitialCondition<FunctionSpace>(space), PF_component_(PF_component), num_fracs_(num_fracs), pressure0_(1.0)
            {

            }

            void read(Input &in) override
            {
                in.get("num_fracs", num_fracs_);
                in.get("pressure0", pressure0_);
            }



            void init(PetscVector &x) override
            {
                using CoeffVector = utopia::StaticVector<Scalar, NNodes>;
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

                if(mpi_world_rank()==0){
                    std::cout<<"width: "<< width << "  \n";
                }

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

    };


    template<class FunctionSpace>
    class InitialCondidtionPFFracNet3D : public InitialCondition<FunctionSpace>
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

            InitialCondidtionPFFracNet3D(FunctionSpace &space, const SizeType & PF_component, const SizeType & num_fracs=10):
            InitialCondition<FunctionSpace>(space), PF_component_(PF_component), num_fracs_(num_fracs)//, pressure0_(1.0)
            {

            }

            void read(Input &in) override
            {
                in.get("num_fracs", num_fracs_);
                // in.get("pressure0", pressure0_);
            }



            void init(PetscVector &x) override
            {
                using CoeffVector = utopia::StaticVector<Scalar, NNodes>;
                // un-hard-code
                auto C = this->space_.subspace(PF_component_);

                auto width =  3.0 * this->space_.mesh().min_spacing();

                if(mpi_world_rank()==0){
                    std::cout<<"width: "<< width << "  \n";
                }

                std::vector<Paralleloid<Scalar>> paralleloids;

                for(auto r=0; r < num_fracs_; r++){
                    paralleloids.push_back(Paralleloid<Scalar>(width));
                }

                auto sampler = utopia::sampler(C, [&paralleloids](const Point &x) -> Scalar {

                    for(auto r=0; r < paralleloids.size(); r++){
                        if(paralleloids[r].belongs_to_paralleloid(x[0], x[1], x[2]))
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

        private:
            SizeType PF_component_;
            SizeType num_fracs_;
    };


    template<class FunctionSpace>
    class InitialCondidtionPFSneddon : public InitialCondition<FunctionSpace>
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


            InitialCondidtionPFSneddon(FunctionSpace &space, const SizeType & PF_component):   InitialCondition<FunctionSpace>(space),
                                                                                                PF_component_(PF_component)
            {

            }


            void init(PetscVector &x) override
            {
                // un-hard-code
                auto C = this->space_.subspace(PF_component_);

                auto sampler = utopia::sampler(C, UTOPIA_LAMBDA(const Point &x) -> Scalar {
                    Scalar f = 0.0;
                    Scalar h = this->space_.mesh().min_spacing();

                    if( x[0] > 0.55 && x[0] < 0.7  && x[1]  > (0.4)  && x[1]  < (0.6) && x[2] < (0.6+h) && x[2] > (0.6-h)){
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

                    Dev::parallel_for(this->space_.element_range(), UTOPIA_LAMBDA(const SizeType &i) {
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
    class Mixed: public InitialCondition<FunctionSpace>
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

            Mixed(FunctionSpace &space, const SizeType & PF_component):
            InitialCondition<FunctionSpace>(space), PF_component_(PF_component), pressure0_(1.0)
            {

            }

            void init(PetscVector &x) override
            {
                using CoeffVector = utopia::StaticVector<Scalar, NNodes>;
                // un-hard-code
                auto C = this->space_.subspace(PF_component_);

                auto width =  3.0 * this->space_.mesh().min_spacing();
                // auto width = 0.1;

                if(mpi_world_rank()==0){
                    std::cout<<"width: "<< width << "  \n";
                }

                std::vector<Rectangle<Scalar>> rectangles;


                Point2D<Scalar> A;
                A.x = 38.1;
                A.y = 82.55;
                rectangles.push_back(Rectangle<Scalar>(A, -12.7, width, -45));


                Point2D<Scalar> B;
                B.x = 38.1;
                B.y = 69.85;
                rectangles.push_back(Rectangle<Scalar>(B, 12.7, width, -45));



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


        private:
            SizeType PF_component_;
            SizeType num_fracs_;
            Scalar pressure0_;

    };
    
    
    
    template<class FunctionSpace>
    class AsphaltTension: public InitialCondition<FunctionSpace>
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
        
        AsphaltTension(FunctionSpace &space, const SizeType & PF_component):
        InitialCondition<FunctionSpace>(space), PF_component_(PF_component), pressure0_(1.0)
        {
            
        }
        
        void init(PetscVector &x) override
        {
            using CoeffVector = utopia::StaticVector<Scalar, NNodes>;
            // un-hard-code
            auto C = this->space_.subspace(PF_component_);
            
            auto width =  3.0 * this->space_.mesh().min_spacing();
            // auto width = 0.1;
            
            if(mpi_world_rank()==0){
                std::cout<<"width: "<< width << "  \n";
            }
            
            std::vector<Rectangle<Scalar>> rectangles;
            
            
            Point2D<Scalar> A;
            A.x = 12.00;
            A.y = 21.50;
            rectangles.push_back(Rectangle<Scalar>(A, -3.5355, width, 45));
            
            
            Point2D<Scalar> B;
            B.x = 15;
            B.y = 20.00;
            rectangles.push_back(Rectangle<Scalar>(B, 5.000, width, 0.0));
            
            
            
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
        
        
    private:
        SizeType PF_component_;
        SizeType num_fracs_;
        Scalar pressure0_;
        
    };



}

#endif

#endif
