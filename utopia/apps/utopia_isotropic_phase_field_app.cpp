
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

#include <random>
#include <cmath>
#include <chrono>

namespace utopia {

    template<class FunctionSpace>
    static void init_phase_field_tension(FunctionSpace &space, PetscVector &x)
    {
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

        // un-hard-code
        auto C = space.subspace(0);

        auto sampler = utopia::sampler(C, UTOPIA_LAMBDA(const Point &x) -> Scalar {
            Scalar f = 0.0;
                if(  x[0] > (0.5-space.mesh().min_spacing()) && x[0] < (0.5 + space.mesh().min_spacing())  && x[1]  < 0.5 ){
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
            auto x_view       = space.assembly_view_device(x);

            Dev::parallel_for(space.local_element_range(), UTOPIA_LAMBDA(const SizeType &i) {
                ElemViewScalar e;
                C_view.elem(i, e);

                StaticVector<Scalar, NNodes> s;
                sampler_view.assemble(e, s);
                C_view.set_vector(e, s, x_view);
            });
        }
    }

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
                unsigned seed = 1.0; 
                static std::default_random_engine generator (seed);

                // this one needs to be replaced 
                std::uniform_real_distribution<> distr_point(0.0, 1.0);                 
                std::uniform_int_distribution<> distr_angle(0.0, 180); 
                std::uniform_real_distribution<> distr_length(0.0, 1.0);                 

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


    // this code does not run in parallel correctly 
    template<class FunctionSpace>
    static void init_phase_field_frac_net(FunctionSpace &space, PetscVector &x)
    {
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

        // un-hard-code
        auto C = space.subspace(0);

        auto width =  3.0 * space.mesh().min_spacing(); 

        std::cout<<"width: "<< width << "  \n"; 
        std::vector<Rectangle<Scalar>> rectangles; 

        for(auto r=0; r < 30; r++){
            rectangles.push_back(Rectangle<Scalar>(width)); 
        }

        auto sampler = utopia::sampler(C, [&rectangles](const Point &x) -> Scalar {

            for(auto r=0; r< rectangles.size(); r++){
                if(rectangles[r].belongs_to_rectangle(x[0], x[1]))
                    return 1.0; 
            }
            return 0.0; 

            // return 1.0; 

        });

        {
            auto C_view       = C.view_device();
            auto sampler_view = sampler.view_device();
            auto x_view       = space.assembly_view_device(x);

            Dev::parallel_for(space.local_element_range(), UTOPIA_LAMBDA(const SizeType &i) {
                ElemViewScalar e;
                C_view.elem(i, e);

                StaticVector<Scalar, NNodes> s;
                sampler_view.assemble(e, s);
                C_view.set_vector(e, s, x_view);
            });
        }
    }




    template<class FunctionSpace>
    static void enforce_BC_time_dependent(FunctionSpace &space, const typename FunctionSpace::Scalar & disp,  const typename FunctionSpace::Scalar & t){

        static const int Dim   = FunctionSpace::Dim;

        using Point          = typename FunctionSpace::Point;
        using Scalar         = typename FunctionSpace::Scalar;

        space.emplace_dirichlet_condition(
            SideSet::left(),
            UTOPIA_LAMBDA(const Point &p) -> Scalar {
                return 0.0;
            },
            1
            );

        space.emplace_dirichlet_condition(
            SideSet::right(),
            UTOPIA_LAMBDA(const Point &p) -> Scalar {
                return t*disp;
            },
            1
            );

        for(int d = 2; d < Dim + 1; ++d) {
            space.emplace_dirichlet_condition(
                SideSet::left(),
                UTOPIA_LAMBDA(const Point &p) -> Scalar {
                return 0.0;
                },
                d
            );

            space.emplace_dirichlet_condition(
                SideSet::right(),
                UTOPIA_LAMBDA(const Point &p) -> Scalar {
                    return 0.0;
                },
                d
                );
            }
    }

    template<class FunctionSpace>
    static void enforce_BC_all_disp_fixed(FunctionSpace &space){

        static const int Dim   = FunctionSpace::Dim;

        using Point          = typename FunctionSpace::Point;
        using Scalar         = typename FunctionSpace::Scalar;

        for(int d = 1; d < Dim + 1; ++d) {
            space.emplace_dirichlet_condition(
                SideSet::left(),
                UTOPIA_LAMBDA(const Point &p) -> Scalar {
                return 0.0;
                },
                d
            );

            space.emplace_dirichlet_condition(
                SideSet::right(),
                UTOPIA_LAMBDA(const Point &p) -> Scalar {
                    return 0.0;
                },
                d
                );

            space.emplace_dirichlet_condition(
                SideSet::top(),
                UTOPIA_LAMBDA(const Point &p) -> Scalar {
                    return 0.0;
                },
                d
                ); 

            space.emplace_dirichlet_condition(
                SideSet::bottom(),
                UTOPIA_LAMBDA(const Point &p) -> Scalar {
                    return 0.0;
                },
                d
                );
        }
    }


    
    template<class FunctionSpace>
    static void build_irreversility_constraint(const PetscVector &x_old, PetscVector &x_new, const typename FunctionSpace::SizeType & comp)
    {   
        static const int Dim   = FunctionSpace::Dim;
        using Scalar         = typename FunctionSpace::Scalar;
        using SizeType       = typename FunctionSpace::SizeType;

        {
            auto d_x_old = const_device_view(x_old);

            parallel_transform(x_new, UTOPIA_LAMBDA(const SizeType &i, const Scalar &xi) -> Scalar 
            {
                if(i%(Dim+1)==comp){
                    return d_x_old.get(i); 
                }
                else{
                    return -9e15; 
                }                    
            });
        }
    }    


    template<class FunctionSpace>
    static void isotropic_phase_field_fracture_sim(
        FunctionSpace &space,
        MPITimeStatistics &stats,
        Input &in)
    {
        static const int Dim   = FunctionSpace::Dim;
        static const int NVars = FunctionSpace::Dim + 1;

        //expose inner types
        using Comm           = typename FunctionSpace::Comm;
        using Mesh           = typename FunctionSpace::Mesh;
        using Elem           = typename FunctionSpace::Shape;
        using ElemView       = typename FunctionSpace::ViewDevice::Elem;
        using SizeType       = typename FunctionSpace::SizeType;
        using Scalar         = typename FunctionSpace::Scalar;
        using Dev            = typename FunctionSpace::Device;
        using Point          = typename FunctionSpace::Point;
        using ElemViewScalar = typename utopia::FunctionSpace<Mesh, 1, Elem>::ViewDevice::Elem;

        static const int NNodes = Elem::NNodes;

        using FEFunction     = utopia::FEFunction<FunctionSpace>;
        using Quadrature     = utopia::Quadrature<Elem, 2>;
        using Parameters     = typename IsotropicPhaseFieldForBrittleFractures<FunctionSpace>::Parameters;

        auto &mesh = space.mesh();

        Scalar disp = 0.001;

        in.get("disp", disp);

        bool with_damage = true;
        in.get("with_damage", with_damage);

        bool with_BC = true;
        in.get("with_BC", with_BC);

        bool pressure_test = false;
        in.get("pressure_test", pressure_test);      


        bool tension_test = false;
        in.get("tension_test", tension_test);     


        std::cout<<"pressure_test: "<< pressure_test << "   \n"; 
        std::cout<<"tension_test: "<< tension_test << "   \n"; 


        Scalar dt = 1e-5;
        in.get("dt", dt);   

        Scalar num_time_steps = 1000;  
        in.get("num_time_steps", num_time_steps);   


        stats.start();


        stats.stop_collect_and_restart("BC");

        Scalar pressure0; 
        
        PetscVector x;
        space.create_vector(x);
        x.set(0.0);

        if(with_damage) {

            if(tension_test){
                init_phase_field_tension(space, x); 
            }

            if(pressure_test){
                init_phase_field_frac_net(space, x); 
                in.get("pressure", pressure0);                
            }
        }
        
        stats.stop_collect_and_restart("phase-field-init");



        Scalar time_= dt; 
        std::string output_path = "isotropic_PFfrac_test";
        in.get("output-path", output_path);
        // print IG 
        rename("X", x);

        PetscVector irreversibility_constraint = x; 

        // as space gets copied, we need to instantiate PF problem every time BC changes ... 
        IsotropicPhaseFieldForBrittleFractures<FunctionSpace> pp(space);
        pp.read(in);     

        // testing IC ... 
        space.apply_constraints(x);                       
        space.write(output_path+"_"+std::to_string(0.0)+".vtk", x);     
        // exit(0); 

        if(with_BC) {
            if(pressure_test){
                enforce_BC_all_disp_fixed(space); 
            }        
        }

        for (auto t=1; t < num_time_steps; t++)
        {
            std::cout<<"Time-step: "<< t << "  \n"; 
     
            if(with_BC) {
                if(tension_test){
                    space.reset_bc(); 
                    enforce_BC_time_dependent(space, disp, time_);              
                }
            }

            space.apply_constraints(x);    

            std::cout<<"pressure: "<< pressure0 * time_ << "  \n"; 
            pp.set_pressure(pressure0 * time_);     

                                                                                        // PF component=0
            build_irreversility_constraint<FunctionSpace>(x, irreversibility_constraint, 0); 


            // auto linear_solver = std::make_shared<Factorization<PetscMatrix, PetscVector>>();
            // auto linear_solver = std::make_shared<GMRES<PetscMatrix, PetscVector>>();
            // // linear_solver->max_it(200); 
            // linear_solver->pc_type("jacobi");
            // Newton<PetscMatrix, PetscVector> solver(linear_solver);
            // in.get("solver", solver);

            // MPRGP sucks as a solver, as it can not be preconditioned easily ... 
            auto qp_solver = std::make_shared<utopia::MPGRP<PetscMatrix, PetscVector> >();
            qp_solver->max_it(1000); 

            // tao seems to be even slower... 
            // auto linear_solver = std::make_shared<GMRES<PetscMatrix, PetscVector>>();
            // // linear_solver->max_it(200); 
            // linear_solver->pc_type("jacobi");            
            // auto qp_solver =  std::make_shared<utopia::TaoQPSolver<PetscMatrix, PetscVector> >(linear_solver);
            TrustRegionVariableBound<PetscMatrix, PetscVector> solver(qp_solver);
            auto box = make_lower_bound_constraints(make_ref(irreversibility_constraint));
            solver.set_box_constraints(box);
            in.get("solver", solver);
            solver.solve(pp, x);
            

            // auto hessian_approx   = std::make_shared<LBFGS<PetscVector> >(15);
            // auto qp_solver = std::make_shared<MPGRP<PetscMatrix, PetscVector> >();
            // qp_solver->max_it(5000); 
            
            // // // PetscMatrix H; 
            // // // std::function< void(const PetscVector &, PetscVector &) > H0_action_fun =
            // // // [&pp, &H](const PetscVector &x, PetscVector & result){ 
            // // //     pp.hessian(x, H); 
            // // //     result = H*x; 
            // // // };
            // // // hessian_approx->H0_action(H0_action_fun); 

            // QuasiTrustRegionVariableBound<PetscVector> solver(hessian_approx, qp_solver);
            // in.get("solver", solver);
            // solver.solve(pp, x);


            rename("X", x);
            space.write(output_path+"_"+std::to_string(time_)+".vtk", x);       

            // increment time step 
            time_+=dt;
        }
 
        stats.stop_collect_and_restart("solve+assemble");

        stats.stop_and_collect("output");
        stats.describe(std::cout);

    }

    static void petsc_tension_isotropic_phase_field_2(Input &in)
    {
        static const int Dim = 2;
        static const int NVars = Dim + 1;

        using Comm           = utopia::PetscCommunicator;
        using Mesh           = utopia::PetscDM<Dim>;
        using Elem           = utopia::PetscUniformQuad4;
        using FunctionSpace  = utopia::FunctionSpace<Mesh, NVars, Elem>;
        using SizeType       = FunctionSpace::SizeType;

        Comm world;

        MPITimeStatistics stats(world);
        stats.start();

        SizeType scale = (world.size() + 1);
        SizeType nx = scale * 4;
        SizeType ny = scale * 4;

        in.get("nx", nx);
        in.get("ny", ny);

        FunctionSpace space;

        space.build(
            world,
            {nx, ny},
            {0.0, 0.0},
            {1.0, 1.0}
            );

        space.mesh().set_field_name(0, "c");
        space.mesh().set_field_name(1, "disp_x");
        space.mesh().set_field_name(2, "disp_y");

        stats.stop_and_collect("space-creation");

        isotropic_phase_field_fracture_sim(
            space,
            stats,
            in
        );
    }

    UTOPIA_REGISTER_APP(petsc_tension_isotropic_phase_field_2);

    static void petsc_tension_isotropic_phase_field_3(Input &in)
    {
        static const int Dim = 3;
        static const int NVars = Dim + 1;

        using Comm           = utopia::PetscCommunicator;
        using Mesh           = utopia::PetscDM<Dim>;
        using Elem           = utopia::PetscUniformHex8;
        using FunctionSpace  = utopia::FunctionSpace<Mesh, NVars, Elem>;
        using SizeType       = FunctionSpace::SizeType;

        Comm world;

        MPITimeStatistics stats(world);
        stats.start();

        SizeType scale = (world.size() + 1);
        SizeType nx = scale * 4;
        SizeType ny = scale * 4;
        SizeType nz = scale * 4;

        in.get("nx", nx);
        in.get("ny", ny);
        in.get("nz", nz);

        FunctionSpace space;

        space.build(
            world,
            {nx, ny, nz},
            {0.0, 0.0, 0.0},
            {1.0, 1.0, 1.0}
            );

        space.mesh().set_field_name(0, "c");
        space.mesh().set_field_name(1, "disp_x");
        space.mesh().set_field_name(2, "disp_y");
        space.mesh().set_field_name(3, "disp_z");

        stats.stop_and_collect("space-creation");

        isotropic_phase_field_fracture_sim(
            space,
            stats,
            in);
    }

    UTOPIA_REGISTER_APP(petsc_tension_isotropic_phase_field_3);
}

