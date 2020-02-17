
#include "utopia_Base.hpp"

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
#include "utopia_PhaseField.hpp"
#include "utopia_FEFunction.hpp"
#include "utopia_SampleView.hpp"

#include <cmath>

namespace utopia {

    static void petsc_dm_multivar()
    {
        std::cout << "excuting: petsc_dm_multivar" << std::endl;
        static const int Dim = 2;
        static const int NVars = 2;

        using DMDA             = utopia::PetscDM<Dim>;
        using Elem             = utopia::PetscUniformQuad4;
        using FunctionSpace    = utopia::FunctionSpace<DMDA, NVars, Elem>;
        using SizeType         = DMDA::SizeType;
        using Scalar           = DMDA::Scalar;
        // using Point            = DMDA::Point;

        PetscCommunicator world;

        SizeType scale = (world.size() + 1);
        SizeType nx = scale * 2;
        SizeType ny = scale * 2;

        FunctionSpace space;
        space.build(
            world,
            {nx, ny},
            {0.0, 0.0},
            {1.0, 1.0}
        );

        PetscVector v;
        space.create_vector(v);

        each_write(v, [](const SizeType &i) -> Scalar {
            return i % NVars;
        });


        disp(space.mesh().n_nodes());
        disp(size(v));

        rename("U", v);
        space.write("prova.vtk", v);
    }

    UTOPIA_REGISTER_APP(petsc_dm_multivar);

    static void petsc_local_vec_view()
    {
        static const int Dim = 3;
        static const int NNodes = 8;

        using Mesh             = utopia::PetscDM<Dim>;
        using Elem             = utopia::PetscUniformHex8;
        using FunctionSpace    = utopia::FunctionSpace<Mesh, 1, Elem>;
        using SizeType         = Mesh::SizeType;
        // using Scalar           = Mesh::Scalar;
        // using Point            = Mesh::Point;

        PetscCommunicator world;

        SizeType scale = (world.size() + 1);
        SizeType nx = scale * 2;
        SizeType ny = scale * 2;
        SizeType nz = scale * 2;

        Mesh mesh(
            world,
            {nx, ny, nz},
            {0.0, 0.0, 0.0},
            {1.0, 1.0, 1.0}
        );

        FunctionSpace space(mesh);

        PetscVector vec;
        mesh.create_local_vector(vec);
        LocalViewDevice<const PetscVector, 1> view(vec);
    }

    UTOPIA_REGISTER_APP(petsc_local_vec_view);

    static void petsc_bratu(Input &in)
    {
        static const int Dim = 3;
        static const int NNodes = 8;

        using Mesh             = utopia::PetscDM<Dim>;
        using Elem             = utopia::PetscUniformHex8;
        using FunctionSpace    = utopia::FunctionSpace<Mesh, 1, Elem>;
        using SizeType         = Mesh::SizeType;
        using Scalar           = Mesh::Scalar;
        using Point            = Mesh::Point;

        PetscCommunicator world;

        MPITimeStatistics stats(world);

        ///////////////////////////////////////

        stats.start();

        SizeType scale = (world.size() + 1);
        SizeType nx = scale * 5;
        SizeType ny = scale * 5;
        SizeType nz = scale * 5;

        Scalar min_x = 0.0, max_x = 1.0;
        Scalar min_y = 0.0, max_y = 1.0;
        Scalar min_z = 0.0, max_z = 1.0;


        Scalar decay = 50.0;
        Scalar amplitude = 1.0;

        Point c1; c1.set(0.25);
        Point c2; c2.set(0.75);

        std::string output_path = "./fe.vtk";
        std::string model = "bratu";

        in.get("nx", nx);
        in.get("ny", ny);
        in.get("nz", nz);
        in.get("decay", decay);
        in.get("amplitude", amplitude);
        in.get("c1-x", c1[0]);
        in.get("c1-y", c1[1]);
        in.get("c1-z", c1[2]);

        in.get("c2-x", c2[0]);
        in.get("c2-y", c2[1]);
        in.get("c2-z", c2[2]);


        in.get("min-x", min_x);
        in.get("min-y", min_y);
        in.get("min-z", min_z);

        in.get("max-x", max_x);
        in.get("max-y", max_y);
        in.get("max-z", max_z);

        in.get("model", model);

        Mesh mesh(
            world,
            {nx, ny, nz},
            {min_x, min_y, min_z},
            {max_x, max_y, max_z}
        );

        stats.stop_and_collect("mesh-gen");

        ///////////////////////////////////////

        stats.start();

        FunctionSpace space(mesh);

        //boundary conditions
        space.emplace_dirichlet_condition(
            SideSet::left(),
            UTOPIA_LAMBDA(const Point &p) -> Scalar {
                return 0.0;
            }
        );

        space.emplace_dirichlet_condition(
            SideSet::right(),
            UTOPIA_LAMBDA(const Point &p) -> Scalar {
                return 0.0;
            }
        );

        stats.stop_and_collect("space+bc");
        ///////////////////////////////////////
        stats.start();

        std::shared_ptr<Function<PetscMatrix, PetscVector>> fun;

        auto ff = UTOPIA_LAMBDA(const Point &p) {
                return amplitude * device::exp(-decay * norm2(c1 - p)) +
                       amplitude * device::exp(-decay * norm2(c2 - p));
            };

        if(model == "poisson") {
            auto fe_model = std::make_shared<PoissonFE<FunctionSpace>>(space);
            fe_model->init_forcing_function(ff);
            fun = fe_model;
        } else {
            auto fe_model = std::make_shared<BratuFE<FunctionSpace>>(space);
            fe_model->init_forcing_function(ff);
            fun = fe_model;
        }

        stats.stop_and_collect("projection");

        stats.start();

        PetscVector x;
        space.create_vector(x);

        Newton<PetscMatrix, PetscVector> newton;
        newton.verbose(false);

        auto linear_solver = std::make_shared<ConjugateGradient<PetscMatrix, PetscVector, HOMEMADE>>();
        auto prec = std::make_shared<InvDiagPreconditioner<PetscMatrix, PetscVector>>();
        linear_solver->set_preconditioner(prec);

        newton.set_linear_solver(linear_solver);

        in.get("newton", newton);

        space.apply_constraints(x);
        newton.solve(*fun, x);

        stats.stop_and_collect("solve+assemblies");

        stats.start();

        rename("x", x);
        space.write(output_path, x);

        stats.stop_and_collect("write");
        stats.describe(std::cout);
    }

    UTOPIA_REGISTER_APP(petsc_bratu);

    static void petsc_dm_app()
    {
        using Mesh = utopia::PetscDM<2>;
        using SizeType = Mesh::SizeType;
        SizeType nx = 10;
        SizeType ny = 10;

        PetscCommunicator world;
        Mesh dm(
            world,
            {nx, ny},
            {0.0, 0.0},
            {1.0, 1.0}
        );

        PetscMatrix mat;
        dm.create_matrix(mat);

        dm.each_element([](const Mesh::Elem &e) {
            // std::cout << e.idx() << std::endl;
        });

        std::stringstream ss;

        dm.each_node([&ss](const Mesh::Node &node) {
            assert(!node.is_ghost());
        });

        dm.each_node_with_ghosts([&ss](const Mesh::Node &node) {
            ss << "(" << node.idx() <<  ", " << node.is_ghost() << ") ";
        });

        int size = world.size();
        int rank = world.rank();

        world.barrier();

        for(int i = 0; i < size; ++i) {
            if(i == rank) {
                std::cout << "--------------------------------------------\n";
                std::cout << ss.str() << std::endl;
                std::cout << "--------------------------------------------\n";
            }

            world.barrier();
        }
    }

    UTOPIA_REGISTER_APP(petsc_dm_app);

    template<class FunctionSpace>
    void compression_only_bc(FunctionSpace &space)
    {
        using Point  = typename FunctionSpace::Point;
        using Scalar = typename FunctionSpace::Scalar;
        ////////////////////////////////////

        //HORIZONTAL
        space.emplace_dirichlet_condition(
            SideSet::left(),
            UTOPIA_LAMBDA(const Point &) -> Scalar {
                return 0.0;
            },
            1
        );

        space.emplace_dirichlet_condition(
            SideSet::left(),
            UTOPIA_LAMBDA(const Point &) -> Scalar {
                return 0.0;
            },
            2
        );

        space.emplace_dirichlet_condition(
            SideSet::right(),
            UTOPIA_LAMBDA(const Point &) -> Scalar {
                return 0.0;
            },
            1
        );

        space.emplace_dirichlet_condition(
            SideSet::right(),
            UTOPIA_LAMBDA(const Point &) -> Scalar {
                return 0.0;
            },
            2
        );

        ////////////////////////////////////

        //VERTICAL

        space.emplace_dirichlet_condition(
            SideSet::top(),
            UTOPIA_LAMBDA(const Point &) -> Scalar {
                return 0.0;
            },
            0
        );

        space.emplace_dirichlet_condition(
            SideSet::top(),
            UTOPIA_LAMBDA(const Point &) -> Scalar {
                return 0.0;
            },
            2
        );

        space.emplace_dirichlet_condition(
            SideSet::bottom(),
            UTOPIA_LAMBDA(const Point &) -> Scalar {
                return 0.0;
            },
            0
        );

        space.emplace_dirichlet_condition(
            SideSet::bottom(),
            UTOPIA_LAMBDA(const Point &) -> Scalar {
                return 0.0;
            },
            2
        );

        ////////////////////////////////////


       //DEPTH
       space.emplace_dirichlet_condition(
           SideSet::front(),
           UTOPIA_LAMBDA(const Point &) -> Scalar {
               return 0.0;
           },
           0
       );

       space.emplace_dirichlet_condition(
           SideSet::front(),
           UTOPIA_LAMBDA(const Point &) -> Scalar {
               return 0.0;
           },
           1
       );

       space.emplace_dirichlet_condition(
           SideSet::back(),
           UTOPIA_LAMBDA(const Point &) -> Scalar {
               return 0.0;
           },
           0
       );

       space.emplace_dirichlet_condition(
           SideSet::back(),
           UTOPIA_LAMBDA(const Point &) -> Scalar {
               return 0.0;
           },
           1
       );

       ////////////////////////////////////


       space.emplace_dirichlet_condition(
           SideSet::left(),
           UTOPIA_LAMBDA(const Point &) -> Scalar {
               return -0.1;
           },
           0
       );

       space.emplace_dirichlet_condition(
           SideSet::right(),
           UTOPIA_LAMBDA(const Point &) -> Scalar {
               return 0.1;
           },
           0
       );

       space.emplace_dirichlet_condition(
           SideSet::bottom(),
           UTOPIA_LAMBDA(const Point &) -> Scalar {
               return -0.1;
           },
           1
       );

       space.emplace_dirichlet_condition(
           SideSet::top(),
           UTOPIA_LAMBDA(const Point &) -> Scalar {
               return 0.1;
           },
           1
       );


       space.emplace_dirichlet_condition(
           SideSet::back(),
           UTOPIA_LAMBDA(const Point &) -> Scalar {
               return -0.1;
           },
           2
       );

       space.emplace_dirichlet_condition(
           SideSet::front(),
           UTOPIA_LAMBDA(const Point &) -> Scalar {
               return 0.1;
           },
           2
       );
    }

    static void petsc_fe_function()
    {
        static const int Dim = 3;
        static const int NVars = Dim;

        using Comm           = utopia::PetscCommunicator;
        using Mesh           = utopia::PetscDM<Dim>;
        using Elem           = utopia::PetscUniformHex8;
        using FunctionSpace  = utopia::FunctionSpace<Mesh, NVars, Elem>;
        using ElemView       = FunctionSpace::ViewDevice::Elem;
        using SizeType       = FunctionSpace::SizeType;
        using Scalar         = FunctionSpace::Scalar;
        using Quadrature     = utopia::Quadrature<Elem, 2>;
        using Dev            = FunctionSpace::Device;
        using FEFunction     = utopia::FEFunction<FunctionSpace>;

        //BEGIN: Host context
        Comm world;

        SizeType scale = (world.size() + 1);
        SizeType nx = scale * 2;
        SizeType ny = scale * 2;
        SizeType nz = scale * 2;

        FunctionSpace space;

        space.build(
            world,
            {nx, ny, nz},
            {0.0, 0.0, 0.0},
            {1.0, 1.0, 1.0}
        );

        FEFunction fun(space);

        Quadrature q;

        auto coeff = fun.coefficient();

        auto f = fun.value(q);
        auto g = fun.gradient(q);

        auto shape      = space.shape(q);
        auto shape_grad = space.shape_grad(q);

        //custom operator can be create with factory functions
        auto lapl       = laplacian(space, q);

        //END: Host context

        {
            //BEGIN: Device context
            auto space_view = space.view_device();
            auto coeff_view = coeff.view_device();
            auto f_view     = f.view_device();
            auto g_view     = g.view_device();

            auto shape_view      = shape.view_device();
            auto shape_grad_view = shape_grad.view_device();


            // Device Kernel (GPU or CPU) (this should be hidden better)
            Dev::parallel_for(space.local_element_range(), UTOPIA_LAMBDA(const SizeType &idx) {
                ElemView e;
                space_view.elem(idx, e);

                auto s_grad = shape_grad_view.make(e);

            });

            //END: Device context
        }

    }

    UTOPIA_REGISTER_APP(petsc_fe_function);

    static void petsc_sample()
    {
        static const int Dim = 2;
        static const int NVars = 1;

        using Comm           = utopia::PetscCommunicator;
        using Mesh           = utopia::PetscDM<Dim>;
        using Elem           = utopia::PetscUniformQuad4;
        using FunctionSpace  = utopia::FunctionSpace<Mesh, NVars, Elem>;
        using ElemView       = FunctionSpace::ViewDevice::Elem;
        using SizeType       = FunctionSpace::SizeType;
        using Scalar         = FunctionSpace::Scalar;
        using Quadrature     = utopia::Quadrature<Elem, 2>;
        using Dev            = FunctionSpace::Device;
        using FEFunction     = utopia::FEFunction<FunctionSpace>;
        using Point          = typename FunctionSpace::Point;

        Comm world;

        SizeType scale = (world.size() + 1);
        SizeType nx = scale * 10;
        SizeType ny = scale * 10;

        FunctionSpace space;

        space.build(
            world,
            {nx, ny},
            {0.0, 0.0},
            {1.0, 1.0}
        );

        PetscVector x;
        space.create_vector(x);

        x.set(0.0);

        auto sampler = utopia::sampler(space, UTOPIA_LAMBDA(const Point &x) {
            auto dist_x = 0.5 - x[0];
            return device::exp(-10.0 * dist_x * dist_x);
        });

        {
            auto space_view   = space.view_device();
            auto x_view       = space.assembly_view_device(x);
            auto sampler_view = sampler.view_device();

            Dev::parallel_for(space.local_element_range(), UTOPIA_LAMBDA(const SizeType &i) {
                ElemView e;
                space_view.elem(i, e);

                StaticVector<Scalar, 4> s;
                sampler_view.assemble(e, s);
                space_view.set_vector(e, s, x_view);
            });

        }

        rename("C", x);
        space.write("sample.vtk", x);
    }

    UTOPIA_REGISTER_APP(petsc_sample);

}

