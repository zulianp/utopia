
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
    static void poisson_problem(FunctionSpace &space, const bool use_direct_solver, const bool debug_matrices = false)
    {
        // using Mesh             = typename FunctionSpace::Mesh;
        using Elem             = typename FunctionSpace::Elem;
        using Dev              = typename FunctionSpace::Device;
        using Vector           = typename FunctionSpace::Vector;
        using Matrix           = typename FunctionSpace::Matrix;
        using Comm             = typename FunctionSpace::Comm;

        static const int Dim    = Elem::Dim;
        // static const int NNodes = Elem::NNodes;
        static const int NFunctions = Elem::NFunctions;

        // using DevFunctionSpace = typename FunctionSpace::ViewDevice;
        using Point            = typename FunctionSpace::Point;
        using Scalar           = typename FunctionSpace::Scalar;
        using SizeType         = typename FunctionSpace::SizeType;
        using Quadrature       = utopia::Quadrature<Elem, 2>;
        using ElementMatrix    = utopia::StaticMatrix<Scalar, NFunctions, NFunctions>;
        // using ElementVector    = utopia::StaticVector<Scalar, NFunctions>;

        Comm &comm = space.comm();

        MPITimeStatistics stats(comm);

        Quadrature quadrature;
        auto &&space_view = space.view_device();

        comm.barrier();
        stats.start();

        Matrix mat, mass_mat;
        space.create_matrix(mat);
        // space.create_matrix(mass_mat);
        //copying is faster than create_matrix
        mass_mat = mat;

        Vector rhs;
        space.create_vector(rhs);

        stats.stop_and_collect("create-matrix");

        for(int c = 0; c < space.n_components(); ++c) {
            space.emplace_dirichlet_condition(
                SideSet::left(),
                UTOPIA_LAMBDA(const Point &) -> Scalar {
                    return 1.0;
                },
                c
            );

            space.emplace_dirichlet_condition(
                SideSet::right(),
                UTOPIA_LAMBDA(const Point &) -> Scalar {
                    return -1.0;
                },
                c
            );
        }

        auto diffusivity = UTOPIA_LAMBDA(const Point &p) -> Scalar {
            Scalar dist = 0.0;
            for(SizeType i = 0; i < p.size(); ++i) {
                Scalar v = p[i] - 0.5;
                dist += v*v;
            }

            if(device::sqrt(dist) > 0.2) {
                return 1.0;
            } else {
                return 1e-4;
            }
        };

        stats.start();


        auto lapl   = laplacian(space, quadrature);
        auto mass_m = mass_matrix(space, quadrature);

        {
            //GPU assembly mock-prototype
            auto mat_view      = space.assembly_view_device(mat);
            auto mass_mat_view = space.assembly_view_device(mass_mat);
            auto rhs_view      = space.assembly_view_device(rhs);

            auto l_view = lapl.view_device();
            auto m_view = mass_m.view_device();

            if(debug_matrices) {
                disp("lapl");
                l_view.describe();

                disp("mass_m");
                m_view.describe();
            }

            Dev::parallel_for(
                space.local_element_range(),
                UTOPIA_LAMBDA(const SizeType &i)
            {
                Elem e;

                //FIXME this is too big for GPU stack memory for hexas
                ElementMatrix el_mat;
                Point c;
                space_view.elem(i, e);

                //Assemble local laplacian
                el_mat.set(0.0);
                l_view.assemble(i, e, el_mat);
                e.centroid(c);
                el_mat *= diffusivity(c);
                space_view.add_matrix(e, el_mat, mat_view);

                //Assemble local mass-matrix and reuse el_mat
                el_mat.set(0.0);
                m_view.assemble(i, e, el_mat);
                space_view.add_matrix(e, el_mat, mass_mat_view);
            });
        }

        stats.stop_and_collect("assemblies");

        stats.start();
        rhs = mass_mat * rhs;

        // rename("m",  mass_mat);
        // write("M.m", mass_mat);

        // rename("a",  mat);
        // write("A.m", mat);

        Scalar vol = sum(mass_mat);
        std::cout << "vol: " << vol << std::endl;

        Scalar zero = sum(mat);
        std::cout << "zero: " << zero << std::endl;

        space.apply_constraints(mat, rhs);

        stats.stop_and_collect("boundary conditions ");

        stats.start();
        Vector x = rhs;
        x.set(0.0);

        if(use_direct_solver) {
            Factorization<Matrix, Vector> solver;
            solver.solve(mat, rhs, x);
        }  else {
            ConjugateGradient<Matrix, Vector, HOMEMADE> cg;
            auto prec = std::make_shared<InvDiagPreconditioner<Matrix, Vector>>();
            cg.set_preconditioner(prec);
            cg.verbose(true);

            const SizeType n_iter = space.n_dofs();

            assert(n_iter > 0);

            cg.max_it(n_iter);
            cg.rtol(1e-8);
            cg.solve(mat, rhs, x);
        }

        stats.stop_and_collect("solve");

        stats.start();
        rename("x", x);
        space.write("X.vtk", x);
        stats.stop_and_collect("write");

        if(comm.rank() == 0) std::cout << "n_dofs: " << space.n_dofs() << std::endl;

        stats.describe(std::cout);
    }

    static void petsc_dm_assemble_2()
    {
        static const int Dim = 2;
        static const int NNodes = 4;

        using Mesh             = utopia::PetscDM<Dim>;
        using Elem             = utopia::PetscUniformQuad4;
        using FunctionSpace    = utopia::FunctionSpace<Mesh, 1, Elem>;
        using SizeType         = Mesh::SizeType;

        PetscCommunicator world;

        SizeType scale = (world.size() + 1);
        SizeType nx = scale * 10;
        SizeType ny = scale * 10;
        SizeType nz = 10;

        Mesh mesh(
            world,
            {nx, ny},
            {0.0, 0.0},
            {1.0, 1.0}
        );

        FunctionSpace space(mesh);
        poisson_problem(space, true);
    }

    UTOPIA_REGISTER_APP(petsc_dm_assemble_2);

    static void petsc_dm_assemble_3()
    {
        static const int Dim = 3;
        static const int NNodes = 8;

        using Mesh             = utopia::PetscDM<Dim>;
        using Elem             = utopia::PetscUniformHex8;
        using FunctionSpace    = utopia::FunctionSpace<Mesh, 1, Elem>;
        using SizeType         = Mesh::SizeType;

        PetscCommunicator world;

        SizeType scale = (world.size() + 1);
        SizeType nx = scale * 15;
        SizeType ny = scale * 15;
        SizeType nz = scale * 15;

        Mesh mesh(
            world,
            {nx, ny, nz},
            {0.0, 0.0, 0.0},
            {1.0, 1.0, 1.0}
        );

        FunctionSpace space(mesh);
        poisson_problem(space, false);
    }

    UTOPIA_REGISTER_APP(petsc_dm_assemble_3);


    static void petsc_dm_mvar_poisson_2()
    {
        static const int Dim = 2;
        static const int NVars = Dim;

        using Mesh             = utopia::PetscDM<Dim>;
        using Elem             = utopia::PetscUniformQuad4;
        using FunctionSpace    = utopia::FunctionSpace<Mesh, NVars, Elem>;
        using SizeType         = Mesh::SizeType;

        PetscCommunicator world;

        SizeType scale = (world.size() + 1);
        SizeType nx = scale * 30;
        SizeType ny = scale * 30;

        FunctionSpace space;
        space.build(
            world,
            {nx, ny},
            {0.0, 0.0},
            {1.0, 1.0}
        );

        poisson_problem(space, true, true);
    }

    UTOPIA_REGISTER_APP(petsc_dm_mvar_poisson_2);

    static void petsc_dm_mvar_poisson_3()
    {
        static const int Dim = 3;
        static const int NVars = Dim;

        using Mesh             = utopia::PetscDM<Dim>;
        using Elem             = utopia::PetscUniformHex8;
        using FunctionSpace    = utopia::FunctionSpace<Mesh, NVars, Elem>;
        using SizeType         = Mesh::SizeType;

        PetscCommunicator world;

        SizeType scale = (world.size() + 1);
        SizeType nx = scale * 20;
        SizeType ny = scale * 20;
        SizeType nz = scale * 20;

        FunctionSpace space;

        space.build(
            world,
            {nx, ny, nz},
            {0.0, 0.0, 0.0},
            {1.0, 1.0, 1.0}
        );

        poisson_problem(space, false);
    }

    UTOPIA_REGISTER_APP(petsc_dm_mvar_poisson_3);

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

    template<class FunctionSpace>
    static void linear_elasticity(FunctionSpace &space, const bool use_direct_solver, const bool debug_matrices = false)
    {
        // using Mesh             = typename FunctionSpace::Mesh;
        using Elem             = typename FunctionSpace::Elem;
        using Dev              = typename FunctionSpace::Device;
        using Vector           = typename FunctionSpace::Vector;
        using Matrix           = typename FunctionSpace::Matrix;
        using Comm             = typename FunctionSpace::Comm;

        static const int Dim    = Elem::Dim;
        // static const int NNodes = Elem::NNodes;
        static const int NFunctions = Elem::NFunctions;

        // using DevFunctionSpace = typename FunctionSpace::ViewDevice;
        using Point            = typename FunctionSpace::Point;
        using Scalar           = typename FunctionSpace::Scalar;
        using SizeType         = typename FunctionSpace::SizeType;
        using Quadrature       = utopia::Quadrature<Elem, 2>;
        using ElementMatrix    = utopia::StaticMatrix<Scalar, NFunctions, NFunctions>;
        // using ElementVector    = utopia::StaticVector<Scalar, NFunctions>;

        Comm &comm = space.comm();

        MPITimeStatistics stats(comm);

        Quadrature quadrature;
        auto &&space_view = space.view_device();

        comm.barrier();
        stats.start();

        Matrix mat, mass_mat;
        space.create_matrix(mat);
        // space.create_matrix(mass_mat);
        //copying is cheaper than create_matrix
        mass_mat = mat;

        Vector rhs;
        space.create_vector(rhs);

        stats.stop_and_collect("create-matrix");

        // compression_only_bc(space);

        for(int c = 0; c < space.n_components(); ++c) {
            space.emplace_dirichlet_condition(
                SideSet::top(),
                UTOPIA_LAMBDA(const Point &) -> Scalar {
                    return 0.1 * Scalar(c==1);
                },
                c
            );

            space.emplace_dirichlet_condition(
                SideSet::bottom(),
                UTOPIA_LAMBDA(const Point &) -> Scalar {
                    return -0.1*Scalar(c==1);
                },
                c
            );
        }

        stats.start();

        LinearElasticity<FunctionSpace, Quadrature> elast(space, quadrature);
        MassMatrix<FunctionSpace, Quadrature> mass_matrix(space, quadrature);

        {
            auto mat_view      = space.assembly_view_device(mat);
            auto mass_mat_view = space.assembly_view_device(mass_mat);
            auto rhs_view      = space.assembly_view_device(rhs);

            auto elast_view = elast.view_device();
            auto m_view     = mass_matrix.view_device();

            if(debug_matrices) {
                disp("elast");
                elast_view.describe();

                disp("mass_matrix");
                m_view.describe();
            }

            Dev::parallel_for(
                space.local_element_range(),
                UTOPIA_LAMBDA(const SizeType &i)
            {
                Elem e;

                //FIXME this is too big for GPU stack memory for hexas
                ElementMatrix el_mat;
                space_view.elem(i, e);

                //Assemble local elast
                el_mat.set(0.0);
                elast_view.assemble(i, e, el_mat);
                space_view.add_matrix(e, el_mat, mat_view);

                //Assemble local mass-matrix and reuse el_mat
                el_mat.set(0.0);
                m_view.assemble(i, e, el_mat);
                space_view.add_matrix(e, el_mat, mass_mat_view);
            });
        }

        stats.stop_and_collect("assemblies");

        stats.start();
        rhs = mass_mat * rhs;

        Scalar vol = sum(mass_mat);
        std::cout << "vol: " << vol << std::endl;

        Scalar zero = sum(mat);
        std::cout << "zero: " << zero << std::endl;

        space.apply_constraints(mat, rhs);

        stats.stop_and_collect("boundary conditions ");

        stats.start();
        Vector x = rhs;
        x.set(0.0);

        disp("Solving...");

        if(use_direct_solver) {
            Factorization<Matrix, Vector> solver;
            solver.solve(mat, rhs, x);
        }  else {
            ConjugateGradient<Matrix, Vector, HOMEMADE> cg;
            auto prec = std::make_shared<InvDiagPreconditioner<Matrix, Vector>>();
            cg.set_preconditioner(prec);
            cg.verbose(true);

            const SizeType n_iter = space.n_dofs();

            assert(n_iter > 0);

            cg.max_it(n_iter);
            cg.rtol(1e-6);
            cg.atol(1e-6);
            cg.solve(mat, rhs, x);
        }

        stats.stop_and_collect("solve");

        stats.start();
        compute_strain_energy_splitting(space, x);
        stats.stop_and_collect("splitting");


        stats.start();
        rename("x", x);
        space.write("X.vtk", x);
        stats.stop_and_collect("write");

        // rename("rhs", rhs);
        // space.write("R.vtk", rhs);

        if(comm.rank() == 0) std::cout << "n_dofs: " << space.n_dofs() << std::endl;

        stats.describe(std::cout);
    }


    static void petsc_elasticity_3()
    {
        static const int Dim = 3;
        static const int NVars = Dim;

        using Mesh             = utopia::PetscDM<Dim>;
        using Elem             = utopia::PetscUniformHex8;
        using FunctionSpace    = utopia::FunctionSpace<Mesh, NVars, Elem>;
        using SizeType         = Mesh::SizeType;

        PetscCommunicator world;

        SizeType scale = (world.size() + 1);
        SizeType nx = scale * 20;
        SizeType ny = scale * 20;
        SizeType nz = scale * 20;

        FunctionSpace space;

        space.build(
            world,
            {nx, ny, nz},
            {0.0, 0.0, 0.0},
            {1.0, 1.0, 1.0}
        );

        linear_elasticity(space, false, false);

    }

    UTOPIA_REGISTER_APP(petsc_elasticity_3);

    static void petsc_strain()
    {
        static const int Dim = 3;
        static const int NVars = Dim;

        using Mesh           = utopia::PetscDM<Dim>;
        using Elem           = utopia::PetscUniformHex8;
        using FunctionSpace  = utopia::FunctionSpace<Mesh, NVars, Elem>;
        using ElemView       = FunctionSpace::ViewDevice::Elem;
        using SizeType       = Mesh::SizeType;
        using Scalar         = Mesh::Scalar;
        using Quadrature     = utopia::Quadrature<Elem, 2>;
        using Dev            = FunctionSpace::Device;
        using VectorD        = utopia::StaticVector<Scalar, Dim>;
        // using MatrixDxD        = utopia::StaticMatrix<Scalar, Dim, Dim>;

        PetscCommunicator world;

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

        PetscVector u;
        space.create_vector(u);
        u.set(0.1);

        compute_strain_energy_splitting(space, u);
    }

    UTOPIA_REGISTER_APP(petsc_strain);


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



    static void petsc_phase_field(Input &in)
    {
        static const int Dim = 3;
        static const int NVars = Dim + 1;

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
        using Point          = typename FunctionSpace::Point;

        Comm world;

        SizeType scale = (world.size() + 1);
        SizeType nx = scale * 5;
        SizeType ny = scale * 5;
        SizeType nz = scale * 5;

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

        PhaseFieldForBrittleFractures<FunctionSpace> pp(space);

        PetscMatrix H;
        PetscVector x, g;
        Scalar f;

        space.create_vector(x);

        x.set(0.0);


        auto C = space.subspace(0);
        auto sampler = utopia::sampler(C, UTOPIA_LAMBDA(const Point &x) {
            auto dist_x = 0.5 - x[0];
            return device::exp(-5.0 * dist_x * dist_x);
        });

        {
            auto C_view       = C.view_device();
            auto sampler_view = sampler.view_device();
            auto x_view       = space.assembly_view_device(x);

            Dev::parallel_for(space.local_element_range(), UTOPIA_LAMBDA(const SizeType &i) {
                utopia::FunctionSpace<Mesh, 1, Elem>::ViewDevice::Elem e;
                C_view.elem(i, e);

                StaticVector<Scalar, 8> s;
                sampler_view.assemble(e, s);
                C_view.set_vector(e, s, x_view);
            });

        }

        std::string output_path = "phase_field.vtu";

        in.get("output-path", output_path);

        rename("X", x);


        C.write(output_path, x);

        pp.assemble(x, H, g, f);
    }

    UTOPIA_REGISTER_APP(petsc_phase_field);
}

