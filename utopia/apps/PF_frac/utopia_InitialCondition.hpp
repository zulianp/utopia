#ifndef UTOPIA_INITIAL_CONDITION_PF_HPP
#define UTOPIA_INITIAL_CONDITION_PF_HPP

#include "utopia_Base.hpp"

// FIXME: this file causes nvcc to fail
#ifndef KOKKOS_ENABLE_CUDA

#include "utopia_Base.hpp"
#include "utopia_RangeDevice.hpp"

// include edsl components
#include "utopia_AppRunner.hpp"
#include "utopia_AssemblyView.hpp"
#include "utopia_ConjugateGradient.hpp"
#include "utopia_Core.hpp"
#include "utopia_DeviceView.hpp"
#include "utopia_FEFunction.hpp"
#include "utopia_FracNetGenerator2D.hpp"
#include "utopia_FracNetGenerator3D.hpp"
#include "utopia_GradInterpolate.hpp"
#include "utopia_InitialConditionBase.hpp"
#include "utopia_IsotropicPhaseField.hpp"
#include "utopia_LinearElasticityView.hpp"
#include "utopia_MPITimeStatistics.hpp"
#include "utopia_MassMatrixView.hpp"
#include "utopia_PrincipalStrainsView.hpp"
#include "utopia_SampleView.hpp"
#include "utopia_TrivialPreconditioners.hpp"

#include "utopia_petsc.hpp"
#include "utopia_petsc_DM.hpp"
#include "utopia_petsc_DMDA_FunctionSpace.hpp"
#include "utopia_petsc_DirichletBoundaryConditions.hpp"
#include "utopia_petsc_Matrix.hpp"
#include "utopia_petsc_impl.hpp"

#include <chrono>
#include <cmath>
#include <random>
#include <algorithm>

namespace utopia {

    template <class FunctionSpace>
    class InitialCondidtionPFSneddon : public InitialCondition<FunctionSpace> {
    public:
        // using Comm           = typename FunctionSpace::Comm;
        using Mesh = typename FunctionSpace::Mesh;
        using Elem = typename FunctionSpace::Shape;
        using ElemView = typename FunctionSpace::ViewDevice::Elem;
        using SizeType = typename FunctionSpace::SizeType;
        using Scalar = typename FunctionSpace::Scalar;
        using Dev = typename FunctionSpace::Device;
        using Point = typename FunctionSpace::Point;
        using ElemViewScalar = typename utopia::FunctionSpace<Mesh, 1, Elem>::ViewDevice::Elem;
        static const int NNodes = Elem::NNodes;
        static const int Dim = FunctionSpace::Dim;

        InitialCondidtionPFSneddon(FunctionSpace &space, const SizeType &PF_component)
            : InitialCondition<FunctionSpace>(space), PF_component_(PF_component), l_0_(1.0) {}

        void read(Input &in) override { in.get("l_0", l_0_); }

        void init(PetscVector &x) override {
            auto C = this->space_.subspace(PF_component_);

            auto sampler = utopia::sampler(
                C, UTOPIA_LAMBDA(const Point &x)->Scalar {
                    const Scalar thickness = 3.0 * this->space_.mesh().min_spacing();

                    Scalar r_squared = 0.0;

                    if (Dim == 2) {
                        r_squared = x[0] * x[0];
                    } else {
                        r_squared = (x[0] * x[0]) + (x[2] * x[2]);
                    }

                    Scalar f = 0.0;
                    if ((r_squared <= l_0_ * l_0_) && ((device::abs(2.0 * x[1]) <= thickness))) {
                        f = 1.0;
                    } else {
                        f = 0.0;
                    }
                    return f;
                });

            {
                auto C_view = C.view_device();
                auto sampler_view = sampler.view_device();
                auto x_view = this->space_.assembly_view_device(x);

                Dev::parallel_for(
                    this->space_.element_range(), UTOPIA_LAMBDA(const SizeType &i) {
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
        Scalar l_0_;
    };

    template <class FunctionSpace>
    class InitialCondidtionPFTension : public InitialCondition<FunctionSpace> {
    public:
        // using Comm           = typename FunctionSpace::Comm;
        using Mesh = typename FunctionSpace::Mesh;
        using Elem = typename FunctionSpace::Shape;
        using ElemView = typename FunctionSpace::ViewDevice::Elem;
        using SizeType = typename FunctionSpace::SizeType;
        using Scalar = typename FunctionSpace::Scalar;
        using Dev = typename FunctionSpace::Device;
        using Point = typename FunctionSpace::Point;
        using ElemViewScalar = typename utopia::FunctionSpace<Mesh, 1, Elem>::ViewDevice::Elem;
        static const int NNodes = Elem::NNodes;

        InitialCondidtionPFTension(FunctionSpace &space, const SizeType &PF_component)
            : InitialCondition<FunctionSpace>(space), PF_component_(PF_component) {}

        void read(Input & in) override {
            in.get("initial_crack_width", initial_crack_width_);
        }

        void init(PetscVector &x) override {
            // un-hard-code
            auto C = this->space_.subspace(PF_component_);

            //model dimensions
            auto xyz_min = this->space_.mesh().box_min();
            auto xyz_max = this->space_.mesh().box_max();

            double dx = this->space_.mesh().min_spacing();
            double y_mid = (xyz_max[1]-xyz_min[1])/2.0 - dx/20.0;


            if (mpi_world_rank() == 0){
                std::cout << "Minimum Mesh Spacing: " <<  dx << std::endl;
            }

            auto sampler = utopia::sampler(
                C, UTOPIA_LAMBDA(const Point &x)->Scalar {
                    Scalar f = 0.0;
                    if (x[1] > (y_mid - initial_crack_width_*dx) &&
                        x[1] < (y_mid + initial_crack_width_*dx) && x[0] < y_mid) {
                        // if (x[0] <= 0.5 && x[1] <= 0.5) {
                        // f = 1.0;
                        f = 1.0;
                        // f = 0.5;
                    } else {
                        f = 0.0;
                        // f = 0.5;
                    }
                    return f;
                });

            {
                auto C_view = C.view_device();
                auto sampler_view = sampler.view_device();
                auto x_view = this->space_.assembly_view_device(x);

                Dev::parallel_for(
                    this->space_.element_range(), UTOPIA_LAMBDA(const SizeType &i) {
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
        double initial_crack_width_;
    };

    template <class FunctionSpace>
    class InitialCondidtionPFTbar : public InitialCondition<FunctionSpace> {
    public:
        // using Comm           = typename FunctionSpace::Comm;
        using Mesh = typename FunctionSpace::Mesh;
        using Elem = typename FunctionSpace::Shape;
        using ElemView = typename FunctionSpace::ViewDevice::Elem;
        using SizeType = typename FunctionSpace::SizeType;
        using Scalar = typename FunctionSpace::Scalar;
        using Dev = typename FunctionSpace::Device;
        using Point = typename FunctionSpace::Point;
        using ElemViewScalar = typename utopia::FunctionSpace<Mesh, 1, Elem>::ViewDevice::Elem;
        static const int NNodes = Elem::NNodes;

        InitialCondidtionPFTbar(FunctionSpace &space, const SizeType &PF_component)
            : InitialCondition<FunctionSpace>(space), PF_component_(PF_component) {}

        void init(PetscVector &x) override {
            // un-hard-code
            auto C = this->space_.subspace(PF_component_);

            auto sampler = utopia::sampler(
                C, UTOPIA_LAMBDA(const Point &x)->Scalar {
                    Scalar f = 0.0;
                    if (x[0] > (0.5 - this->space_.mesh().min_spacing()) &&
                        x[0] < (0.5 + this->space_.mesh().min_spacing()) && x[1] > 0.3 && x[1] < 0.5) {
                        f = 1.0;
                    } else if ((x[0] > 0.4) && (x[0] < 0.6) && (x[1] > 0.7 - this->space_.mesh().min_spacing()) &&
                               (x[1] < 0.7 + this->space_.mesh().min_spacing())) {
                        f = 1.0;
                    } else {
                        f = 0.0;
                    }

                    return f;
                });

            {
                auto C_view = C.view_device();
                auto sampler_view = sampler.view_device();
                auto x_view = this->space_.assembly_view_device(x);

                Dev::parallel_for(
                    this->space_.element_range(), UTOPIA_LAMBDA(const SizeType &i) {
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

    template <class FunctionSpace>
    class InitialCondidtionPFParallelFrac3D : public InitialCondition<FunctionSpace> {
    public:
        // using Comm           = typename FunctionSpace::Comm;
        using Mesh = typename FunctionSpace::Mesh;
        using Elem = typename FunctionSpace::Shape;
        using ElemView = typename FunctionSpace::ViewDevice::Elem;
        using SizeType = typename FunctionSpace::SizeType;
        using Scalar = typename FunctionSpace::Scalar;
        using Dev = typename FunctionSpace::Device;
        using Point = typename FunctionSpace::Point;
        using ElemViewScalar = typename utopia::FunctionSpace<Mesh, 1, Elem>::ViewDevice::Elem;
        static const int NNodes = Elem::NNodes;

        InitialCondidtionPFParallelFrac3D(FunctionSpace &space, const SizeType &PF_component)
            : InitialCondition<FunctionSpace>(space), PF_component_(PF_component) {}

        void init(PetscVector &x) override {
            // un-hard-code
            auto C = this->space_.subspace(PF_component_);

            auto sampler = utopia::sampler(
                C, UTOPIA_LAMBDA(const Point &x)->Scalar {
                    Scalar f = 0.0;
                    Scalar h = this->space_.mesh().min_spacing();

                    if (x[0] > 0.55 && x[0] < 0.70 && x[1] > (0.40 - h) && x[1] < (0.40 + h) && x[2] < (0.6 + h) &&
                        x[2] > (0.6 - h)) {
                        f = 1.0;
                    } else if (x[0] > (0.26 - h) && x[0] < (0.26 + h) && x[1] > 0.38 && x[1] < 0.55 &&
                               x[2] < (0.40 + h) && x[2] > (0.40 - h)) {
                        f = 1.0;
                    } else {
                        f = 0.0;
                    }
                    return f;
                });

            {
                auto C_view = C.view_device();
                auto sampler_view = sampler.view_device();
                auto x_view = this->space_.assembly_view_device(x);

                Dev::parallel_for(
                    this->space_.element_range(), UTOPIA_LAMBDA(const SizeType &i) {
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

    template <class FunctionSpace>
    class Mixed : public InitialCondition<FunctionSpace> {
    public:
        // using Comm           = typename FunctionSpace::Comm;
        using Mesh = typename FunctionSpace::Mesh;
        using Elem = typename FunctionSpace::Shape;
        using ElemView = typename FunctionSpace::ViewDevice::Elem;
        using SizeType = typename FunctionSpace::SizeType;
        using Scalar = typename FunctionSpace::Scalar;
        using Dev = typename FunctionSpace::Device;
        using Point = typename FunctionSpace::Point;
        using ElemViewScalar = typename utopia::FunctionSpace<Mesh, 1, Elem>::ViewDevice::Elem;
        static const int NNodes = Elem::NNodes;

        Mixed(FunctionSpace &space, const SizeType &PF_component)
            : InitialCondition<FunctionSpace>(space), PF_component_(PF_component), pressure0_(1.0) {}

        void init(PetscVector &x) override {
            using CoeffVector = utopia::StaticVector<Scalar, NNodes>;
            // un-hard-code
            auto C = this->space_.subspace(PF_component_);

            auto width = 3.0 * this->space_.mesh().min_spacing();
            // auto width = 0.1;

            if (mpi_world_rank() == 0) {
                utopia::out() << "width: " << width << "  \n";
            }

            std::vector<Rectangle<Scalar>> rectangles;

            Point2D<Scalar> A{};
            A.x = 38.1;
            A.y = 82.55;
            rectangles.push_back(Rectangle<Scalar>(A, -12.7, width, -45));

            Point2D<Scalar> B{};
            B.x = 38.1;
            B.y = 69.85;
            rectangles.push_back(Rectangle<Scalar>(B, 12.7, width, -45));

            auto sampler = utopia::sampler(C, [&rectangles](const Point &x) -> Scalar {
                for (std::size_t r = 0; r < rectangles.size(); r++) {
                    if (rectangles[r].belongs_to_rectangle(x[0], x[1])) return 1.0;
                }
                return 0.0;
            });

            {
                auto C_view = C.view_device();
                auto sampler_view = sampler.view_device();
                auto x_view = this->space_.assembly_view_device(x);

                Dev::parallel_for(
                    this->space_.element_range(), UTOPIA_LAMBDA(const SizeType &i) {
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

    template <class FunctionSpace>
    class AsphaltTension : public InitialCondition<FunctionSpace> {
    public:
        // using Comm           = typename FunctionSpace::Comm;
        using Mesh = typename FunctionSpace::Mesh;
        using Elem = typename FunctionSpace::Shape;
        using ElemView = typename FunctionSpace::ViewDevice::Elem;
        using SizeType = typename FunctionSpace::SizeType;
        using Scalar = typename FunctionSpace::Scalar;
        using Dev = typename FunctionSpace::Device;
        using Point = typename FunctionSpace::Point;
        using ElemViewScalar = typename utopia::FunctionSpace<Mesh, 1, Elem>::ViewDevice::Elem;
        static const int NNodes = Elem::NNodes;

        AsphaltTension(FunctionSpace &space, const SizeType &PF_component)
            : InitialCondition<FunctionSpace>(space), PF_component_(PF_component), pressure0_(1.0) {}

        void init(PetscVector &x) override {
            using CoeffVector = utopia::StaticVector<Scalar, NNodes>;
            // un-hard-code
            auto C = this->space_.subspace(PF_component_);

            auto width = 3.0 * this->space_.mesh().min_spacing();
            // auto width = 0.2;

            if (mpi_world_rank() == 0) {
                utopia::out() << "width: " << width << "  \n";
            }

            std::vector<Rectangle<Scalar>> rectangles;

            Point2D<Scalar> A{};
            A.x = 12.00;
            A.y = 21.50;
            rectangles.push_back(Rectangle<Scalar>(A, -3.5355, width, 45));

            Point2D<Scalar> B{};
            B.x = 15;
            B.y = 20.00;
            rectangles.push_back(Rectangle<Scalar>(B, 5.000, width, 0.0));

            auto sampler = utopia::sampler(C, [&rectangles](const Point &x) -> Scalar {
                for (std::size_t r = 0; r < rectangles.size(); r++) {
                    if (rectangles[r].belongs_to_rectangle(x[0], x[1])) return 1.0;
                }
                return 0.0;
            });

            {
                auto C_view = C.view_device();
                auto sampler_view = sampler.view_device();
                auto x_view = this->space_.assembly_view_device(x);

                Dev::parallel_for(
                    this->space_.element_range(), UTOPIA_LAMBDA(const SizeType &i) {
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

    template <class FunctionSpace>
    class SingleFault: public InitialCondition<FunctionSpace> {
    public:
        // using Comm           = typename FunctionSpace::Comm;
        using Mesh = typename FunctionSpace::Mesh;
        using Elem = typename FunctionSpace::Shape;
        using ElemView = typename FunctionSpace::ViewDevice::Elem;
        using SizeType = typename FunctionSpace::SizeType;
        using Scalar = typename FunctionSpace::Scalar;
        using Dev = typename FunctionSpace::Device;
        using Point = typename FunctionSpace::Point;
        using ElemViewScalar = typename utopia::FunctionSpace<Mesh, 1, Elem>::ViewDevice::Elem;
        static const int NNodes = Elem::NNodes;

        SingleFault(FunctionSpace &space, const SizeType &PF_component)
            : InitialCondition<FunctionSpace>(space), PF_component_(PF_component){}

        void read(Input & in) override {
            in.get("init_x", init_x);
            in.get("init_y", init_y);
            in.get("fault_thickness", fault_thickness);
            in.get("fault_length", fault_length );
            in.get("number_fractures", number_fractures);
            in.get("fracture_length", fracture_length );
        }

        void init(PetscVector &x) override {
            using CoeffVector = utopia::StaticVector<Scalar, NNodes>;
            // un-hard-code
            auto C = this->space_.subspace(PF_component_);

            auto width = 2.0 * this->space_.mesh().min_spacing();
            // auto width = 0.2;

            if (mpi_world_rank() == 0) {
                utopia::out() << "width: " << width << "  \n";
            }



            double fracture_orientation = 0.0;


            std::vector<double> x_loc;
            std::vector<double> y_loc;
            std::vector<double> f_len;


            for (int i = 0; i < number_fractures; i++){
                x_loc.push_back(init_x + static_cast<double>(i)/static_cast<double>(number_fractures)*fault_length );
                y_loc.push_back(init_y + static_cast<double>(i)/static_cast<double>(number_fractures)*fault_thickness);
                f_len.push_back(fracture_length/2.0 + std::pow( (static_cast<double>(i)/static_cast<double>(number_fractures)),3.0)/(0.5*fracture_length) );
            }

            auto rng = std::default_random_engine {};
            auto rng2 = std::default_random_engine {};
            auto rng3 = std::default_random_engine {};

            rng.seed(10);
            rng2.seed(100);
            rng3.seed(1000);


            std::shuffle(std::begin(x_loc), std::end(x_loc), rng);
            std::shuffle(std::begin(y_loc), std::end(y_loc), rng2);
            std::shuffle(std::begin(f_len), std::end(f_len), rng3);

            std::vector<Rectangle<Scalar>> rectangles;

            for (int i = 0; i < number_fractures; i++){
                Point2D<Scalar> XY{};
                XY.x = x_loc[i];
                XY.y = y_loc[i];
                rectangles.push_back(Rectangle<Scalar>(XY, f_len[i], width, fracture_orientation));

            }

//            Point2D<Scalar> XY{};
//            XY.x = 0.7;
//            XY.y = 5.5;
//            rectangles.push_back(Rectangle<Scalar>(XY, 1.8, width, fracture_orientation));

            auto sampler = utopia::sampler(C, [&rectangles](const Point &x) -> Scalar {
                for (std::size_t r = 0; r < rectangles.size(); r++) {
                    if (rectangles[r].belongs_to_rectangle(x[0], x[1])) return 1.0;
                }
                return 0.0;
            });

            {
                auto C_view = C.view_device();
                auto sampler_view = sampler.view_device();
                auto x_view = this->space_.assembly_view_device(x);

                Dev::parallel_for(
                    this->space_.element_range(), UTOPIA_LAMBDA(const SizeType &i) {
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
        Scalar init_x, init_y, fault_thickness, fault_length, fracture_length, number_fractures;
    };


    template <class FunctionSpace>
    class SlantedCrack2D : public InitialCondition<FunctionSpace> {
    public:
        // using Comm           = typename FunctionSpace::Comm;
        using Mesh = typename FunctionSpace::Mesh;
        using Elem = typename FunctionSpace::Shape;
        using ElemView = typename FunctionSpace::ViewDevice::Elem;
        using SizeType = typename FunctionSpace::SizeType;
        using Scalar = typename FunctionSpace::Scalar;
        using Dev = typename FunctionSpace::Device;
        using Point = typename FunctionSpace::Point;
        using ElemViewScalar = typename utopia::FunctionSpace<Mesh, 1, Elem>::ViewDevice::Elem;
        static const int NNodes = Elem::NNodes;

        SlantedCrack2D(FunctionSpace &space, const SizeType &PF_component)
            : InitialCondition<FunctionSpace>(space), PF_component_(PF_component){}

        void init(PetscVector &x) override {
            using CoeffVector = utopia::StaticVector<Scalar, NNodes>;
            // un-hard-code
            auto C = this->space_.subspace(PF_component_);

            auto width = 3.0 * this->space_.mesh().min_spacing();

            auto xyz_min = this->space_.mesh().box_min();

            // auto width = 0.2;

            if (mpi_world_rank() == 0) {
                utopia::out() << "width: " << width << "  \n";
            }

            std::vector<Rectangle<Scalar>> rectangles;

            Point2D<Scalar> A{};
            A.x = 12.00;
            A.y = 21.50;
            rectangles.push_back(Rectangle<Scalar>(A, -3.5355, width, 45));

            Point2D<Scalar> B{};
            B.x = 15;
            B.y = 20.00;
            //rectangles.push_back(Rectangle<Scalar>(B, 5.000, width, 0.0));

            auto sampler = utopia::sampler(C, [&rectangles](const Point &x) -> Scalar {
                for (std::size_t r = 0; r < rectangles.size(); r++) {
                    if (rectangles[r].belongs_to_rectangle(x[0], x[1])) return 1.0;
                }
                return 0.0;
            });

            {
                auto C_view = C.view_device();
                auto sampler_view = sampler.view_device();
                auto x_view = this->space_.assembly_view_device(x);

                Dev::parallel_for(
                    this->space_.element_range(), UTOPIA_LAMBDA(const SizeType &i) {
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
    };




    template <class FunctionSpace>
    class RandomlyDistributed : public InitialCondition<FunctionSpace> {
    public:
        // using Comm           = typename FunctionSpace::Comm;
        using Mesh = typename FunctionSpace::Mesh;
        using Elem = typename FunctionSpace::Shape;
        using ElemView = typename FunctionSpace::ViewDevice::Elem;
        using SizeType = typename FunctionSpace::SizeType;
        using Scalar = typename FunctionSpace::Scalar;
        using Dev = typename FunctionSpace::Device;
        using Point = typename FunctionSpace::Point;
        using ElemViewScalar = typename utopia::FunctionSpace<Mesh, 1, Elem>::ViewDevice::Elem;
        static const int NNodes = Elem::NNodes;

        RandomlyDistributed(FunctionSpace &space, const SizeType &PF_component)
            : InitialCondition<FunctionSpace>(space), PF_component_(PF_component){}

        void init(PetscVector &x) override {
            using CoeffVector = utopia::StaticVector<Scalar, NNodes>;
            // un-hard-code
            auto C = this->space_.subspace(PF_component_);
            int  total_nodes = this->space_.mesh().n_nodes();
            auto min_elem = this->space_.mesh().min_spacing();

            auto xyz_min = this->space_.mesh().box_min();
            auto xyz_max = this->space_.mesh().box_max();

            auto height = xyz_max(1) - xyz_min(1);
            auto width  = xyz_max(0) - xyz_min(0);

            double width_x = width / 6.0;
            double height_min = 4.0*height/10.0 + xyz_min(1) + 3.0 ;
            double height_max = 5.0*height/10.0 + xyz_min(1) - 3.0;
            double damage_threshold = 1.0;

            std::default_random_engine generator;
            std::poisson_distribution<int> distribution(4);

            if (mpi_world_rank() == 0) {
                utopia::out() << "rand: " << double(rand())/RAND_MAX << "  \n";
            }

            std::vector<Rectangle<Scalar>> rectangles;

            int seed = mpi_world_rank()*total_nodes;
            generator.seed(seed);

            auto sampler = utopia::sampler(C, [&generator, &distribution, xyz_min, xyz_max, width_x, height_min, height_max, total_nodes, damage_threshold](const Point &x) -> Scalar {
                if ( x(1) < height_max && x(1) > height_min &&
                     x(0) > xyz_min(0) + width_x && x(0) < xyz_max(0) - width_x ) {
                    double random_damage = double(distribution(generator))/10.0;
                    if (random_damage >= damage_threshold) return 1.0;
                    else return 0.0;
                } else return 0.0;
                });

            {
                auto C_view = C.view_device();
                auto sampler_view = sampler.view_device();
                auto x_view = this->space_.assembly_view_device(x);

                Dev::parallel_for(
                    this->space_.element_range(), UTOPIA_LAMBDA(const SizeType &i) {
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
    };


    template <class FunctionSpace>
    class UniformSpacingOnLine: public InitialCondition<FunctionSpace> {
    public:
        // using Comm           = typename FunctionSpace::Comm;
        using Mesh = typename FunctionSpace::Mesh;
        using Elem = typename FunctionSpace::Shape;
        using ElemView = typename FunctionSpace::ViewDevice::Elem;
        using SizeType = typename FunctionSpace::SizeType;
        using Scalar = typename FunctionSpace::Scalar;
        using Dev = typename FunctionSpace::Device;
        using Point = typename FunctionSpace::Point;
        using ElemViewScalar = typename utopia::FunctionSpace<Mesh, 1, Elem>::ViewDevice::Elem;
        static const int NNodes = Elem::NNodes;

        UniformSpacingOnLine(FunctionSpace &space, const SizeType &PF_component)
            : InitialCondition<FunctionSpace>(space), PF_component_(PF_component){}

        void read(Input & in) override {
            in.get("IC_number", number_);
            in.get("IC_ycoord", ycoord_);
        }

        void init(PetscVector &x) override {
            using CoeffVector = utopia::StaticVector<Scalar, NNodes>;
            // un-hard-code
            auto C = this->space_.subspace(PF_component_);
            int  total_nodes = this->space_.mesh().n_nodes();
            auto min_elem_x = this->space_.mesh().min_spacing_x();
            auto min_elem_y = this->space_.mesh().min_spacing_y();

            auto xyz_min = this->space_.mesh().box_min();
            auto xyz_max = this->space_.mesh().box_max();

            double height = xyz_max(1) - xyz_min(1);
            double width  = xyz_max(0) - xyz_min(0);
            double width_x = width / 9.0;

            double ycoord = ycoord_;
            double number = number_;

            if (mpi_world_rank() == 0){
                std::cout << "min_elem_x" << "   min_elem_y   min_spacing" << std::endl;
                std::cout << min_elem_x << "   " << min_elem_y << "    " <<  this->space_.mesh().min_spacing() << std::endl;
            }

            auto sampler = utopia::sampler(C, [width, ycoord, min_elem_x, min_elem_y, number, xyz_min, xyz_max, width_x ](const Point &x) -> Scalar {
                if ( x(1) <= ycoord + min_elem_y/2.0 && x(1) >= ycoord - min_elem_y/2.0 &&
                     x(0) > xyz_min(0) + width_x && x(0) < xyz_max(0) - width_x &&
                     fmod( x(0) , width/number) <= 0.5*width/number + 0.501*min_elem_x &&
                     fmod( x(0) , width/number) >= 0.5*width/number - 0.501*min_elem_x  ){
                    //std::cout << x(0) << "  " << x(1) << std::endl;
                    return 1.0;
                }
                else return 0.0;

                });

            {
                auto C_view = C.view_device();
                auto sampler_view = sampler.view_device();
                auto x_view = this->space_.assembly_view_device(x);

                Dev::parallel_for(
                    this->space_.element_range(), UTOPIA_LAMBDA(const SizeType &i) {
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
        double number_;
        double ycoord_;
    };


    template <class FunctionSpace>
    class DamagedSedimentaryLayers: public InitialCondition<FunctionSpace> {
    public:
        // using Comm           = typename FunctionSpace::Comm;
        using Mesh = typename FunctionSpace::Mesh;
        using Elem = typename FunctionSpace::Shape;
        using ElemView = typename FunctionSpace::ViewDevice::Elem;
        using SizeType = typename FunctionSpace::SizeType;
        using Scalar = typename FunctionSpace::Scalar;
        using Dev = typename FunctionSpace::Device;
        using Point = typename FunctionSpace::Point;
        using ElemViewScalar = typename utopia::FunctionSpace<Mesh, 1, Elem>::ViewDevice::Elem;
        static const int NNodes = Elem::NNodes;

        DamagedSedimentaryLayers(FunctionSpace &space, const SizeType &PF_component)
            : InitialCondition<FunctionSpace>(space), PF_component_(PF_component){}

        void read(Input & in) override {
            in.get("top_layer_height", ymax_);
            in.get("bottom_layer_height", ymin_);
            in.get("random_damage", random_damage_);
            in.get("uniform_damage", uniform_damage_);
            in.get("central_damage", central_damage_);
            in.get("IC_number", number_);

        }

        void init(PetscVector &x) override {
            using CoeffVector = utopia::StaticVector<Scalar, NNodes>;
            // un-hard-code

            ///RANDOM DAMAGE /// ============================
            auto C = this->space_.subspace(PF_component_);
            int  total_nodes = this->space_.mesh().n_nodes();

            auto xyz_min = this->space_.mesh().box_min();
            auto xyz_max = this->space_.mesh().box_max();
            auto width  = xyz_max(0) - xyz_min(0);

            double width_x = width / 10.0;
            double height_min = ymin_;
            double height_max = ymax_;
            double height = ymax_ - ymin_;

            double damage_threshold = 1.;

            std::default_random_engine generator;
            std::poisson_distribution<int> distribution(3);

            int seed = mpi_world_rank()*total_nodes;
            generator.seed(seed);
            bool random_damage = random_damage_;

            /// Uniform Damage /// ========================
            auto min_elem_x = this->space_.mesh().min_spacing_x();
            auto min_elem_y = this->space_.mesh().min_spacing_y();

            double number = number_;
            bool uniform_damage = uniform_damage_;

            bool central_damage = central_damage_;
            std::vector<Rectangle<Scalar>> rectangles;

            Point2D<Scalar> A{};
            A.x = xyz_min(0) + width/2.0;
            A.y = height_max - height/3.0;
            auto crack_width = 1.0 * this->space_.mesh().min_spacing();
            rectangles.push_back(Rectangle<Scalar>(A, height/3.0, crack_width, -90.));


            auto sampler = utopia::sampler(C, [&generator, &distribution, &rectangles, central_damage, uniform_damage, random_damage, width, height_max, height_min,
                                             min_elem_x, min_elem_y, number, xyz_min, xyz_max, width_x, height,damage_threshold](const Point &x) -> Scalar {

                if (central_damage){
                    for (std::size_t r = 0; r < rectangles.size(); r++) {
                        if (rectangles[r].belongs_to_rectangle(x[0], x[1])) return 1.0;
                    }
                    return 0.0;
                }

                if (uniform_damage){
                    if ( x(1) <= height_max - height/3.0  && x(1) >= height_max - 2.0*height/3.0 &&
                         x(0) > xyz_min(0) + width && x(0) < xyz_max(0) - width_x &&
                         fmod( x(0) , width/number) <= 0.5*width/number + 0.501*min_elem_x &&
                         fmod( x(0) , width/number) >= 0.5*width/number - 0.501*min_elem_x  ){
                        //std::cout << x(0) << "  " << x(1) << std::endl;
                        return 1.0;
                    }
                    else return 0.0;
                }
                if (random_damage){
                    if ( x(1) < height_max && x(1) > height_min &&
                         x(0) > xyz_min(0) + width_x && x(0) < xyz_max(0) - width_x ) {
                        double rand_damage = double(distribution(generator))/10.0;
                        if (rand_damage >= damage_threshold) return 1.0;
                        else return 0.0;
                    } else return 0.0;
                }
                return 0.0;
                });

            {
                auto C_view = C.view_device();
                auto sampler_view = sampler.view_device();
                auto x_view = this->space_.assembly_view_device(x);

                Dev::parallel_for(
                    this->space_.element_range(), UTOPIA_LAMBDA(const SizeType &i) {
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
        Scalar ymin_, ymax_;
        bool random_damage_, uniform_damage_, central_damage_;
        Scalar number_;
    };



    template <class FunctionSpace>
    class HomogeneousBar: public InitialCondition<FunctionSpace> {
    public:
        // using Comm           = typename FunctionSpace::Comm;
        using Mesh = typename FunctionSpace::Mesh;
        using Elem = typename FunctionSpace::Shape;
        using ElemView = typename FunctionSpace::ViewDevice::Elem;
        using SizeType = typename FunctionSpace::SizeType;
        using Scalar = typename FunctionSpace::Scalar;
        using Dev = typename FunctionSpace::Device;
        using Point = typename FunctionSpace::Point;
        using ElemViewScalar = typename utopia::FunctionSpace<Mesh, 1, Elem>::ViewDevice::Elem;
        static const int NNodes = Elem::NNodes;

        HomogeneousBar(FunctionSpace &space, const SizeType &PF_component)
            : InitialCondition<FunctionSpace>(space), PF_component_(PF_component){}

        void read(Input & ) override {
        }

        void init(PetscVector &x) override {
            using CoeffVector = utopia::StaticVector<Scalar, NNodes>;
            // un-hard-code
            auto C = this->space_.subspace(PF_component_);

            auto sampler = utopia::sampler(C, [](const Point &) -> Scalar {
                return 0.0;
                });

            {
                auto C_view = C.view_device();
                auto sampler_view = sampler.view_device();
                auto x_view = this->space_.assembly_view_device(x);

                Dev::parallel_for(
                    this->space_.element_range(), UTOPIA_LAMBDA(const SizeType &i) {
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

    };


    template <class FunctionSpace>
    class SmoothQuadraticPhaseField: public InitialCondition<FunctionSpace> {
    public:
        // using Comm           = typename FunctionSpace::Comm;
        using Mesh = typename FunctionSpace::Mesh;
        using Elem = typename FunctionSpace::Shape;
        using ElemView = typename FunctionSpace::ViewDevice::Elem;
        using SizeType = typename FunctionSpace::SizeType;
        using Scalar = typename FunctionSpace::Scalar;
        using Dev = typename FunctionSpace::Device;
        using Point = typename FunctionSpace::Point;
        using ElemViewScalar = typename utopia::FunctionSpace<Mesh, 1, Elem>::ViewDevice::Elem;
        static const int NNodes = Elem::NNodes;

        SmoothQuadraticPhaseField(FunctionSpace &space, const SizeType &PF_component)
            : InitialCondition<FunctionSpace>(space), PF_component_(PF_component){}

        void read(Input & ) override {
        }

        void init(PetscVector &x) override {
            using CoeffVector = utopia::StaticVector<Scalar, NNodes>;
            // un-hard-code
            auto C = this->space_.subspace(PF_component_);

            auto xyz_min = this->space_.mesh().box_min();
            auto xyz_max = this->space_.mesh().box_max();


            auto sampler = utopia::sampler(C, [xyz_max](const Point &x) -> Scalar {
                double norm = xyz_max[0]*xyz_max[0] + xyz_max[1]*xyz_max[1];
                double alpha = x[0]*x[0] + x[1]*x[1];
                std::cout << alpha/norm << std::endl;
                return alpha/norm;
                });

            {
                auto C_view = C.view_device();
                auto sampler_view = sampler.view_device();
                auto x_view = this->space_.assembly_view_device(x);

                Dev::parallel_for(
                    this->space_.element_range(), UTOPIA_LAMBDA(const SizeType &i) {
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

    };

    template <class FunctionSpace>
    class SmoothLinearPhaseField: public InitialCondition<FunctionSpace> {
    public:
        // using Comm           = typename FunctionSpace::Comm;
        using Mesh = typename FunctionSpace::Mesh;
        using Elem = typename FunctionSpace::Shape;
        using ElemView = typename FunctionSpace::ViewDevice::Elem;
        using SizeType = typename FunctionSpace::SizeType;
        using Scalar = typename FunctionSpace::Scalar;
        using Dev = typename FunctionSpace::Device;
        using Point = typename FunctionSpace::Point;
        using ElemViewScalar = typename utopia::FunctionSpace<Mesh, 1, Elem>::ViewDevice::Elem;
        static const int NNodes = Elem::NNodes;

        SmoothLinearPhaseField(FunctionSpace &space, const SizeType &PF_component)
            : InitialCondition<FunctionSpace>(space), PF_component_(PF_component){}

        void read(Input & ) override {
        }

        void init(PetscVector &x) override {
            using CoeffVector = utopia::StaticVector<Scalar, NNodes>;
            // un-hard-code
            auto C = this->space_.subspace(PF_component_);

            auto sampler = utopia::sampler(C, [](const Point &x) -> Scalar {
                double alpha = x[0] + x[1];
                return alpha;
                });

            {
                auto C_view = C.view_device();
                auto sampler_view = sampler.view_device();
                auto x_view = this->space_.assembly_view_device(x);

                Dev::parallel_for(
                    this->space_.element_range(), UTOPIA_LAMBDA(const SizeType &i) {
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

    };



    template <class FunctionSpace>
    class FracPlateIC : public InitialCondition<FunctionSpace> {
    public:
        // using Comm           = typename FunctionSpace::Comm;
        using Mesh = typename FunctionSpace::Mesh;
        using Elem = typename FunctionSpace::Shape;
        using ElemView = typename FunctionSpace::ViewDevice::Elem;
        using SizeType = typename FunctionSpace::SizeType;
        using Scalar = typename FunctionSpace::Scalar;
        using Dev = typename FunctionSpace::Device;
        using Point = typename FunctionSpace::Point;
        using ElemViewScalar = typename utopia::FunctionSpace<Mesh, 1, Elem>::ViewDevice::Elem;
        static const int NNodes = Elem::NNodes;

        using Coord = Point2D<Scalar>;

        FracPlateIC(FunctionSpace &space, const SizeType &PF_component)
            : InitialCondition<FunctionSpace>(space), PF_component_(PF_component) {}

        void init(PetscVector &x) override {
            using CoeffVector = utopia::StaticVector<Scalar, NNodes>;
            // un-hard-code
            auto C = this->space_.subspace(PF_component_);

            auto width = 3.0 * this->space_.mesh().min_spacing();
            // auto width = 0.1;

            if (mpi_world_rank() == 0) {
                utopia::out() << "width: " << width << "  \n";
            }

            std::vector<Rectangle<Scalar>> rectangles;

            rectangles.push_back(Rectangle<Scalar>(Coord(0.308514, 1.531184), Coord(0.488788, 1.711458), width));
            rectangles.push_back(Rectangle<Scalar>(Coord(0.605291, 1.511563), Coord(0.713210, 1.332516), width));
            rectangles.push_back(Rectangle<Scalar>(Coord(1.128942, 1.518921), Coord(1.359495, 1.694289), width));
            rectangles.push_back(Rectangle<Scalar>(Coord(1.517694, 1.412228), Coord(1.673441, 1.229502), width));
            rectangles.push_back(Rectangle<Scalar>(Coord(0.268045, 0.819902), Coord(0.411528, 0.991591), width));
            rectangles.push_back(Rectangle<Scalar>(Coord(0.829713, 0.903294), Coord(1.087246, 0.957253), width));
            rectangles.push_back(Rectangle<Scalar>(Coord(1.456377, 0.994043), Coord(1.592502, 0.819902), width));
            rectangles.push_back(Rectangle<Scalar>(Coord(0.326910, 0.368605), Coord(0.482656, 0.520673), width));
            rectangles.push_back(Rectangle<Scalar>(Coord(0.908199, 0.465487), Coord(1.090925, 0.346531), width));
            rectangles.push_back(Rectangle<Scalar>(Coord(1.436755, 0.364926), Coord(1.624387, 0.493693), width));

            auto sampler = utopia::sampler(C, [&rectangles](const Point &x) -> Scalar {
                for (std::size_t r = 0; r < rectangles.size(); r++) {
                    if (rectangles[r].belongs_to_rectangle(x[0], x[1])) return 1.0;
                }
                return 0.0;
            });

            {
                auto C_view = C.view_device();
                auto sampler_view = sampler.view_device();
                auto x_view = this->space_.assembly_view_device(x);

                Dev::parallel_for(
                    this->space_.element_range(), UTOPIA_LAMBDA(const SizeType &i) {
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
    };

}  // namespace utopia

#endif

#endif
