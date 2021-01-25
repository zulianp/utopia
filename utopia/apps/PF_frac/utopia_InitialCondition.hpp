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

namespace utopia {

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

        void init(PetscVector &x) override {
            // un-hard-code
            auto C = this->space_.subspace(PF_component_);

            auto sampler = utopia::sampler(
                C, UTOPIA_LAMBDA(const Point &x)->Scalar {
                    Scalar f = 0.0;
                    if (x[0] > (0.5 - this->space_.mesh().min_spacing()) &&
                        x[0] < (0.5 + this->space_.mesh().min_spacing()) && x[1] < 0.25) {
                        // if (x[0] <= 0.5 && x[1] <= 0.5) {
                        // f = 1.0;
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
                        x[0] < (0.5 + this->space_.mesh().min_spacing()) && x[1] < 0.5 && x[1] > 0.3) {
                        f = 1.0;
                    } else if ((x[0] > 0.3) && (x[0] < 0.7) && (x[1] > 0.7 - this->space_.mesh().min_spacing()) &&
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

        InitialCondidtionPFSneddon(FunctionSpace &space, const SizeType &PF_component)
            : InitialCondition<FunctionSpace>(space), PF_component_(PF_component) {}

        void init(PetscVector &x) override {
            // un-hard-code
            auto C = this->space_.subspace(PF_component_);

            auto sampler = utopia::sampler(
                C, UTOPIA_LAMBDA(const Point &x)->Scalar {
                    Scalar f = 0.0;
                    Scalar h = this->space_.mesh().min_spacing();

                    if (x[0] > 0.55 && x[0] < 0.7 && x[1] > (0.4) && x[1] < (0.6) && x[2] < (0.6 + h) &&
                        x[2] > (0.6 - h)) {
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
