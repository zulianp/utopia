#ifndef UTOPIA_BC_SETUP_HPP
#define UTOPIA_BC_SETUP_HPP

#include "utopia_Base.hpp"
#include "utopia_RangeDevice.hpp"

// include edsl components
#include "utopia_AppRunner.hpp"
#include "utopia_AssemblyView.hpp"
#include "utopia_Backtracking.hpp"
#include "utopia_BratuFE.hpp"
#include "utopia_ConjugateGradient.hpp"
#include "utopia_Core.hpp"
#include "utopia_DeviceView.hpp"
#include "utopia_FEFunction.hpp"
#include "utopia_GradInterpolate.hpp"
#include "utopia_InitialCondition.hpp"
#include "utopia_IsotropicPhaseField.hpp"
#include "utopia_LBFGS.hpp"
#include "utopia_LaplacianView.hpp"
#include "utopia_LinearElasticityView.hpp"
#include "utopia_MPITimeStatistics.hpp"
#include "utopia_MPRGP.hpp"
#include "utopia_MassMatrixView.hpp"
#include "utopia_PoissonFE.hpp"
#include "utopia_PrincipalStrainsView.hpp"
#include "utopia_QuasiNewtonBound.hpp"
#include "utopia_QuasiTrustRegionVariableBound.hpp"
#include "utopia_SampleView.hpp"
#include "utopia_TrivialPreconditioners.hpp"
#include "utopia_TrustRegionVariableBound.hpp"

#include "utopia_petsc.hpp"
#include "utopia_petsc_DM.hpp"
#include "utopia_petsc_DMDA_FunctionSpace.hpp"
#include "utopia_petsc_DirichletBoundaryConditions.hpp"
#include "utopia_petsc_Matrix.hpp"
#include "utopia_petsc_impl.hpp"

#include "utopia_petsc_Constraints.hpp"

#include <chrono>
#include <cmath>
#include <random>

namespace utopia {

    template <class FunctionSpace>
    class BCSetup : virtual public Configurable {
    public:
        using Scalar = typename FunctionSpace::Scalar;

        BCSetup(FunctionSpace &space) : space_(space) {}

        ~BCSetup() override = default;

        void read(Input &) override {}

        virtual void emplace_BC(){};
        virtual void emplace_time_dependent_BC(const Scalar & /*time*/){};

    protected:
        FunctionSpace &space_;
    };

    template <class FunctionSpace>
    class PFFracFixAllDisp : public BCSetup<FunctionSpace> {
    public:
        using Scalar = typename FunctionSpace::Scalar;
        using Vector = typename FunctionSpace::Vector;

        PFFracFixAllDisp(FunctionSpace &space) : BCSetup<FunctionSpace>(space) {}

        void emplace_time_dependent_BC(const Scalar & /*time*/) override {
            static const int Dim = FunctionSpace::Dim;

            using Point = typename FunctionSpace::Point;
            this->space_.reset_bc();

            for (int d = 1; d < Dim + 1; ++d) {
                // this->space_.emplace_dirichlet_condition(
                //     SideSet::left(), UTOPIA_LAMBDA(const Point &)->Scalar { return 0.0;
                //     }, d);

                // this->space_.emplace_dirichlet_condition(
                //     SideSet::right(), UTOPIA_LAMBDA(const Point &)->Scalar { return
                //     0.0; }, d);

                if (d == 1) {
                    this->space_.emplace_dirichlet_condition(
                        SideSet::left(), UTOPIA_LAMBDA(const Point &)->Scalar { return 0.000; }, d);
                } else {
                    this->space_.emplace_dirichlet_condition(
                        SideSet::left(), UTOPIA_LAMBDA(const Point &)->Scalar { return 0.0; }, d);
                }

                if (d == 1) {
                    this->space_.emplace_dirichlet_condition(
                        SideSet::right(), UTOPIA_LAMBDA(const Point &)->Scalar { return 0.001; }, d);
                } else {
                    this->space_.emplace_dirichlet_condition(
                        SideSet::right(), UTOPIA_LAMBDA(const Point &)->Scalar { return 0.0; }, d);
                }

                if (Dim == 3) {
                    this->space_.emplace_dirichlet_condition(
                        SideSet::front(), UTOPIA_LAMBDA(const Point &)->Scalar { return 0.0; }, d);

                    this->space_.emplace_dirichlet_condition(
                        SideSet::back(), UTOPIA_LAMBDA(const Point &)->Scalar { return 0.0; }, d);
                }
            }
        }
    };

    template <class FunctionSpace>
    class PFFracFixAllDisp4Sides : public BCSetup<FunctionSpace> {
    public:
        using Scalar = typename FunctionSpace::Scalar;
        using Vector = typename FunctionSpace::Vector;

        PFFracFixAllDisp4Sides(FunctionSpace &space) : BCSetup<FunctionSpace>(space) {}

        void emplace_time_dependent_BC(const Scalar & /*time*/) override {
            static const int Dim = FunctionSpace::Dim;

            using Point = typename FunctionSpace::Point;
            this->space_.reset_bc();

            for (int d = 1; d < Dim + 1; ++d) {
                this->space_.emplace_dirichlet_condition(
                    SideSet::left(), UTOPIA_LAMBDA(const Point &)->Scalar { return 0.0; }, d);

                this->space_.emplace_dirichlet_condition(
                    SideSet::right(), UTOPIA_LAMBDA(const Point &)->Scalar { return 0.0; }, d);

                this->space_.emplace_dirichlet_condition(
                    SideSet::top(), UTOPIA_LAMBDA(const Point &)->Scalar { return 0.0; }, d);

                this->space_.emplace_dirichlet_condition(
                    SideSet::bottom(), UTOPIA_LAMBDA(const Point &)->Scalar { return 0.0; }, d);
            }
        }
    };

    template <class FunctionSpace>
    class PFFracFixAllDisp3D : public BCSetup<FunctionSpace> {
    public:
        using Scalar = typename FunctionSpace::Scalar;
        using Vector = typename FunctionSpace::Vector;

        PFFracFixAllDisp3D(FunctionSpace &space) : BCSetup<FunctionSpace>(space) {}

        void emplace_time_dependent_BC(const Scalar & /*time*/) override {
            static const int Dim = FunctionSpace::Dim;

            using Point = typename FunctionSpace::Point;
            this->space_.reset_bc();

            for (int d = 1; d < Dim + 1; ++d) {
                this->space_.emplace_dirichlet_condition(
                    SideSet::left(), UTOPIA_LAMBDA(const Point &)->Scalar { return 0.0; }, d);

                this->space_.emplace_dirichlet_condition(
                    SideSet::right(), UTOPIA_LAMBDA(const Point &)->Scalar { return 0.0; }, d);

                this->space_.emplace_dirichlet_condition(
                    SideSet::front(), UTOPIA_LAMBDA(const Point &)->Scalar { return 0.0; }, d);

                this->space_.emplace_dirichlet_condition(
                    SideSet::back(), UTOPIA_LAMBDA(const Point &)->Scalar { return 0.0; }, d);

                this->space_.emplace_dirichlet_condition(
                    SideSet::top(), UTOPIA_LAMBDA(const Point &)->Scalar { return 0.0; }, d);

                this->space_.emplace_dirichlet_condition(
                    SideSet::bottom(), UTOPIA_LAMBDA(const Point &)->Scalar { return 0.0; }, d);
            }
        }
    };

    template <class FunctionSpace>
    class PFFracFixAllDispComp2D : public BCSetup<FunctionSpace> {
    public:
        using Scalar = typename FunctionSpace::Scalar;
        using Vector = typename FunctionSpace::Vector;

        PFFracFixAllDispComp2D(FunctionSpace &space) : BCSetup<FunctionSpace>(space) {}

        void emplace_time_dependent_BC(const Scalar & /*time*/) override {
            static const int Dim = FunctionSpace::Dim;

            using Point = typename FunctionSpace::Point;
            this->space_.reset_bc();

            for (int d = 0; d < Dim + 1; ++d) {
                if (d == 1) {
                    this->space_.emplace_dirichlet_condition(
                        SideSet::left(), UTOPIA_LAMBDA(const Point &)->Scalar { return 1e-5; }, d);
                } else {
                    this->space_.emplace_dirichlet_condition(
                        SideSet::left(), UTOPIA_LAMBDA(const Point &)->Scalar { return 0.0; }, d);
                }

                this->space_.emplace_dirichlet_condition(
                    SideSet::top(), UTOPIA_LAMBDA(const Point &)->Scalar { return 0.0; }, d);

                this->space_.emplace_dirichlet_condition(
                    SideSet::right(), UTOPIA_LAMBDA(const Point &)->Scalar { return 0.0; }, d);

                this->space_.emplace_dirichlet_condition(
                    SideSet::bottom(), UTOPIA_LAMBDA(const Point &)->Scalar { return 0.0; }, d);
            }
        }
    };

    template <class FunctionSpace>
    class PFFracTension : public BCSetup<FunctionSpace> {
    public:
        using Scalar = typename FunctionSpace::Scalar;
        using Vector = typename FunctionSpace::Vector;

        PFFracTension(FunctionSpace &space, const Scalar &disp_y = 1.0)
            : BCSetup<FunctionSpace>(space), disp_y_(disp_y) {}

        void read(Input &in) override {
            in.get("disp_y", disp_y_);
            in.get("fix_phase_field_on_sides", fix_phase_field_);
        }

        void emplace_time_dependent_BC(const Scalar &time) override {
            //            static const int Dim = FunctionSpace::Dim;

            using Point = typename FunctionSpace::Point;
            this->space_.reset_bc();

            this->space_.emplace_dirichlet_condition(
                SideSet::bottom(), UTOPIA_LAMBDA(const Point &)->Scalar { return 0.0; }, 1);

            this->space_.emplace_dirichlet_condition(
                SideSet::bottom(), UTOPIA_LAMBDA(const Point &)->Scalar { return 0.0; }, 2);

            this->space_.emplace_dirichlet_condition(
                SideSet::top(), UTOPIA_LAMBDA(const Point &)->Scalar { return disp_y_ * time; }, 2);

            this->space_.emplace_dirichlet_condition(
                SideSet::top(), UTOPIA_LAMBDA(const Point &)->Scalar { return 0.0; }, 1);  // fix x displacement

            // Fixing top and bottom boundary to no have any damage on them
            if (fix_phase_field_) {
                this->space_.emplace_dirichlet_condition(
                    SideSet::bottom(),
                    UTOPIA_LAMBDA(const Point &)->Scalar { return 0.0; },
                    0  // alpha
                );
                this->space_.emplace_dirichlet_condition(
                    SideSet::top(),
                    UTOPIA_LAMBDA(const Point &)->Scalar { return 0.0; },
                    0  // alpha
                );
            }
        }

    private:
        Scalar disp_y_;
        bool fix_phase_field_;
    };

    template <class FunctionSpace>
    class PFFracShear2D : public BCSetup<FunctionSpace> {
    public:
        using Scalar = typename FunctionSpace::Scalar;
        using Vector = typename FunctionSpace::Vector;

        PFFracShear2D(FunctionSpace &space, const Scalar &disp_x = 1.0)
            : BCSetup<FunctionSpace>(space), disp_x_(disp_x) {}

        void read(Input &in) override {
            in.get("disp_x", disp_x_);
            in.get("fix_phase_field_on_sides", fix_phase_field_);
        }

        void emplace_time_dependent_BC(const Scalar &time) override {
            // static const int Dim = FunctionSpace::Dim;

            using Point = typename FunctionSpace::Point;
            this->space_.reset_bc();

            this->space_.emplace_dirichlet_condition(
                SideSet::bottom(), UTOPIA_LAMBDA(const Point &)->Scalar { return 0.0; }, 1);  // fix x

            this->space_.emplace_dirichlet_condition(
                SideSet::bottom(), UTOPIA_LAMBDA(const Point &)->Scalar { return 0.0; }, 2);  // fix y

            this->space_.emplace_dirichlet_condition(
                SideSet::top(), UTOPIA_LAMBDA(const Point &)->Scalar { return disp_x_ * time; }, 1);  // shear x

            this->space_.emplace_dirichlet_condition(
                SideSet::top(), UTOPIA_LAMBDA(const Point &)->Scalar { return 0.0; }, 2);  // fix y

            // Fixing top and bottom boundary to no have any damage on them
            if (fix_phase_field_) {
                this->space_.emplace_dirichlet_condition(
                    SideSet::bottom(),
                    UTOPIA_LAMBDA(const Point &)->Scalar { return 0.0; },
                    0  // alpha
                );
                this->space_.emplace_dirichlet_condition(
                    SideSet::top(),
                    UTOPIA_LAMBDA(const Point &)->Scalar { return 0.0; },
                    0  // alpha
                );
            }
        }

    private:
        Scalar disp_x_;
        bool fix_phase_field_;
    };

    template <class FunctionSpace>
    class AsphaltTension2D : public BCSetup<FunctionSpace> {
    public:
        using Scalar = typename FunctionSpace::Scalar;
        using Vector = typename FunctionSpace::Vector;

        AsphaltTension2D(FunctionSpace &space, const Scalar &disp_y = 1.0)
            : BCSetup<FunctionSpace>(space), disp_y_(disp_y) {}

        void read(Input &in) override { in.get("disp_y", disp_y_); }

        void emplace_time_dependent_BC(const Scalar &time) override {
            // static const int Dim = FunctionSpace::Dim;

            using Point = typename FunctionSpace::Point;
            this->space_.reset_bc();

            this->space_.emplace_dirichlet_condition(
                SideSet::bottom(),
                UTOPIA_LAMBDA(const Point &)->Scalar { return 0.0; },
                1  // disp_x
            );

            this->space_.emplace_dirichlet_condition(
                SideSet::bottom(),
                UTOPIA_LAMBDA(const Point &)->Scalar { return 0.0; },
                2  // disp_y
            );

            this->space_.emplace_dirichlet_condition(
                SideSet::top(),
                UTOPIA_LAMBDA(const Point &)->Scalar { return disp_y_ * time; },
                2  // disp_y
            );
        }

    private:
        Scalar disp_y_;
    };

    template <class FunctionSpace>
    class BiaxialLoading2D : public BCSetup<FunctionSpace> {
    public:
        using Scalar = typename FunctionSpace::Scalar;
        using Vector = typename FunctionSpace::Vector;

        BiaxialLoading2D(FunctionSpace &space, const Scalar &Top_disp_y = 1.0, const Scalar &Right_disp_x = 1.0)
            : BCSetup<FunctionSpace>(space), disp_y_(Top_disp_y), disp_x_(Right_disp_x) {}

        void read(Input &in) override {
            in.get("disp_y", disp_y_);
            in.get("disp_x", disp_x_);
            in.get("fix_phase_field_on_sides", fix_phase_field_);
        }

        void emplace_time_dependent_BC(const Scalar &time) override {
            // static const int Dim = FunctionSpace::Dim;

            using Point = typename FunctionSpace::Point;
            this->space_.reset_bc();

            this->space_.emplace_dirichlet_condition(
                SideSet::bottom(),
                UTOPIA_LAMBDA(const Point &)->Scalar { return 0.0; },
                2  // disp_y
            );

            this->space_.emplace_dirichlet_condition(
                SideSet::left(),
                UTOPIA_LAMBDA(const Point &)->Scalar { return 0.0; },
                1  // disp_x
            );

            this->space_.emplace_dirichlet_condition(
                SideSet::top(),
                UTOPIA_LAMBDA(const Point &)->Scalar { return disp_y_ * time; },
                2  // disp_y
            );

            this->space_.emplace_dirichlet_condition(
                SideSet::right(),
                UTOPIA_LAMBDA(const Point &)->Scalar { return disp_x_ * time; },
                1  // disp_x
            );
        }

    private:
        Scalar disp_y_;
        Scalar disp_x_;
        bool fix_phase_field_;
    };

    template <class FunctionSpace>
    class UniaxialLoading2D : public BCSetup<FunctionSpace> {
    public:
        using Scalar = typename FunctionSpace::Scalar;
        using Vector = typename FunctionSpace::Vector;

        UniaxialLoading2D(FunctionSpace &space, const Scalar &Top_disp_y = 1.0, const Scalar &Right_disp_x = 1.0)
            : BCSetup<FunctionSpace>(space), disp_y_(Top_disp_y), disp_x_(Right_disp_x) {}

        void read(Input &in) override {
            in.get("disp_y", disp_y_);
            in.get("disp_x", disp_x_);
            in.get("fix_phase_field_on_sides", fix_phase_field_);
        }

        void emplace_time_dependent_BC(const Scalar &time) override {
            // static const int Dim = FunctionSpace::Dim;

            using Point = typename FunctionSpace::Point;
            this->space_.reset_bc();

            this->space_.emplace_dirichlet_condition(
                SideSet::left(),
                UTOPIA_LAMBDA(const Point &)->Scalar { return 0.0; },
                2  // disp_y
            );

            this->space_.emplace_dirichlet_condition(
                SideSet::left(),
                UTOPIA_LAMBDA(const Point &)->Scalar { return 0.0; },
                1  // disp_x
            );

            this->space_.emplace_dirichlet_condition(
                SideSet::right(),
                UTOPIA_LAMBDA(const Point &)->Scalar { return disp_x_ * time; },
                1  // disp_x
            );

            // Fixing left and right boundary to no have any damage on them
            if (fix_phase_field_) {
                this->space_.emplace_dirichlet_condition(
                    SideSet::left(),
                    UTOPIA_LAMBDA(const Point &)->Scalar { return 0.0; },
                    0  // alpha
                );

                this->space_.emplace_dirichlet_condition(
                    SideSet::right(),
                    UTOPIA_LAMBDA(const Point &)->Scalar { return 0.0; },
                    0  // alpha
                );
            }
        }

    private:
        Scalar disp_y_;
        Scalar disp_x_;
        bool fix_phase_field_{false};
    };

    template <class FunctionSpace>
    class SedimentaryLayers : public BCSetup<FunctionSpace> {
    public:
        using Scalar = typename FunctionSpace::Scalar;
        using Vector = typename FunctionSpace::Vector;

        SedimentaryLayers(FunctionSpace &space, const Scalar &Top_disp_y = 1.0, const Scalar &Right_disp_x = 1.0)
            : BCSetup<FunctionSpace>(space), disp_y_(Top_disp_y), disp_x_(Right_disp_x) {}

        void read(Input &in) override {
            in.get("disp_y", disp_y_);
            in.get("disp_x", disp_x_);
            in.get("fix_phase_field_on_sides", fix_phase_field_);
            in.get("set_damage_to_zero_in_softer_layers", set_damage_to_zero_in_softer_layers_ );
        }

        void emplace_time_dependent_BC(const Scalar &time) override {
            // static const int Dim = FunctionSpace::Dim;

            using Point = typename FunctionSpace::Point;
            this->space_.reset_bc();

            this->space_.emplace_dirichlet_condition(
                SideSet::left(),
                UTOPIA_LAMBDA(const Point &)->Scalar { return 0.0; },
                1  // disp_x
            );

            this->space_.emplace_dirichlet_condition(
                SideSet::bottom(),
                UTOPIA_LAMBDA(const Point &)->Scalar { return 0.0; },
                2  // disp_y
            );

            this->space_.emplace_dirichlet_condition(
                SideSet::right(),
                UTOPIA_LAMBDA(const Point &)->Scalar { return disp_x_ * time; },
                1  // disp_x
            );

            // Fixing left and right boundary to no have any damage on them
            if (fix_phase_field_) {
                this->space_.emplace_dirichlet_condition(
                    SideSet::left(),
                    UTOPIA_LAMBDA(const Point &)->Scalar { return 0.0; },
                    0  // alpha
                );

                this->space_.emplace_dirichlet_condition(
                    SideSet::right(),
                    UTOPIA_LAMBDA(const Point &)->Scalar { return 0.0; },
                    0  // alpha
                );
            }
        }

    private:
        Scalar disp_y_;
        Scalar disp_x_;
        bool   fix_phase_field_{false};
        bool   set_damage_to_zero_in_softer_layers_{false};
    };

    template <class FunctionSpace>
    class PFMixed2D : public BCSetup<FunctionSpace> {
    public:
        using Scalar = typename FunctionSpace::Scalar;
        using Vector = typename FunctionSpace::Vector;

        PFMixed2D(FunctionSpace &space, const Scalar &disp_y = 1.0) : BCSetup<FunctionSpace>(space), disp_y_(disp_y) {}

        void read(Input &in) override { in.get("disp_y", disp_y_); }

        void emplace_time_dependent_BC(const Scalar &time) override {
            // static const int Dim = FunctionSpace::Dim;

            using Point = typename FunctionSpace::Point;
            this->space_.reset_bc();

            this->space_.emplace_dirichlet_condition(
                SideSet::bottom(),
                UTOPIA_LAMBDA(const Point &)->Scalar { return 0.0; },
                1  // disp_x
            );

            this->space_.emplace_dirichlet_condition(
                SideSet::bottom(),
                UTOPIA_LAMBDA(const Point &)->Scalar { return 0.0; },
                2  // disp_y
            );

            this->space_.emplace_dirichlet_condition(
                SideSet::top(),
                UTOPIA_LAMBDA(const Point &)->Scalar { return disp_y_ * time; },
                2  // disp_y
            );
        }

    private:
        Scalar disp_y_;
    };

    template <class FunctionSpace>
    class FracPlateBC : public BCSetup<FunctionSpace> {
    public:
        using Scalar = typename FunctionSpace::Scalar;
        using Vector = typename FunctionSpace::Vector;

        FracPlateBC(FunctionSpace &space, const Scalar &disp_x = 1.0, const Scalar &disp_y = 1.0)
            : BCSetup<FunctionSpace>(space), disp_x_(disp_x), disp_y_(disp_y) {}

        void read(Input &in) override {
            in.get("disp_x", disp_x_);
            in.get("disp_y", disp_y_);
        }

        void emplace_time_dependent_BC(const Scalar &time) override {
            // static const int Dim = FunctionSpace::Dim;

            using Point = typename FunctionSpace::Point;
            this->space_.reset_bc();

            // this->space_.emplace_dirichlet_condition(
            //     SideSet::bottom(),
            //     UTOPIA_LAMBDA(const Point &p) -> Scalar {
            //         return 0.0;
            //     },
            //     1 // disp_x
            //     );

            this->space_.emplace_dirichlet_condition(
                SideSet::bottom(),
                UTOPIA_LAMBDA(const Point &)->Scalar { return 0.0; },
                2  // disp_y
            );

            this->space_.emplace_dirichlet_condition(
                SideSet::left(),
                UTOPIA_LAMBDA(const Point &)->Scalar { return 0.0; },
                1  // disp_x
            );

            // this->space_.emplace_dirichlet_condition(
            //     SideSet::left(),
            //     UTOPIA_LAMBDA(const Point &p) -> Scalar {
            //         return 0.0;
            //     },
            //     2 // disp_y
            //     );

            /////////////////////////////////////////////////////////////////

            // this->space_.emplace_dirichlet_condition(
            //     SideSet::top(),
            //     UTOPIA_LAMBDA(const Point &p) -> Scalar {
            //         return 0.0;
            //     },
            //     1 // disp_x
            //     );

            this->space_.emplace_dirichlet_condition(
                SideSet::top(),
                UTOPIA_LAMBDA(const Point &)->Scalar { return disp_y_ * time; },
                2  // disp_y
            );

            this->space_.emplace_dirichlet_condition(
                SideSet::right(),
                UTOPIA_LAMBDA(const Point &)->Scalar { return disp_x_ * time; },
                1  // disp_x
            );

            // this->space_.emplace_dirichlet_condition(
            //     SideSet::right(),
            //     UTOPIA_LAMBDA(const Point &p) -> Scalar {
            //         return 0.0;
            //     },
            //     2 // disp_y
            //     );
        }

    private:
        Scalar disp_x_;
        Scalar disp_y_;
    };

    template <class FunctionSpace>
    class FixedSubdomain2D : public BCSetup<FunctionSpace> {
    public:
        using Scalar = typename FunctionSpace::Scalar;
        using Vector = typename FunctionSpace::Vector;

        FixedSubdomain2D(FunctionSpace &space,
                         const Scalar &x_min = 0,
                         const Scalar &y_min = 0,
                         const Scalar &x_max = 1,
                         const Scalar &y_max = 1,
                         const Scalar &val = 0)
            : BCSetup<FunctionSpace>(space), x_min_(x_min), y_min_(y_min), x_max_(x_max), y_max_(y_max), val_(val) {}

        void read(Input &in) override {
            in.get("fixed_x_min", x_min_);
            in.get("fixed_y_min", y_min_);
            in.get("fixed_x_max", x_max_);
            in.get("fixed_y_max", y_max_);
        }

        void emplace_time_dependent_BC(const Scalar &) override {
            using Point = typename FunctionSpace::Point;
            this->space_.reset_constraints();

            this->space_.emplace_constraint(
                UTOPIA_LAMBDA(const Point &p)->bool {
                    return p(0) >= x_min_ && p(1) >= y_min_ && p(0) <= x_max_ && p(1) <= y_max_;
                },
                UTOPIA_LAMBDA(const Point &)->Scalar { return val_; },
                0);
        }

    private:
        Scalar x_min_, y_min_;
        Scalar x_max_, y_max_;
        Scalar val_;
    };

    template <class FunctionSpace>
    class LayeredSubdomain : public BCSetup<FunctionSpace> {
    public:
        using Scalar = typename FunctionSpace::Scalar;
        using Vector = typename FunctionSpace::Vector;

        LayeredSubdomain(FunctionSpace &space,
                         const Scalar &y_min = 0,
                         const Scalar &y_max = 1,
                         const Scalar &val = 0)
            : BCSetup<FunctionSpace>(space), y_min_(y_min), y_max_(y_max), val_(val) {}

        void read(Input &in) override {
            in.get("bottom_layer_height", y_min_);
            in.get("top_layer_height"   , y_max_);
            in.get("bottom_layer_height2", y_min2_);
            in.get("top_layer_height2"   , y_max2_);
        }

        void emplace_time_dependent_BC(const Scalar &) override {
            using Point = typename FunctionSpace::Point;
            this->space_.reset_constraints();

            this->space_.emplace_constraint(
                UTOPIA_LAMBDA(const Point &p)->bool {
                    return  ( !(p(1) >= y_min_  && p(1) <= y_max_) &&   //Not in layer 1  and
                              !(p(1) > y_min2_  && p(1) < y_max2_) );    //Not in layer 2
                },
                UTOPIA_LAMBDA(const Point &)->Scalar { return val_; },
                0);
        }

    private:
        Scalar  y_min_, y_max_;
        Scalar  y_min2_, y_max2_;
        Scalar val_;
    };

    template <class FunctionSpace, class DirichletBC, class Constraints>
    class DirichletAndVolConstraints : public BCSetup<FunctionSpace> {
    public:
        using Scalar = typename FunctionSpace::Scalar;
        using Vector = typename FunctionSpace::Vector;

        DirichletAndVolConstraints(FunctionSpace &space)
            : BCSetup<FunctionSpace>(space), bc_(space), constraints_(space) {}

        void read(Input &in) override {
            bc_.read(in);
            constraints_.read(in);
        }

        void emplace_time_dependent_BC(const Scalar &t) override {
            bc_.emplace_time_dependent_BC(t);
            constraints_.emplace_time_dependent_BC(t);
        }

    private:
        DirichletBC bc_;
        Constraints constraints_;
    };

}  // namespace utopia

#endif  // UTOPIA_BC_SETUP_HPP
