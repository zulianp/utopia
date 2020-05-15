#ifndef UTOPIA_FRAC_NET_GENERATOR3D_PF_HPP
#define UTOPIA_FRAC_NET_GENERATOR3D_PF_HPP

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
#include "utopia_GradInterpolate.hpp"
#include "utopia_InitialConditionBase.hpp"
#include "utopia_IsotropicPhaseField.hpp"
#include "utopia_LinearElasticityView.hpp"
#include "utopia_MPITimeStatistics.hpp"
#include "utopia_MassMatrixView.hpp"
#include "utopia_PetscDM.hpp"
#include "utopia_PrincipalStrainsView.hpp"
#include "utopia_SampleView.hpp"
#include "utopia_TrivialPreconditioners.hpp"
#include "utopia_petsc.hpp"
#include "utopia_petsc_DirichletBoundaryConditions.hpp"
#include "utopia_petsc_Matrix.hpp"
#include "utopia_petsc_dma_FunctionSpace.hpp"

#include <chrono>
#include <cmath>
#include <random>

namespace utopia {

    template <class T>
    class FracNetSamplerParams3D : public FracNetSamplerParams2D<T> {
    public:
        FracNetSamplerParams3D()
            : z_min(0.0),
              z_max(1.0),
              beta_min(0.0),
              beta_max(180),
              uni_min(0.0),
              uni_max(0.0),
              uniform_width(0.0),
              uniform_length(0.0),
              linear_relation(0.0),
              pow_dist_coef_width(2.8),
              min_width(0.1) {}

        void read(Input &in) override {
            FracNetSamplerParams2D<T>::read(in);

            in.get("z_min", z_min);
            in.get("z_max", z_max);

            in.get("beta_min", beta_min);
            in.get("beta_max", beta_max);

            in.get("uni_min", uni_min);
            in.get("uni_max", uni_max);

            in.get("uniform_width", uniform_width);
            in.get("uniform_length", uniform_length);

            in.get("linear_relation", linear_relation);

            in.get("pow_dist_coef_width", pow_dist_coef_width);
            in.get("min_width", min_width);
        }

    public:
        T z_min;
        T z_max;

        T beta_min;
        T beta_max;

        T uni_min;
        T uni_max;

        T uniform_width;
        T uniform_length;

        T linear_relation;

        T pow_dist_coef_width;
        T min_width;
    };

    template <class T>
    struct Point3D {
        Point3D() = default;

        Point3D(const T &xx, const T &yy, const T &&zz) : x(xx), y(yy), z(zz) {}

        T x;
        T y;
        T z;

        void describe() {
            std::cout << "( " << x << " , " << y << " , " << z << " )"
                      << " \n";
        }
    };

    template <class T>
    class Paralleloid {
    public:
        Paralleloid(const T &width, const FracNetSamplerParams3D<T> &params) {
            points_.resize(8);
            randomly_generate(width, params);
        }

        Paralleloid(const Point3D<T> &A, const T &length, const T &width, const T &theta, const T &gamma) {
            points_.resize(8);
            points_[0] = A;
            this->generate_paralleloid(length, width, theta, gamma);
        }

        bool belongs_to_paralleloid(const T &x_coord, const T &y_coord, const T &z_coord) {
            Point3D<T> M{};
            M.x = x_coord;
            M.y = y_coord;
            M.z = z_coord;

            return belongs_to_paralleloid(M);
        }

        bool belongs_to_paralleloid(Point3D<T> M) {
            Point3D<T> u{}, v{}, w{};
            build_vector(points_[0], points_[4], u);
            build_vector(points_[0], points_[1], v);
            build_vector(points_[0], points_[2], w);

            const T dotMu = vec_dot(M, u);
            const T dotMv = vec_dot(M, v);
            const T dotMw = vec_dot(M, w);

            bool flg1 = (vec_dot(u, points_[4]) <= dotMu) && (dotMu <= vec_dot(u, points_[0]));
            bool flg2 = (vec_dot(v, points_[1]) <= dotMv) && (dotMv <= vec_dot(v, points_[0]));
            bool flg3 = (vec_dot(w, points_[2]) <= dotMw) && (dotMw <= vec_dot(w, points_[0]));

            return flg1 && flg2 && flg3;
        }

        void describe() {
            std::cout << "A: " << points_[0].x << " " << points_[0].y << " " << points_[0].z << " \n";
            std::cout << "B: " << points_[1].x << " " << points_[1].y << " " << points_[1].z << " \n";
            std::cout << "C: " << points_[2].x << " " << points_[2].y << " " << points_[2].z << " \n";
            std::cout << "D: " << points_[3].x << " " << points_[3].y << " " << points_[3].z << " \n \n";

            std::cout << "E: " << points_[4].x << " " << points_[4].y << "  " << points_[4].z << "  \n";
            std::cout << "F: " << points_[5].x << " " << points_[5].y << "  " << points_[5].z << "  \n";
            std::cout << "G: " << points_[6].x << " " << points_[6].y << "  " << points_[6].z << "  \n";
            std::cout << "H: " << points_[7].x << " " << points_[7].y << "  " << points_[7].z << "  \n";
            std::cout << "------------------------------  \n";
        }

    private:
        void randomly_generate(const T &depth, const FracNetSamplerParams3D<T> &params) {
            static std::default_random_engine generator(params.seed);

            std::uniform_real_distribution<> distr_point_x(params.uni_min, params.uni_max);
            std::uniform_real_distribution<> distr_point_y(params.uni_min, params.uni_max);
            std::uniform_real_distribution<> distr_point_z(params.uni_min, params.uni_max);

            std::uniform_int_distribution<> distr_angle_alpha(params.alpha_min, params.alpha_max);
            std::uniform_int_distribution<> distr_angle_beta(params.beta_min, params.beta_max);

            points_[0].x = distr_point_x(generator);
            points_[0].y = distr_point_y(generator);
            points_[0].z = distr_point_z(generator);

            const T alpha = distr_angle_alpha(generator);
            const T beta = distr_angle_beta(generator);

            auto max_l = std::max(std::max(params.x_min, params.y_min), params.z_min);
            auto min_l = std::min(std::min(params.x_max, params.y_max), params.z_max);

            // std::uniform_real_distribution<> distr_length(max_l, min_l);
            std::uniform_real_distribution<> distr_length(0, 1);

            T length;
            if (params.uniform_length == 0) {
                const T x_min = (params.min_length == 0) ? 3.0 * depth : params.min_length;
                const T r = distr_length(generator);
                length = x_min * std::pow((1. - r), (-1. / (params.pow_dist_coef_length - 1.)));
            } else {
                length = params.uniform_length;
            }

            T width;
            if (params.uniform_width == 0 && params.linear_relation == 0) {
                const T x_min = (params.min_width == 0) ? 3.0 * depth : params.min_width;
                const T r = distr_length(generator);
                width = x_min * std::pow((1. - r), (-1. / (params.pow_dist_coef_width - 1.)));
            } else if (params.linear_relation != 0) {
                std::cout << "linear_relation" << params.linear_relation << std::endl;

                width = params.linear_relation * length;
            } else {
                width = params.uniform_width;
            }

            std::uniform_int_distribution<> distr_dir(0.0, 6);
            const int i = distr_dir(generator);

            if (i == 0) {
                generate_paralleloid(length, depth, width, alpha, beta, params);
            } else if (i == 1) {
                generate_paralleloid(length, width, depth, alpha, beta, params);
            } else if (i == 3) {
                generate_paralleloid(depth, length, width, alpha, beta, params);
            } else if (i == 4) {
                generate_paralleloid(depth, width, length, alpha, beta, params);
            } else if (i == 5) {
                generate_paralleloid(length, depth, width, alpha, beta, params);
            } else {
                generate_paralleloid(width, depth, length, alpha, beta, params);
            }
        }

        void generate_paralleloid(const T &a,
                                  const T &b,
                                  const T &c,
                                  const T &theta,
                                  const T &gamma,
                                  const FracNetSamplerParams3D<T> &params) {
            const T pi = std::acos(-1.0);
            const T theta_rad = theta * pi / 180.0;
            const T gamma_rad = gamma * pi / 180.0;

            const auto cos_theta = std::cos(theta_rad);
            const auto sin_theta = std::sin(theta_rad);

            const auto cos_theta_pi = std::cos(theta_rad + pi / 2.0);
            const auto sin_theta_pi = std::sin(theta_rad + pi / 2.0);

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

            if (gamma_rad == 0) {
                // as no rotation is required
                return;
            }

            // rotate
            const auto cos_gamma = std::cos(gamma_rad);
            const auto sin_gamma = std::sin(gamma_rad);

            T max_x = 0;
            T max_z = 0;
            T min_x = 0;
            T min_z = 0;

            for (auto p = 0; p < points_.size(); p++) {
                auto ax = points_[p].x;
                auto az = points_[p].z;
                points_[p].x = (ax * cos_gamma) - (az * sin_gamma);
                points_[p].z = (az * cos_gamma) + (ax * sin_gamma);

                max_x = (max_x > points_[p].x) ? max_x : points_[p].x;
                max_z = (max_z > points_[p].z) ? max_z : points_[p].z;

                min_x = (min_x < points_[p].x) ? min_x : points_[p].x;
                min_z = (min_z < points_[p].z) ? min_z : points_[p].z;
            }

            static std::default_random_engine generator(params.seed);

            if (max_x > params.x_max) {
                const auto move_by_p = max_x - params.x_max;
                const auto move_by_m = max_x - params.x_min;

                std::uniform_real_distribution<> distr_move_by(move_by_p, move_by_m);
                auto moving_factor = distr_move_by(generator);

                for (auto p = 0; p < points_.size(); p++) {
                    points_[p].x = points_[p].x - moving_factor;
                }
            }

            if (min_x < params.x_min) {
                const auto move_by_p = params.x_max - min_x;
                const auto move_by_m = params.x_min - min_x;

                std::uniform_real_distribution<> distr_move_by(move_by_m, move_by_p);
                auto moving_factor = distr_move_by(generator);

                for (auto p = 0; p < points_.size(); p++) {
                    points_[p].x = points_[p].x + moving_factor;
                }
            }

            if (max_z > params.z_max) {
                const auto move_by_p = max_z - params.z_max;
                const auto move_by_m = max_z - params.z_min;

                std::uniform_real_distribution<> distr_move_by(move_by_p, move_by_m);
                auto moving_factor = distr_move_by(generator);

                for (auto p = 0; p < points_.size(); p++) {
                    points_[p].z = points_[p].z - moving_factor;
                }
            }

            if (min_z < params.z_min) {
                const auto move_by_p = params.z_max - min_z;
                const auto move_by_m = params.z_min - min_z;

                std::uniform_real_distribution<> distr_move_by(move_by_m, move_by_p);
                auto moving_factor = distr_move_by(generator);

                for (auto p = 0; p < points_.size(); p++) {
                    points_[p].z = points_[p].z + moving_factor;
                }
            }
        }

        void build_vector(const Point3D<T> &A, const Point3D<T> &B, Point3D<T> &result) {
            result.x = A.x - B.x;
            result.y = A.y - B.y;
            result.z = A.z - B.z;
        }

        T vec_dot(const Point3D<T> &A, const Point3D<T> &B) { return (A.x * B.x) + (A.y * B.y) + (A.z * B.z); }

        void vec_cross(const Point3D<T> &A, const Point3D<T> &B, Point3D<T> &result) {
            result.x = A.y * B.z - A.z * B.y;
            result.y = A.z * B.x - A.x * B.z;
            result.z = A.x * B.y - A.y * B.x;
        }

    private:
        std::vector<Point3D<T>> points_;
    };

    template <class FunctionSpace>
    class InitialCondidtionPFFracNet3D : public InitialCondition<FunctionSpace> {
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

        InitialCondidtionPFFracNet3D(FunctionSpace &space, const SizeType &PF_component, const SizeType &num_fracs = 10)
            : InitialCondition<FunctionSpace>(space),
              PF_component_(PF_component),
              num_fracs_(num_fracs)  //, pressure0_(1.0)
        {}

        void read(Input &in) override {
            in.get("num_fracs", num_fracs_);
            // in.get("pressure0", pressure0_);
            sampler_params_.read(in);
        }

        void init(PetscVector &x) override {
            using CoeffVector = utopia::StaticVector<Scalar, NNodes>;
            // un-hard-code
            auto C = this->space_.subspace(PF_component_);

            auto width = 3.0 * this->space_.mesh().min_spacing();

            if (mpi_world_rank() == 0) {
                std::cout << "width: " << width << "  \n";
            }

            std::vector<Paralleloid<Scalar>> paralleloids;

            for (auto r = 0; r < num_fracs_; r++) {
                paralleloids.push_back(Paralleloid<Scalar>(width, sampler_params_));
            }

            auto sampler = utopia::sampler(C, [&paralleloids](const Point &x) -> Scalar {
                for (auto r = 0; r < paralleloids.size(); r++) {
                    if (paralleloids[r].belongs_to_paralleloid(x[0], x[1], x[2])) return 1.0;
                }
                return 0.0;
            });

            {
                auto C_view = C.view_device();
                auto sampler_view = sampler.view_device();
                auto x_view = this->space_.assembly_view_device(x);

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
        FracNetSamplerParams3D<Scalar> sampler_params_;
    };

}  // namespace utopia

#endif

#endif
