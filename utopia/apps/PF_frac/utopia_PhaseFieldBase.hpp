#ifndef UTOPIA_PHASE_FIELD_BASE_HPP
#define UTOPIA_PHASE_FIELD_BASE_HPP

#include "utopia_CoefStrainView.hpp"
#include "utopia_DeviceTensorContraction.hpp"
#include "utopia_DeviceTensorProduct.hpp"
#include "utopia_DiffController.hpp"
#include "utopia_ExtendedFunction.hpp"
#include "utopia_FEFunction.hpp"
#include "utopia_GradInterpolate.hpp"
#include "utopia_LinearElasticityView.hpp"
#include "utopia_PrincipalShapeStressView.hpp"
#include "utopia_PrincipalStrainsView.hpp"
#include "utopia_StrainView.hpp"
#include "utopia_TensorView4.hpp"
#include "utopia_Tracer.hpp"
#include "utopia_Views.hpp"
#include "utopia_petsc_NeumannBoundaryConditions.hpp"

#include "utopia_make_unique.hpp"

#include "utopia_ui.hpp"

#include <glob.h>
#include <cmath>
#include <cstddef>
#include <cstdio>
#include <fstream>
#include <random>

#define UNROLL_FACTOR 4
#define U_MIN(a, b) ((a) < (b) ? (a) : (b))

namespace utopia {

    template <class FunctionSpace, int Dim = FunctionSpace::Dim>
    class PFFracParameters : public Configurable {
    public:
        using Point = typename FunctionSpace::Point;
        using Scalar = typename FunctionSpace::Scalar;
        using HeteroParamsFunction =
            std::function<void(const Point &, Scalar &, Scalar &, Scalar &, Scalar &, Scalar &, Scalar &)>;

        void read(Input &in) override {
            // Getting length of model (it is a material parameter too! - used in penalization)
            Scalar xyzmin, xyzmax;
            in.get("x_min", xyzmin);
            in.get("x_max", xyzmax);
            Length_x = xyzmax - xyzmin;
            in.get("y_min", xyzmin);
            in.get("y_max", xyzmax);
            Length_y = xyzmax - xyzmin;
            in.get("z_min", xyzmin);
            in.get("z_max", xyzmax);
            Length_z = xyzmax - xyzmin;

            in.get("a", a);
            in.get("b", b);
            in.get("d", d);
            in.get("f", f);
            in.get("regularization", regularization);
            in.get("pressure", pressure);
            in.get("use_pressure", use_pressure);

            in.get("use_penalty_irreversibility", use_penalty_irreversibility);
            in.get("penalty_tol", penalty_tol);
            in.get("penalty_tol_non_neg", penalty_tol_non_neg);

            in.get("use_crack_set_irreversibiblity", use_crack_set_irreversibiblity);
            in.get("crack_set_tol", crack_set_tol);

            in.get("l_0", l_0);
            in.get("pressure0", pressure0);

            in.get("turn_off_uc_coupling", turn_off_uc_coupling);
            in.get("turn_off_cu_coupling", turn_off_cu_coupling);

            // Checking derivatives
            in.get("check_elastic_energy", check_elastic_energy);
            in.get("check_fracture_energy", check_fracture_energy);

            in.get("mobility", mobility);
            in.get("use_mobility", use_mobility);

            // MATERIAL PARAMETERS =============================================
            in.get("length_scale", length_scale);
            in.get("mu", mu);
            in.get("lambda", lambda);
            in.get("nu", nu);
            in.get("E", E);
            in.get("fracture_toughness", fracture_toughness);
            in.get("tensile_strength", tensile_strength);

            // Initialising other parameters
            if (nu != 0.0 && E != 0.0) {
                initialise_Lame_parameters();
            } else {
                initialise_Young_Poisson_parameters();
            }
            // Must be done after lambda and mu
            fill_in_isotropic_elast_tensor();
            // END OF MATERIAL PARAMETERS =============================================

            std::string type;
            in.get("hetero_params", type);

            if (type == "HomogeneousBar") {
                Scalar nx;  // mesh resoultion in x
                in.get("nx", nx);

                Scalar xmin, xmax, ymin, ymax;
                in.get("x_min", xmin);
                in.get("x_max", xmax);
                in.get("y_min", ymin);
                in.get("y_max", ymax);

                hetero_params = [](const Point &, Scalar &, Scalar &, Scalar &, Scalar &, Scalar &, Scalar &) {};

            } else if (type == "SingleSedimentaryLayer") {
                number_of_layers_ = 1;

                Scalar bottom_layer_height_;
                Scalar top_layer_height_;
                Scalar interface_regularisation_length;
                bool include_interface_layer{false};

                in.get("bottom_layer_height", bottom_layer_height_);
                in.get("top_layer_height", top_layer_height_);
                in.get("interface_regularisation_length", interface_regularisation_length);
                bottom_layer_heights_.push_back(bottom_layer_height_);
                top_layer_heights_.push_back(top_layer_height_);

                Scalar E1, E2, nu1, nu2, Gc1, Gc2, Gc_int, l_, ft1{1}, ft2{1}, ft_int{1};
                in.get("E_1", E1);
                in.get("E_2", E2);
                in.get("nu_1", nu1);
                in.get("nu_2", nu2);
                in.get("Gc_1", Gc1);
                in.get("Gc_2", Gc2);
                l_ = length_scale;

                in.get("ft_1", ft1);
                in.get("ft_2", ft2);
                tensile_strength = ft1;

                in.get("include_interface_layer", include_interface_layer);
                if (include_interface_layer) {
                    in.get("Gc_int", Gc_int);
                    in.get("ft_int", ft_int);
                } else {  // interface layer same toughness as outer layer
                    Gc_int = Gc2;
                    ft_int = ft2;
                }

                E = E1;
                nu = nu1;
                fracture_toughness = Gc1;

                // Random variation
                bool random_variation{false};
                Scalar toughness_deviation{1.0};
                Scalar youngs_deviation{1.0};
                in.get("random_variation", random_variation);
                in.get("random_toughness_deviation", toughness_deviation);
                in.get("random_youngs_deviation", youngs_deviation);

                double use_random = random_variation ? 1.0 : 0.0;
                std::normal_distribution<double> distribution_G(0., toughness_deviation);
                std::normal_distribution<double> distribution_E(0., youngs_deviation);
                std::default_random_engine generator;

                Scalar xmin, xmax, ymin, ymax;
                in.get("x_min", xmin);
                in.get("x_max", xmax);
                in.get("y_min", ymin);
                in.get("y_max", ymax);

                hetero_params = [E1,
                                 E2,
                                 nu1,
                                 nu2,
                                 Gc1,
                                 Gc2,
                                 Gc_int,
                                 ft1,
                                 ft2,
                                 ft_int,
                                 l_,
                                 bottom_layer_height_,
                                 top_layer_height_,
                                 xmin,
                                 xmax,
                                 ymin,
                                 ymax,
                                 use_random,
                                 distribution_E,
                                 distribution_G,
                                 generator,
                                 this](const Point &p,
                                       Scalar &mu_out,
                                       Scalar &lambda_out,
                                       Scalar &fracture_toughness_out,
                                       Scalar &tensile_strength,
                                       Scalar &E,
                                       Scalar &nu) mutable {
                    if (p[1] < bottom_layer_height_ ||
                        p[1] > top_layer_height_) {  // Shale (stronger and more compliant)
                        E = E2;
                        nu = nu2;
                        lambda_out = initialise_lambda(E2, nu2);
                        mu_out = E2 / (2. * (1. + nu2));
                        if (p[1] < bottom_layer_height_ - 1.5 * l_ || p[1] > top_layer_height_ + 1.5 * l_) {
                            fracture_toughness_out = Gc2;
                            tensile_strength = ft2;
                        } else {
                            fracture_toughness_out = Gc_int;
                            tensile_strength = ft_int;
                        }
                    } else {  // Dolostone (weaker and stiffer)

                        generator.seed(2.0 * (1e6 * p[1] * p[1] + 1e6 * p[0] * p[0]));

                        double noise_E = distribution_E(generator);
                        distribution_E.reset();

                        E = E1 + use_random * noise_E;
                        nu = nu1;
                        lambda_out = initialise_lambda(E, nu1);
                        mu_out = E / (2. * (1. + nu1));

                        double noise_G = distribution_G(generator);
                        distribution_G.reset();

                        fracture_toughness_out = Gc1 + use_random * noise_G;
                        tensile_strength = ft1;
                    }
                };
            } else if (type == "DoubleSedimentaryLayer") {
                number_of_layers_ = 2;

                Scalar bottom_layer_height_;
                Scalar top_layer_height_;
                Scalar bottom_layer_height2_;
                Scalar top_layer_height2_;
                bool include_interface_layer{false};

                in.get("bottom_layer_height", bottom_layer_height_);
                in.get("top_layer_height", top_layer_height_);
                bottom_layer_heights_.push_back(bottom_layer_height_);
                top_layer_heights_.push_back(top_layer_height_);
                in.get("bottom_layer_height2", bottom_layer_height2_);
                in.get("top_layer_height2", top_layer_height2_);
                bottom_layer_heights_.push_back(bottom_layer_height2_);
                top_layer_heights_.push_back(top_layer_height2_);

                Scalar E1, E2, nu1, nu2, Gc1, Gc2, Gc_int, l_, ft1{1}, ft2{1}, ft_int{1};
                in.get("E_1", E1);
                in.get("E_2", E2);
                in.get("nu_1", nu1);
                in.get("nu_2", nu2);
                in.get("Gc_1", Gc1);
                in.get("Gc_2", Gc2);
                l_ = length_scale;

                in.get("ft_1", ft1);
                in.get("ft_2", ft2);
                tensile_strength = ft1;

                in.get("include_interface_layer", include_interface_layer);
                if (include_interface_layer) {
                    in.get("Gc_int", Gc_int);
                    in.get("ft_int", ft_int);
                } else {  // interface layer same toughness as outer layer
                    Gc_int = Gc2;
                    ft_int = ft2;
                }

                E = E1;
                nu = nu1;
                fracture_toughness = Gc1;

                // Random variation
                bool random_variation{false};
                Scalar toughness_deviation{1.0};
                Scalar youngs_deviation{1.0};
                in.get("random_variation", random_variation);
                in.get("random_toughness_deviation", toughness_deviation);
                in.get("random_youngs_deviation", youngs_deviation);

                double use_random = random_variation ? 1.0 : 0.0;
                std::normal_distribution<double> distribution_G(0., toughness_deviation);
                std::normal_distribution<double> distribution_E(0., youngs_deviation);
                std::default_random_engine generator;

                Scalar xmin, xmax, ymin, ymax;
                in.get("x_min", xmin);
                in.get("x_max", xmax);
                in.get("y_min", ymin);
                in.get("y_max", ymax);

                hetero_params = [E1,
                                 E2,
                                 nu1,
                                 nu2,
                                 Gc1,
                                 Gc2,
                                 Gc_int,
                                 ft1,
                                 ft2,
                                 ft_int,
                                 l_,
                                 bottom_layer_height_,
                                 top_layer_height_,
                                 bottom_layer_height2_,
                                 top_layer_height2_,
                                 xmin,
                                 xmax,
                                 ymin,
                                 ymax,
                                 use_random,
                                 distribution_E,
                                 distribution_G,
                                 generator,
                                 this](const Point &p,
                                       Scalar &mu_out,
                                       Scalar &lambda_out,
                                       Scalar &fracture_toughness_out,
                                       Scalar &tensile_strength,
                                       Scalar &E,
                                       Scalar &nu) mutable {
                    if ((p[1] <= top_layer_height_ && p[1] >= bottom_layer_height_) ||
                        (p[1] <= top_layer_height2_ && p[1] >= bottom_layer_height2_)) {
                        generator.seed(2.0 * (1e6 * p[1] * p[1] + 1e6 * p[0] * p[0]));

                        double noise_E = distribution_E(generator);
                        distribution_E.reset();

                        E = E1 + use_random * noise_E;
                        nu = nu1;
                        lambda_out = initialise_lambda(E, nu1);
                        mu_out = E / (2. * (1. + nu1));

                        double noise_G = distribution_G(generator);
                        distribution_G.reset();

                        fracture_toughness_out = Gc1 + use_random * noise_G;
                        tensile_strength = ft1;

                    } else {  // Shale (stronger and more compliant)
                        E = E2;
                        nu = nu2;
                        lambda_out = initialise_lambda(E2, nu2);
                        mu_out = E2 / (2. * (1. + nu2));
                        if (p[1] < bottom_layer_height_ - 1.5 * l_ || p[1] > top_layer_height_ + 1.5 * l_) {
                            fracture_toughness_out = Gc2;
                            tensile_strength = ft2;
                        } else {
                            fracture_toughness_out = Gc_int;
                            tensile_strength = ft_int;
                        }
                    }  // end of double sedimentary layer function
                };
            } else if (type == "MultiSedimentaryLayer") {
                in.get("number_of_layers", number_of_layers_);

                Scalar bottom_layer_height;
                Scalar average_layer_spacing;
                Scalar average_layer_height;
                Scalar spacing_variation, height_variation;
                bool include_interface_layer{false};

                in.get("bottom_layer_height", bottom_layer_height);
                in.get("average_layer_spacing", average_layer_spacing);
                in.get("average_layer_height", average_layer_height);
                in.get("layer_spacing_variation", spacing_variation);
                in.get("layer_height_variation", height_variation);

                bottom_layer_heights_.push_back(bottom_layer_height);
                top_layer_heights_.push_back(bottom_layer_height + average_layer_height);
                for (size_t i = 0; i < number_of_layers_ - 1; ++i) {
                    bottom_layer_heights_.push_back(top_layer_heights_[i] + average_layer_spacing +
                                                    std::pow(-1.0, double(i)) * spacing_variation);
                    top_layer_heights_.push_back(bottom_layer_heights_[i + 1] + average_layer_height +
                                                 std::pow(-1.0, double(i)) * height_variation);
                }

                Scalar E1, E2, nu1, nu2, Gc1, Gc2, Gc_int, l_, ft1{1}, ft2{1}, ft_int{1};
                in.get("E_1", E1);
                in.get("E_2", E2);
                in.get("nu_1", nu1);
                in.get("nu_2", nu2);
                in.get("Gc_1", Gc1);
                in.get("Gc_2", Gc2);
                l_ = length_scale;

                in.get("ft_1", ft1);
                in.get("ft_2", ft2);
                tensile_strength = ft1;

                in.get("include_interface_layer", include_interface_layer);
                if (include_interface_layer) {
                    in.get("Gc_int", Gc_int);
                    in.get("ft_int", ft_int);
                } else {  // interface layer same toughness as outer layer
                    Gc_int = Gc2;
                    ft_int = ft2;
                }

                E = E1;
                nu = nu1;
                fracture_toughness = Gc1;

                // Random variation
                bool random_variation{false};
                Scalar toughness_deviation{1.0};
                Scalar youngs_deviation{1.0};
                in.get("random_variation", random_variation);
                in.get("random_toughness_deviation", toughness_deviation);
                in.get("random_youngs_deviation", youngs_deviation);

                double use_random = random_variation ? 1.0 : 0.0;
                std::normal_distribution<double> distribution_G(0., toughness_deviation);
                std::normal_distribution<double> distribution_E(0., youngs_deviation);
                std::default_random_engine generator;

                Scalar xmin, xmax, ymin, ymax;
                in.get("x_min", xmin);
                in.get("x_max", xmax);
                in.get("y_min", ymin);
                in.get("y_max", ymax);

                hetero_params = [E1,
                                 E2,
                                 nu1,
                                 nu2,
                                 Gc1,
                                 Gc2,
                                 Gc_int,
                                 ft1,
                                 ft2,
                                 ft_int,
                                 l_,
                                 xmin,
                                 xmax,
                                 ymin,
                                 ymax,
                                 use_random,
                                 distribution_E,
                                 distribution_G,
                                 generator,
                                 this](const Point &p,
                                       Scalar &mu_out,
                                       Scalar &lambda_out,
                                       Scalar &fracture_toughness_out,
                                       Scalar &tensile_strength,
                                       Scalar &E,
                                       Scalar &nu) mutable {
                    bool found_in_layer = false;
                    for (size_t i = 0; i < this->number_of_layers_; ++i) {
                        if (p[1] <= this->top_layer_heights_[i] && p[1] >= this->bottom_layer_heights_[i]) {
                            found_in_layer = true;
                        }
                    }

                    if (found_in_layer) {
                        generator.seed(2.0 * (1e6 * p[1] * p[1] + 1e6 * p[0] * p[0]));

                        double noise_E = distribution_E(generator);
                        distribution_E.reset();

                        E = E1 + use_random * noise_E;
                        nu = nu1;
                        lambda_out = initialise_lambda(E, nu1);
                        mu_out = E / (2. * (1. + nu1));

                        double noise_G = distribution_G(generator);
                        distribution_G.reset();

                        fracture_toughness_out = Gc1 + use_random * noise_G;
                        tensile_strength = ft1;

                    } else {  // Shale (stronger and more compliant)
                        E = E2;
                        nu = nu2;
                        lambda_out = initialise_lambda(E2, nu2);
                        mu_out = E2 / (2. * (1. + nu2));
                        fracture_toughness_out = Gc2;
                        tensile_strength = ft2;
                    }  // end of double sedimentary layer function
                };
            } else if (type == "RegularisedSingleLayer") {
                number_of_layers_ = 1;

                Scalar bottom_layer_height_;
                Scalar top_layer_height_;
                Scalar interface_regularisation_length;
                bool include_interface_layer{false};

                in.get("bottom_layer_height", bottom_layer_height_);
                in.get("top_layer_height", top_layer_height_);
                in.get("interface_regularisation_length", interface_regularisation_length);
                bottom_layer_heights_.push_back(bottom_layer_height_);
                top_layer_heights_.push_back(top_layer_height_);

                Scalar E1, E2, nu1, nu2, Gc1, Gc2, Gc_int, l_, ft1, ft2, ft_int;
                in.get("E_1", E1);
                in.get("E_2", E2);
                in.get("nu_1", nu1);
                in.get("nu_2", nu2);
                in.get("Gc_1", Gc1);
                in.get("Gc_2", Gc2);
                l_ = length_scale;

                in.get("ft_1", ft1);
                in.get("ft_2", ft2);
                tensile_strength = ft1;

                E = E1;
                nu = nu1;
                fracture_toughness = Gc1;

                // Random variation
                bool random_variation{false};
                Scalar toughness_deviation{1.0};
                in.get("random_variation", random_variation);
                in.get("random_standard_deviation", toughness_deviation);

                double use_random = random_variation ? 1.0 : 0.0;
                std::normal_distribution<double> distribution(0., toughness_deviation);
                std::default_random_engine generator;

                Scalar xmin, xmax, ymin, ymax;
                in.get("x_min", xmin);
                in.get("x_max", xmax);
                in.get("y_min", ymin);
                in.get("y_max", ymax);

                hetero_params = [E1,
                                 E2,
                                 nu1,
                                 nu2,
                                 Gc1,
                                 Gc2,
                                 ft1,
                                 ft2,
                                 l_,
                                 bottom_layer_height_,
                                 top_layer_height_,
                                 interface_regularisation_length,
                                 xmin,
                                 xmax,
                                 ymin,
                                 ymax,
                                 use_random,
                                 distribution,
                                 generator,
                                 this](const Point &p,
                                       Scalar &mu_out,
                                       Scalar &lambda_out,
                                       Scalar &fracture_toughness_out,
                                       Scalar &tensile_strength,
                                       Scalar &E,
                                       Scalar &nu) mutable {
                    if (p[1] <= top_layer_height_ - interface_regularisation_length / 2.0 &&
                        p[1] >= bottom_layer_height_ + interface_regularisation_length / 2.0) {  // stiffer layer

                        E = E1;
                        nu = nu1;
                        lambda_out = initialise_lambda(E1, nu1);
                        mu_out = E1 / (2. * (1. + nu1));

                        generator.seed((1e6 * p[1] * p[1] + 1e6 * p[0] * p[0]));
                        double noise = distribution(generator);
                        distribution.reset();
                        fracture_toughness_out = Gc1 + use_random * noise;
                        tensile_strength = ft1;
                    } else if ((p[1] >= top_layer_height_ - interface_regularisation_length / 2.0 &&
                                p[1] <= top_layer_height_) ||
                               (p[1] <= bottom_layer_height_ + interface_regularisation_length / 2.0 &&
                                p[1] >= bottom_layer_height_)) {
                        E = E1;
                        nu = nu1;
                        lambda_out = initialise_lambda(E1, nu1);
                        mu_out = E1 / (2. * (1. + nu1));

                        double dist_to_interf =
                            std::min(std::fabs(p[1] - (bottom_layer_height_ - interface_regularisation_length / 2.0)),
                                     std::fabs(p[1] - (top_layer_height_ + interface_regularisation_length / 2.0)));
                        double Gc_mixed = Gc2 * (1.0 - dist_to_interf / interface_regularisation_length) +
                                          Gc1 * dist_to_interf / interface_regularisation_length;
                        fracture_toughness_out = Gc_mixed;
                        tensile_strength = ft2;
                    } else if ((p[1] <= top_layer_height_ + interface_regularisation_length / 2.0 &&
                                p[1] >= top_layer_height_) ||
                               (p[1] >= bottom_layer_height_ - interface_regularisation_length / 2.0 &&
                                p[1] <= bottom_layer_height_)) {  // Shale (stronger and more compliant)

                        E = E2;
                        nu = nu2;
                        lambda_out = initialise_lambda(E2, nu2);
                        mu_out = E2 / (2. * (1. + nu2));
                        double dist_to_interf =
                            std::min(std::fabs(p[1] - (bottom_layer_height_ - interface_regularisation_length / 2.0)),
                                     std::fabs(p[1] - (top_layer_height_ + interface_regularisation_length / 2.0)));
                        double Gc_mixed = Gc2 * (1.0 - dist_to_interf / interface_regularisation_length) +
                                          Gc1 * dist_to_interf / interface_regularisation_length;
                        fracture_toughness_out = Gc_mixed;
                        tensile_strength = ft2;
                    } else {
                        E = E2;
                        nu = nu2;
                        lambda_out = initialise_lambda(E2, nu2);
                        mu_out = E2 / (2. * (1. + nu2));
                        fracture_toughness_out = Gc2;
                        tensile_strength = ft2;
                    }
                };
            } else if (type == "DoubleLayer") {
                number_of_layers_ = 2;

                Scalar bottom_layer_height_, bottom_layer_height2_;
                Scalar top_layer_height_, top_layer_height2_;
                bool include_interface_layer{false};

                in.get("bottom_layer_height", bottom_layer_height_);
                in.get("bottom_layer_height2", bottom_layer_height2_);
                in.get("top_layer_height", top_layer_height_);
                in.get("top_layer_height2", top_layer_height2_);
                bottom_layer_heights_.push_back(bottom_layer_height_);
                top_layer_heights_.push_back(top_layer_height_);
                bottom_layer_heights_.push_back(bottom_layer_height2_);
                top_layer_heights_.push_back(top_layer_height2_);

                Scalar E1, E2, nu1, nu2, Gc1, Gc2, Gc_int, l_, ft1, ft2, ft_int;
                in.get("E_1", E1);
                in.get("E_2", E2);
                in.get("nu_1", nu1);
                in.get("nu_2", nu2);
                in.get("Gc_1", Gc1);
                in.get("Gc_2", Gc2);
                l_ = length_scale;

                in.get("ft_1", ft1);
                in.get("ft_2", ft2);
                in.get("ft_int", ft_int);
                tensile_strength = ft1;

                in.get("include_interface_layer", include_interface_layer);
                if (include_interface_layer)
                    in.get("Gc_int", Gc_int);
                else  // interface layer same toughness as outer layer
                    Gc_int = Gc2;

                E = E1;
                nu = nu1;
                fracture_toughness = Gc1;

                // Random variation
                bool random_variation{false};
                Scalar toughness_deviation{1.0};
                in.get("random_variation", random_variation);
                in.get("random_standard_deviation", toughness_deviation);

                double use_random = random_variation ? 1.0 : 0.0;
                std::normal_distribution<double> distribution(0., toughness_deviation);
                std::default_random_engine generator;

                bool boundary_protection(false);
                Scalar layer_width(0);
                Scalar xmin, xmax, ymin, ymax;
                in.get("boundary_protection", boundary_protection);
                in.get("layer_width", layer_width);
                in.get("x_min", xmin);
                in.get("x_max", xmax);
                in.get("y_min", ymin);
                in.get("y_max", ymax);

                hetero_params = [E1,
                                 E2,
                                 nu1,
                                 nu2,
                                 Gc1,
                                 Gc2,
                                 Gc_int,
                                 ft1,
                                 ft2,
                                 ft_int,
                                 l_,
                                 bottom_layer_height_,
                                 top_layer_height_,
                                 bottom_layer_height2_,
                                 top_layer_height2_,
                                 boundary_protection,
                                 xmin,
                                 xmax,
                                 ymin,
                                 ymax,
                                 layer_width,
                                 use_random,
                                 distribution,
                                 generator,
                                 this](const Point &p,
                                       Scalar &mu_out,
                                       Scalar &lambda_out,
                                       Scalar &fracture_toughness_out,
                                       Scalar &tensile_strength,
                                       Scalar &E,
                                       Scalar &nu) mutable {
                    if ((p[1] < top_layer_height_ && p[1] > bottom_layer_height_) ||
                        (p[1] < top_layer_height2_ &&
                         p[1] > bottom_layer_height2_)) {  // Dolostone (weaker and stiffer)
                        E = E1;
                        nu = nu1;
                        lambda_out = initialise_lambda(E1, nu1);
                        mu_out = E1 / (2. * (1. + nu1));

                        generator.seed((1e6 * p[1] * p[1] + 1e6 * p[0] * p[0]));
                        double noise = distribution(generator);
                        distribution.reset();
                        fracture_toughness_out = Gc1 + use_random * noise;
                        tensile_strength = ft1;
                    } else {  // Shale (stronger and more compliant)
                        E = E2;
                        nu = nu2;
                        lambda_out = initialise_lambda(E2, nu2);
                        mu_out = E2 / (2. * (1. + nu2));
                        if (p[1] < bottom_layer_height_ - 1.5 * l_ || p[1] > top_layer_height_ + 1.5 * l_) {
                            fracture_toughness_out = Gc2;
                            tensile_strength = ft2;
                        } else {
                            fracture_toughness_out = Gc_int;
                            tensile_strength = ft_int;
                        }
                    }
                };
            }  // end of double layer
            else if (type == "BilinearTest") {
                Scalar xmin, xmax, ymin, ymax;
                in.get("x_min", xmin);
                in.get("x_max", xmax);
                in.get("y_min", ymin);
                in.get("y_max", ymax);

                initialise_Lame_parameters();
                Scalar mu_base = this->mu;
                Scalar lambda_base = this->lambda;
                Scalar FractureTough_base = this->fracture_toughness;

                hetero_params = [xmin, xmax, ymin, ymax, mu_base, lambda_base, FractureTough_base](
                                    const Point &p,
                                    Scalar &mu,
                                    Scalar &lambda,
                                    Scalar &fracture_tough,
                                    Scalar &,
                                    Scalar &,
                                    Scalar &) {
                    mu = mu_base * (1.0 + p(0) / xmax + p(1) / ymax);
                    lambda = lambda_base * (1.0 + p(0) / xmax + p(1) / ymax);
                    fracture_tough = FractureTough_base * (1.0 + p(0) / xmax + p(1) / ymax);
                };

            }  // end of Bilinear test

            // Initialising other parameters
            if (nu != 0.0 && E != 0.0) {
                initialise_Lame_parameters();
            } else {
                initialise_Young_Poisson_parameters();
            }
            if (mpi_world_rank() == 0) {
                utopia::out() << "E: " << E << "  nu: " << nu << "  Gc: " << fracture_toughness << " mu: " << mu
                              << "  lambda: " << lambda << " f_t: " << tensile_strength << "\n";
            }

            // Must be done after lambda and mu
            fill_in_isotropic_elast_tensor();
        }  // end of read

        PFFracParameters()
            : a(1.0),
              b(1.0),
              d(1.0),
              f(1.0),
              length_scale(0.0),
              fracture_toughness(1e-3),
              mu(80.0),
              lambda(120.0),
              // mu(100.0),
              // lambda(100.0),
              nu(0.0),
              E(0.0),
              l_0(1.0),
              pressure0(0.0),
              tensile_strength(0.0),
              regularization(1e-10),
              pressure(0.0),
              penalty_param_irreversible(0.0),
              penalty_param_non_neg(0.0),
              crack_set_tol(0.93),
              penalty_tol(0.01),
              // mobility(1e-5)
              mobility(1e-6),
              Length_x(0),
              Length_y(0),
              Length_z(0),
              use_pressure(false),
              turn_off_uc_coupling(false),
              turn_off_cu_coupling(false)

        {
            initialise_kappa();
        }

        void update(const Point &p, bool update_elastic_tensor) {
            if (hetero_params) {
                hetero_params(p, mu, lambda, fracture_toughness, tensile_strength, E, nu);
                initialise_kappa();
                // check_consistent_material_params();
            }
            if (update_elastic_tensor) fill_in_isotropic_elast_tensor();
        }

        void initialise_Lame_parameters() {
            lambda = E * nu / ((1. + nu) * (1. - static_cast<double>(Dim - 1) * nu));
            mu = E / (2. * (1. + nu));
        }

        Scalar initialise_lambda(Scalar E, Scalar nu) {
            return E * nu / ((1. + nu) * (1. - static_cast<double>(Dim - 1) * nu));
        }

        void initialise_Young_Poisson_parameters() {
            if constexpr (Dim == 3) {
                E = mu * (3.0 * lambda + 2.0 * mu) / (lambda + mu);
                nu = lambda / (2.0 * (lambda + mu));
            } else {
                E = 4.0 * mu * (lambda + mu) / (lambda + 2.0 * mu);
                nu = lambda / (lambda + 2.0 * mu);
            }
        }

        void initialise_kappa() { kappa = lambda + (2.0 * mu / static_cast<double>(Dim)); }

        std::pair<double, double> return_Lame_parameters(double E, double nu) {
            std::pair<double, double> lames;
            lames.first = initialise_lambda(E, nu);
            lames.second = E / (2. * (1. + nu));
            return lames;
        }

        bool kroneckerDelta(const SizeType &i, const SizeType &j) { return (i == j) ? 1.0 : 0.0; }

        void fill_in_isotropic_elast_tensor() {
            for (SizeType i = 0; i < Dim; ++i) {
                for (SizeType j = 0; j < Dim; ++j) {
                    for (SizeType k = 0; k < Dim; ++k) {
                        for (SizeType l = 0; l < Dim; ++l) {
                            Scalar val = this->lambda * kroneckerDelta(i, j) * kroneckerDelta(k, l);
                            val += this->mu * (kroneckerDelta(i, k) * kroneckerDelta(j, l));
                            val += this->mu * (kroneckerDelta(i, l) * kroneckerDelta(j, k));
                            elast_tensor.set(i, j, k, l, val);
                        }
                    }
                }
            }

            I4sym.identity_sym_k();          // E.P Fixed i4sym to correct version for Kappa
            I4shear.identity_shearmod(Dim);  // E.P i4shear for mu identity contribution to elasticity tensor
            initialise_kappa();
        }

        void check_consistent_material_params() {
            // std::cout <<"checking material now"<< std::endl;
            if (Dim == 2) {
                double E__mu_l =
                    std::abs(this->E - (4.0 * this->mu * (this->lambda + this->mu) / (this->lambda + 2.0 * this->mu)));
                double nu__l_mu = std::abs(this->nu - (this->lambda / (this->lambda + 2.0 * this->mu)));
                double l__E_nu = std::abs(this->lambda - (this->E * this->nu / ((1.0 + this->nu) * (1.0 - this->nu))));
                double mu__E_nu = std::abs(this->mu - (this->E / (2.0 * (1.0 + this->nu))));
                double k__l_mu = std::abs(this->kappa - (this->lambda + this->mu));
                double k__E_nu = std::abs(this->kappa - this->E / (2.0 * (1 - this->nu)));
                double k__E_mu = std::abs(this->kappa - this->E * this->mu / (4.0 * this->mu - this->E));
                double k__l_nu = std::abs(this->kappa - this->lambda * (1 + this->nu) / (2.0 * this->nu));
                double k__mu_nu = std::abs(this->kappa - this->mu * (1 + this->nu) / (1.0 - this->nu));

                if (E__mu_l <= 10.0 * std::numeric_limits<double>::epsilon() &&
                    nu__l_mu <= 10.0 * std::numeric_limits<double>::epsilon() &&
                    l__E_nu <= 10.0 * std::numeric_limits<double>::epsilon() &&
                    mu__E_nu <= 10.0 * std::numeric_limits<double>::epsilon() &&
                    k__l_mu <= 10.0 * std::numeric_limits<double>::epsilon() &&
                    k__E_nu <= 10.0 * std::numeric_limits<double>::epsilon() &&
                    k__E_mu <= 10.0 * std::numeric_limits<double>::epsilon() &&
                    k__l_nu <= 10.0 * std::numeric_limits<double>::epsilon() &&
                    k__mu_nu <= 10.0 * std::numeric_limits<double>::epsilon()) {
                    // std::cout << "\nCheck Passed. Material parameters are consistent"<<std::endl;
                } else {
                    std::cout << "\nWARNING!!! Material parameters are not consistent!!\n" << std::endl;
                    utopia::out() << "E: " << this->E << "  nu: " << this->nu << " mu: " << this->mu
                                  << "  lambda: " << this->lambda << "  kappa " << this->kappa << "\n";
                    std::cout << E__mu_l << " " << nu__l_mu << " " << l__E_nu << " " << mu__E_nu << " " << k__l_mu
                              << " " << k__E_nu << " " << k__E_mu << " " << k__l_nu << " " << k__mu_nu << " "
                              << std::endl;
                    exit(1);
                }
            } else {
                std::cout << "checking material for dim = 3 " << std::endl;
            }
        }

        Scalar a, b, d, f, length_scale, fracture_toughness, mu, lambda, kappa, nu, E, l_0, pressure0, tensile_strength;
        Scalar regularization, pressure, penalty_param_irreversible, penalty_param_non_neg, crack_set_tol, penalty_tol,
            penalty_tol_non_neg, mobility;
        Scalar Length_x, Length_y, Length_z;

        int number_of_layers_{0};
        std::vector<double> bottom_layer_heights_, top_layer_heights_;

        // Scalar E1, E2, nu1, nu2, Gc1, Gc2; //for hobbs three layer model
        bool use_penalty_irreversibility{false}, use_crack_set_irreversibiblity{false}, use_pressure{false};
        bool turn_off_uc_coupling{false}, turn_off_cu_coupling{false};
        bool use_mobility{false};

        Tensor4th<Scalar, Dim, Dim, Dim, Dim> elast_tensor;
        Tensor4th<Scalar, Dim, Dim, Dim, Dim> I4sym, I4shear;

        HeteroParamsFunction hetero_params;

        bool check_elastic_energy{false}, check_fracture_energy{false};
    };

    template <class FunctionSpace, int Dim = FunctionSpace::Dim>
    class PhaseFieldFracBase : public ExtendedFunction<typename FunctionSpace::Matrix, typename FunctionSpace::Vector> {
    public:
        using Scalar = typename FunctionSpace::Scalar;
        using SizeType = typename FunctionSpace::SizeType;
        using Vector = typename FunctionSpace::Vector;
        using Matrix = typename FunctionSpace::Matrix;
        using Device = typename FunctionSpace::Device;
        using Point = typename FunctionSpace::Point;

        using USpace = typename FunctionSpace::template Subspace<Dim>;
        using CSpace = typename FunctionSpace::template Subspace<1>;

        using UElem = typename USpace::ViewDevice::Elem;
        using CElem = typename CSpace::ViewDevice::Elem;
        using MixedElem = typename FunctionSpace::ViewDevice::Elem;

        using PFFracParameters = utopia::PFFracParameters<FunctionSpace>;
        using HeteroParamsFunction = typename PFFracParameters::HeteroParamsFunction;

        // FIXME
        using Shape = typename FunctionSpace::Shape;
        using Quadrature = utopia::Quadrature<Shape, 2 * (Shape::Order)>;

        static const int C_NDofs = CSpace::NDofs;
        static const int U_NDofs = USpace::NDofs;

        static const int NQuadPoints = Quadrature::NPoints;

        // using ExtendedFunction<typename FunctionSpace::Matrix, typename
        // FunctionSpace::Vector>::get_eq_constrains_flg; using
        // ExtendedFunction<typename FunctionSpace::Matrix,
        //                        typename
        //                        FunctionSpace::Vector>::get_eq_constrains_values;

        // E.P FIX: HardCoded for AT2 Model
        // this computation follows eq. 50 from "On penalization in variational
        // phase-field models of britlle fracture, Gerasimov, Lorenzis"
        void configure_penalty_term_for_AT2() {
            assert(params_.use_penalty_irreversibility);
            Scalar tol2 = params_.penalty_tol * params_.penalty_tol;
            params_.penalty_param_irreversible = params_.fracture_toughness / params_.length_scale * (1.0 / tol2 - 1.0);
            // if (mpi_world_rank()==0)
            //     utopia::out() << "Lengthscale: " << params_.length_scale << "  Penalty: " <<
            //     params_.penalty_param_irreversible << std::endl;
        }

        void read(Input &in) override {
            // reading parameters
            params_.read(in);

            in.get("checkpoint_each", checkpoint_each_);
            in.get("checkpoint_path", checkpoint_path_);

            // Sets a length scale based on mesh size if it is not set in the Yaml file (NOT GOOD SINCE different
            // material properties for each level
            if (params_.length_scale == 0) {
                params_.length_scale = 2.0 * space_.mesh().min_spacing();
            }
            if (mpi_world_rank() == 0) {
                utopia::out() << "using ls = " << params_.length_scale << "  \n";
            }

            // configuring penalty term
            if (params_.use_penalty_irreversibility) configure_penalty_term_for_AT2();

            in.get("use_dense_hessian", use_dense_hessian_);
            in.get("check_derivatives", check_derivatives_);
            in.get("diff_controller", diff_ctrl_);

            init_force_field(in);
        }

        // this is a bit of hack
        void set_pressure(const Scalar &pressure) { params_.pressure = pressure; }

        void use_crack_set_irreversibiblity(const bool &flg) { params_.use_crack_set_irreversibiblity = flg; }
        void use_penalty_irreversibility(const bool &flg) { params_.use_penalty_irreversibility = flg; }

        void turn_off_uc_coupling(const bool &flg) { params_.turn_off_uc_coupling = flg; }
        void turn_off_cu_coupling(const bool &flg) { params_.turn_off_cu_coupling = flg; }

        PFFracParameters &non_const_params() const { return const_cast<PhaseFieldFracBase *>(this)->params_; }

        int number_of_layers() { return params_.number_of_layers_; }
        double bottom_layer_height(int i) { return params_.bottom_layer_heights_[i]; }
        double top_layer_height(int i) { return params_.top_layer_heights_[i]; }

        void init_force_field(Input &in) {
            UTOPIA_TRACE_SCOPE("PhaseFieldFracBase::init_force_field");

            in.get("neumann_bc", [&](Input &in) {
                in.get_all([&](Input &node) {
                    if (empty(force_field_)) {
                        space_.create_vector(force_field_);
                        force_field_.set(0.0);
                    }

                    auto bc = utopia::make_unique<NeumannBoundaryCondition<FunctionSpace>>(space_);
                    bc->read(node);
                    bc->apply(force_field_);
                    neumann_bcs.push_back(std::move(bc));
                });
            });

            // if (false)
            // if (true) {
            //     space_.write("force_field.vtr", force_field_);
            // }
        }

        PhaseFieldFracBase(FunctionSpace &space) : space_(space) {
            // in case of constant pressure field
            // if(params_.pressure){
            params_.use_pressure = true;
            setup_constant_pressure_field(params_.pressure);
            // }

            // needed for ML setup
            space_.create_vector(this->_x_eq_values);
            space_.create_vector(this->_eq_constrains_flg);

            this->local_x_ = std::make_shared<Vector>();
            space_.create_local_vector(*this->local_x_);

            this->local_pressure_field_ = std::make_shared<Vector>();
            space_.create_local_vector(*this->local_pressure_field_);

            this->local_c_old_ = std::make_shared<Vector>();
            space_.create_local_vector(*this->local_c_old_);
        }

        void set_hetero_params(HeteroParamsFunction hetero_params) { params_.hetero_params = hetero_params; }

        PhaseFieldFracBase(FunctionSpace &space, const PFFracParameters &params) : space_(space), params_(params) {
            this->local_x_ = std::make_shared<Vector>();
            space_.create_local_vector(*this->local_x_);

            this->local_pressure_field_ = std::make_shared<Vector>();
            space_.create_local_vector(*this->local_pressure_field_);

            this->local_c_old_ = std::make_shared<Vector>();
            space_.create_local_vector(*this->local_c_old_);
        }

        void use_dense_hessian(const bool val) { use_dense_hessian_ = val; }

        inline bool initialize_hessian(Matrix &H, Matrix & /*H_pre*/) const override {
            // if(use_dense_hessian_) {
            //     H = local_zeros({space_.n_dofs(), space_.n_dofs()}); //FIXME
            // } else {
            space_.create_matrix(H);
            // }
            return true;
        }

        bool update(const Vector &x) override { return true; }

        ////////////////////////////////////////////////////////////////////////////////////
        virtual bool fracture_energy(const Vector & /*x_const*/, Scalar & /*val*/) const = 0;
        virtual bool elastic_energy(const Vector & /*x_const*/, Scalar & /*val*/) const = 0;

        // elastic energy specified by inherited class
        virtual bool elastic_energy_in_middle_layer(const Vector & /*x_const*/, Scalar & /*val*/, int) const {
            std::cout << "PhaseFieldFracBase::elastic_energy_in_middle_layer(), Please define method in "
                         "derived class!"
                      << std::endl;
            exit(1);
        }

        // elastic energy specified by inherited class
        virtual bool fracture_energy_in_middle_layer(const Vector & /*x_const*/, Scalar & /*val*/, int) const {
            std::cout << "PhaseFieldFracBase::fracture_energy_in_middle_layer(), Please define method in "
                         "derived class!"
                      << std::endl;
            exit(1);
        }

        // compute total crack volume (TCV), so we can compare to exact solution
        // !!! E.P Changed: USED TO RETURN error INSTEAD OF THE TOTAL CRACK VOLUME
        virtual bool compute_tcv(const Vector &x_const, Scalar &computed_tcv) const {
            UTOPIA_TRACE_REGION_BEGIN("PFBase::compute_tcv");
            const Scalar PI = 3.141592653589793238463;

            Scalar tcv_exact = 0.0;
            computed_tcv = 0.0;

            USpace U;
            this->space_.subspace(1, U);
            CSpace C = this->space_.subspace(0);

            CSpace U1 = this->space_.subspace(1);
            CSpace U2 = this->space_.subspace(2);
            CSpace U3 = this->space_.subspace(3);

            ///////////////////////////////////////////////////////////////////////////

            // update local vector x
            this->space_.global_to_local(x_const, *this->local_x_);
            auto u_coeff = std::make_shared<Coefficient<USpace>>(U, this->local_x_);
            auto c_coeff = std::make_shared<Coefficient<CSpace>>(C, this->local_x_);

            auto u_coeff1 = std::make_shared<Coefficient<CSpace>>(U1, this->local_x_);
            auto u_coeff2 = std::make_shared<Coefficient<CSpace>>(U2, this->local_x_);
            auto u_coeff3 = std::make_shared<Coefficient<CSpace>>(U3, this->local_x_);

            FEFunction<CSpace> c_fun(c_coeff);
            FEFunction<USpace> u_fun(u_coeff);

            FEFunction<CSpace> u1_fun(u_coeff1);
            FEFunction<CSpace> u2_fun(u_coeff2);
            FEFunction<CSpace> u3_fun(u_coeff3);

            ////////////////////////////////////////////////////////////////////////////

            Quadrature q;

            auto c_val = c_fun.value(q);
            auto c_grad = c_fun.gradient(q);
            auto u_val = u_fun.value(q);

            auto u1_val = u1_fun.value(q);
            auto u2_val = u2_fun.value(q);
            auto u3_val = u3_fun.value(q);

            auto differential = C.differential(q);

            {
                auto U_view = U.view_device();
                auto C_view = C.view_device();
                auto U1_view = U1.view_device();
                auto U2_view = U2.view_device();
                auto U3_view = U3.view_device();

                auto c_view = c_val.view_device();
                auto c_grad_view = c_grad.view_device();
                auto u_view = u_val.view_device();

                auto u1_view = u1_val.view_device();
                auto u2_view = u2_val.view_device();
                auto u3_view = u3_val.view_device();

                auto differential_view = differential.view_device();

                Device::parallel_reduce(
                    this->space_.element_range(),
                    UTOPIA_LAMBDA(const SizeType &i) {
                        StaticVector<Scalar, NQuadPoints> c;
                        StaticVector<Scalar, NQuadPoints> u1;
                        StaticVector<Scalar, NQuadPoints> u2;
                        StaticVector<Scalar, NQuadPoints> u3;

                        CElem c_e;
                        C_view.elem(i, c_e);
                        c_view.get(c_e, c);

                        CElem u1_e;
                        U1_view.elem(i, u1_e);
                        u1_view.get(u1_e, u1);

                        CElem u2_e;
                        U2_view.elem(i, u2_e);
                        u2_view.get(u2_e, u2);

                        CElem u3_e;
                        U3_view.elem(i, u3_e);
                        u3_view.get(u3_e, u3);

                        auto c_grad_el = c_grad_view.make(c_e);

                        auto dx = differential_view.make(c_e);

                        Scalar tcv_mine = 0.0;

                        // TCV = \int_\Omega u \cdot \nabla \phi
                        for (SizeType qp = 0; qp < NQuadPoints; ++qp) {
                            if (c[qp] > 0.0) {  // E.P Zero gives correct answer! Added to catch only aperture of fully
                                                // cracked fracture
                                if (Dim == 2) {
                                    tcv_mine -= ((u1[qp] * c_grad_el[qp](0)) + (u2[qp] * c_grad_el[qp](1))) * dx(qp);
                                } else {
                                    tcv_mine -= ((u1[qp] * c_grad_el[qp](0)) + (u2[qp] * c_grad_el[qp](1)) +
                                                 (u3[qp] * c_grad_el[qp](2))) *
                                                dx(qp);
                                }
                            }
                        }

                        assert(tcv_mine == tcv_mine);
                        return tcv_mine;
                    },
                    computed_tcv);
            }

            computed_tcv = x_const.comm().sum(computed_tcv);
            assert(computed_tcv == computed_tcv);

            // params depends on initial fracture lenght and pressure
            if (Dim == 2) {
                tcv_exact = 2.0 * this->params_.pressure0 * this->params_.l_0 * this->params_.l_0 *
                            (1.0 - this->params_.nu * this->params_.nu) * PI / this->params_.E;
            } else {
                tcv_exact = 16.0 * this->params_.pressure0 * this->params_.l_0 * this->params_.l_0 * this->params_.l_0 *
                            (1.0 - this->params_.nu * this->params_.nu) / this->params_.E / 3.0;
            }

            // error = device::abs(computed_tcv - tcv_exact);
            //            if (mpi_world_rank() == 0) {
            //                std::cout << "computed_tcv: " << computed_tcv << "  exact: " << tcv_exact << "  error: "
            //                << error
            //                          << "\n ";
            //            }

            UTOPIA_TRACE_REGION_END("PFBase::compute_tcv");
            return true;
        }

        virtual void compute_cod(const Vector &x_const, Scalar &error) const {
            const Scalar PI = 3.141592653589793238463;

            // coordinates of the point at which we compute crack openning
            Scalar x_coord = 0.0;
            Scalar y_coord = 0.0;

            Scalar cod_exact = 4.0 / PI * this->params_.pressure0 * this->params_.l_0 *
                               (1. - (this->params_.nu * this->params_.nu)) / this->params_.E;

            Scalar rho = std::sqrt((x_coord * x_coord) + (y_coord * y_coord));
            Scalar rho_l0 = (rho / this->params_.l_0);

            cod_exact = cod_exact * std::sqrt(1. - rho_l0 * rho_l0);

            // search for max disp_y (should be symmetric)
            Scalar cod_computed = 0.0;

            {
                Read<Vector> r(x_const);

                Range range_w = range(x_const);
                for (SizeType i = range_w.begin(); i != range_w.end(); i++) {
                    if (i % (Dim + 1) == 2 && x_const.get(i) > cod_computed) {
                        cod_computed = x_const.get(i);
                    }
                }
            }

            cod_computed = x_const.comm().max(cod_computed);

            error = device::abs(cod_computed - cod_exact);

            if (mpi_world_rank() == 0) {
                std::cout << "cod_exact: " << cod_exact << "  cod_computed: " << cod_computed << "  error: " << error
                          << "\n ";
            }
        }

        bool export_material_params(std::string output_path) {
            UTOPIA_TRACE_REGION_BEGIN("PhaseFieldFracBase::export_mechanical_params");

            static const int total_components = 6.0;  // E, nu, Gc, lambda, mu, tensile_strength

            using PSpace = typename FunctionSpace::template Subspace<total_components>;
            using SElem = typename PSpace::ViewDevice::Elem;
            using WSpace = typename FunctionSpace::template Subspace<1>;

            Vector w;
            Vector g;

            /// Creating strain subspace
            // cloning mesh
            auto param_mesh = this->space_.mesh().clone(total_components);
            assert(param_mesh->n_components() == total_components);
            // Creating Subspace with cloned mesh

            PSpace S(std::move(param_mesh));
            WSpace C(this->space_.mesh().clone(1));

            S.create_vector(g);
            C.create_vector(w);
            ///////////////////////////////////////////////////////////////////////////

            {
                ////////////////////////////////////////////////////////////////////////////

                auto S_view = S.view_device();
                auto C_view = C.view_device();

                // Preparing the vector for which the parameter function space knows the dimensions (nodes*components),
                // so that we can write on this later
                auto g_view = S.assembly_view_device(g);
                auto w_view = C.assembly_view_device(w);

                Device::parallel_for(
                    this->space_.element_range(), UTOPIA_LAMBDA(const SizeType &i) {
                        StaticVector<Scalar, total_components * C_NDofs> material_params;
                        StaticVector<Scalar, C_NDofs> node_count;
                        material_params.set(.0);
                        node_count.set(1.0);

                        SElem s_e;
                        S_view.elem(i, s_e);  // just needed for add_vector into g (and node coords)

                        CElem c_e;
                        C_view.elem(i, c_e);  // getting element for storing wieghts in CSpace

                        for (SizeType n = 0; n < C_NDofs; n++) {
                            ////////////////////////////////////////////
                            bool update_elast_tensor = false;
                            Point coord;
                            s_e.node(n, coord);
                            non_const_params().update(coord, update_elast_tensor);
                            ////////////////////////////////////////////

                            Scalar Gc = this->params_.fracture_toughness;
                            Scalar mu = this->params_.mu;
                            Scalar l = this->params_.lambda;
                            Scalar E = mu * (3. * l + 2. * mu) / (l + mu);
                            Scalar nu = E / (2.0 * mu) - 1.;
                            Scalar tens_strength =
                                std::pow(3. / 8. * Gc * E / this->params_.length_scale, 0.5);  // AT1 hardcoded

                            material_params[n] = E;
                            material_params[C_NDofs + n] = nu;
                            material_params[C_NDofs * 2 + n] = Gc;
                            material_params[C_NDofs * 3 + n] = l;
                            material_params[C_NDofs * 4 + n] = mu;
                            material_params[C_NDofs * 5 + n] = tens_strength;
                        }

                        C_view.add_vector(c_e, node_count, w_view);
                        S_view.add_vector(s_e, material_params, g_view);
                    });  // end of parallel loop
            }

            {
                auto mat_params = local_view_device(g);
                auto node_count = local_view_device(w);
                auto r = local_range_device(w);  // range of vector w (using primitivo di utopio)
                parallel_for(
                    r, UTOPIA_LAMBDA(int i) {
                        auto wi = node_count.get(i);  // extracts vector component
                        for (int k = 0; k < total_components; k++) {
                            int nodal_offset =
                                i * total_components;  // vector g is 2*straincomponents bigger than vector of weights w
                            auto si = mat_params.get(nodal_offset +
                                                     k);  // get k'th strain corresponding to node i with weight i
                            mat_params.set(nodal_offset + k, si / wi);  // normalise the strain value by the weight wi
                        }
                    });
            }  // incase backed PETSC needs synchronisation (create view in scopes and destroy them when not needed)

            rename("E nu Gc lambda mu st", g);
            output_path += "_params.vtr";
            S.write(output_path, g);  // Function space knows how to write

            UTOPIA_TRACE_REGION_END("PhaseFieldFracBase::export_mechanical_params");
            return true;
        }

        virtual bool export_strain_and_stress(std::string output_path, const Vector &x_const, const Scalar time) const {
            UTOPIA_TRACE_REGION_BEGIN("PhaseFieldFracBase::strain");

            static const int strain_components = (Dim - 1) * 3;
            static const int Total_components = strain_components * 2;

            using WSpace = typename FunctionSpace::template Subspace<1>;
            using SSpace = typename FunctionSpace::template Subspace<Total_components>;
            using SElem = typename SSpace::ViewDevice::Elem;

            Vector w;
            Vector g;
            // Getting displacement subspace
            USpace U;
            this->space_.subspace(1, U);

            WSpace C(this->space_.mesh().clone(1));
            CSpace CC = this->space_.subspace(0);

            /// Creating strain subspace

            // cloning mesh
            auto strain_mesh = this->space_.mesh().clone(Total_components);
            assert(strain_mesh->n_components() == Total_components);
            // Creating Subspace with cloned mesh

            SSpace S(std::move(strain_mesh));

            assert(S.n_dofs() == C.n_dofs() * Total_components);

            S.create_vector(g);
            C.create_vector(w);

            assert(g.size() == w.size() * Total_components);

            ///////////////////////////////////////////////////////////////////////////

            // update local vector x
            this->space_.global_to_local(x_const, *this->local_x_);  // Gets the vector local to the MPI processor
            auto u_coeff = std::make_shared<Coefficient<USpace>>(
                U, this->local_x_);  // Sets stage for getting accessing the element node variables
            auto c_coeff = std::make_shared<Coefficient<CSpace>>(CC, this->local_x_);

            // getting FEFunction Space which contains objects for shape function manipulation
            FEFunction<USpace> u_fun(u_coeff);
            FEFunction<CSpace> c_fun(c_coeff);

            {
                ////////////////////////////////////////////////////////////////////////////

                // Quadrature for shape function integration
                Quadrature q;

                // Creating objects for Nodal and Gradient interpolation
                auto u_val = u_fun.value(q);
                auto u_grad = u_fun.gradient(q);
                auto c_val = c_fun.value(q);

                // What is thAis for ???
                auto differential = C.differential(q);

                // auto v_grad_shape = U.shape_grad(q);
                auto c_shape = C.shape(q);            // Getting shape functions from FunctionSpace
                auto c_grad_shape = C.shape_grad(q);  // Getting derivative of shape functions from FunctionSpace

                CoefStrain<USpace, Quadrature> strain(u_coeff, q);  // displacement coefficients
                // Strain<USpace, Quadrature> ref_strain_u(U, q); //Test strains (just shape functions gradients for
                // strain)

                auto U_view = U.view_device();
                auto C_view = C.view_device();
                auto S_view = S.view_device();
                auto CC_view = CC.view_device();  // EP

                auto c_view = c_val.view_device();  // EP
                auto u_view = u_val.view_device();

                auto strain_view = strain.view_device();
                auto differential_view = differential.view_device();

                // auto v_grad_shape_view = v_grad_shape.view_device();
                auto c_shape_view = c_shape.view_device();  // scalar shape functions
                // auto c_grad_shape_view = c_grad_shape.view_device();

                // Preparing the vector for which the Strain function space nows the dimensions (nodes*components), so
                // that we can write on this later
                auto g_view = S.assembly_view_device(g);
                auto w_view = C.assembly_view_device(w);

                // auto ref_strain_u_view = ref_strain_u.view_device();

                Device::parallel_for(
                    this->space_.element_range(), UTOPIA_LAMBDA(const SizeType &i) {
                        StaticMatrix<Scalar, Dim, Dim> strain_value, stress_value;
                        StaticVector<Scalar, Total_components * C_NDofs> strain_and_stress_el_vec;
                        StaticVector<Scalar, C_NDofs> weight_el_vec;
                        StaticMatrix<Scalar, Dim, Dim> stress;  //, strain_p;

                        strain_and_stress_el_vec.set(0.0);
                        weight_el_vec.set(0.0);

                        ////////////////////////////////////////////

                        UElem u_e;
                        U_view.elem(i, u_e);
                        auto el_strain =
                            strain_view.make(u_e);  // el_strain.strain[qp] gives matrix of strain at int point

                        SElem s_e;
                        S_view.elem(i, s_e);  // just needed for add_vector into g

                        // auto u_grad_shape_el = v_grad_shape_view.make(u_e);
                        // auto &&u_strain_shape_el = ref_strain_u_view.make(u_e);

                        ////////////////////////////////////////////

                        CElem c_e;
                        C_view.elem(i, c_e);  // getting element for storing wieghts in CSpace

                        CElem cc_e;
                        CC_view.elem(i, cc_e);  // getting element for storing wieghts in CSpace

                        StaticVector<Scalar, NQuadPoints> c;
                        c_view.get(cc_e, c);

                        auto dx = differential_view.make(c_e);
                        auto c_shape_fun_el = c_shape_view.make(c_e);  // shape functions (scalar)

                        ////////////////////////////////////////////
                        bool update_elast_tensor = false;
                        Point centroid;
                        c_e.centroid(centroid);
                        non_const_params().update(centroid, update_elast_tensor);
                        ////////////////////////////////////////////

                        // loop over all nodes, and for each node, we integrate the strain at the int point weightwd by
                        // the distance to the node (shape function)
                        for (SizeType n = 0; n < C_NDofs; n++) {
                            strain_value.set(0.0);
                            stress_value.set(0.0);
                            for (SizeType qp = 0; qp < NQuadPoints; ++qp) {
                                auto shape = c_shape_fun_el(n, qp);  // shape function at N and Quadrature point
                                auto weight = dx(qp);                // no need for weights! we want length instead

                                // Calculate strain at quadrature point
                                auto &epsi = el_strain.strain[qp];

                                // calculate stress at quadrature
                                const Scalar tr_strain_u = trace(el_strain.strain[qp]);
                                compute_stress(this->params_,
                                               tr_strain_u,
                                               el_strain.strain[qp],
                                               stress);  // gets stress at quadrature point

                                strain_value +=
                                    epsi * shape * weight;  // matrix of strains added to existing nodal strain (
                                stress_value += quadratic_degradation(this->params_, c[qp]) * stress * shape *
                                                weight;  // Sum stress at integration point

                                // getting nodal weight for normalisation
                                weight_el_vec[n] += shape * weight;
                            }

                            // now we need to accumulate the matrix strain into engineering strain vector
                            int offset = C_NDofs, idx{0};
                            for (int r = 0; r < Dim; ++r) {
                                for (int c = r; c < Dim; c++) {
                                    strain_and_stress_el_vec[idx * offset + n] = stress_value(r, c);
                                    if (strain_components < Total_components)
                                        strain_and_stress_el_vec[(strain_components + idx) * offset + n] =
                                            strain_value(r, c);
                                    idx++;
                                }
                            }
                        }

                        // now adding element contribution to global strain and weight vector
                        S_view.add_vector(s_e, strain_and_stress_el_vec, g_view);
                        C_view.add_vector(c_e, weight_el_vec, w_view);
                    });  // end of parallel for

            }  // destruction of view activates MPI Synchronisation

            //            int weight_index = (i - (i % strain_components) ) / strain_components;

            {
                // disp(g.size());
                // disp(w.size());

                // viewing strain vector we just created
                auto strain_and_stress_view = local_view_device(g);
                auto weight_view = local_view_device(w);
                auto r = local_range_device(w);  // range of vector w (using primitivo di utopio)
                parallel_for(
                    r, UTOPIA_LAMBDA(int i) {
                        auto wi = weight_view.get(i);  // extracts vector component
                        for (int k = 0; k < Total_components; k++) {
                            int nodal_offset =
                                i * Total_components;  // vector g is 2*straincomponents bigger than vector of weights w
                            auto si = strain_and_stress_view.get(
                                nodal_offset + k);  // get k'th strain corresponding to node i with weight i
                            strain_and_stress_view.set(nodal_offset + k,
                                                       si / wi);  // normalise the strain value by the weight wi

                            //                            if (strain_components != total_components) {
                            //                                auto sig_i =
                            //                                    strain_and_stress_view.get(nodal_offset +
                            //                                    strain_components +
                            //                                                               k);  // get stress
                            //                                                               component which is offset
                            //                                                               additionally
                            //                                                                    // in the g vector by
                            //                                                                    the strain components
                            //                                strain_and_stress_view.set(nodal_offset +
                            //                                strain_components + k, sig_i / wi);
                            //                            }
                            // assert( std::signbit(si) == std::signbit(sig_i));
                        }
                    });
            }  // incase backed PETSC needs synchronisation (create view in scopes and destroy them when not needed)

            rename("stress and strain", g);
            output_path += "_strainstress_" + std::to_string(time) + ".vtr";
            S.write(output_path, g);  // Function space knows how to write

            UTOPIA_TRACE_REGION_END("PhaseFieldFracBase::strain");
            return true;
        }

        virtual void update_history_field(const Vector & /*x_const*/) const {}

        ////////////////////////////////////////////////////////////////////////////////////

        template <typename PhaseFieldValue, class StressShape, class Grad>
        UTOPIA_INLINE_FUNCTION static PhaseFieldValue bilinear_uu(const PFFracParameters &params,
                                                                  const PhaseFieldValue &phase_field_value,
                                                                  const StressShape &stress,
                                                                  const Grad &strain_test) {
            const auto gc = ((1.0 - params.regularization) * quadratic_degradation(params, phase_field_value) +
                             params.regularization);
            return inner(gc * stress, strain_test);
        }

        template <class Grad>
        UTOPIA_INLINE_FUNCTION static auto diffusion_c(const PFFracParameters &params,
                                                       const Grad &g_trial,
                                                       const Grad &g_test) {
            return params.fracture_toughness * params.length_scale * inner(g_trial, g_test);
        }

        template <typename FunValue>
        UTOPIA_INLINE_FUNCTION static FunValue reaction_c(const PFFracParameters &params,
                                                          const FunValue &trial,
                                                          const FunValue &test) {
            return (params.fracture_toughness / params.length_scale) * trial * test;
        }

        UTOPIA_INLINE_FUNCTION static Scalar elastic_deriv_cc(const PFFracParameters &params,
                                                              const Scalar &phase_field_value,
                                                              const Scalar &elastic_energy_positive,
                                                              const Scalar &trial,
                                                              const Scalar &test) {
            const Scalar dcc = (1.0 - params.regularization) * quadratic_degradation_deriv2(params, phase_field_value);
            return dcc * trial * elastic_energy_positive * test;
        }

        template <typename PhaseFieldValue, class Grad, typename TestFunction, class GradTest>
        UTOPIA_INLINE_FUNCTION static PhaseFieldValue grad_fracture_energy_wrt_c(
            const PFFracParameters &params,
            const PhaseFieldValue &phase_field_value,
            const Grad &phase_field_grad,
            const TestFunction &test_function,
            const GradTest &grad_test_function) {
            return params.fracture_toughness * ((1. / params.length_scale * phase_field_value * test_function) +
                                                (params.length_scale * inner(phase_field_grad, grad_test_function)));
        }

        template <class Stress, class FullStrain>
        UTOPIA_INLINE_FUNCTION static Scalar bilinear_uc(const PFFracParameters &params,
                                                         const Scalar &phase_field_value,
                                                         const Stress &stress_p,
                                                         const FullStrain &full_strain,
                                                         const Scalar &c_trial_fun) {
            return c_trial_fun * inner(quadratic_degradation_deriv(params, phase_field_value) * stress_p, full_strain);
        }

        template <typename PhaseFieldValue, class Grad>
        UTOPIA_INLINE_FUNCTION static PhaseFieldValue fracture_energy(const PFFracParameters &params,
                                                                      const PhaseFieldValue &phase_field_value,
                                                                      const Grad &phase_field_grad) {
            return params.fracture_toughness *
                   (1. / (2.0 * params.length_scale) * phase_field_value * phase_field_value +
                    params.length_scale / 2.0 * inner(phase_field_grad, phase_field_grad));
        }

        template <class GradShape>
        UTOPIA_INLINE_FUNCTION static Scalar bilinear_cc(const PFFracParameters &params,
                                                         const Scalar &phase_field_value,
                                                         const Scalar &elastic_energy_p,
                                                         const Scalar &shape_trial,
                                                         const Scalar &shape_test,
                                                         const GradShape &grad_trial,
                                                         const GradShape &grad_test) {
            return diffusion_c(params, grad_trial, grad_test) + reaction_c(params, shape_trial, shape_test) +
                   elastic_deriv_cc(params, phase_field_value, elastic_energy_p, shape_trial, shape_test);
        }

        template <typename C>
        UTOPIA_INLINE_FUNCTION static C quadratic_degradation(const PFFracParameters &, const C &c) {
            C imc = 1.0 - c;
            return imc * imc;
        }

        template <typename C>
        UTOPIA_INLINE_FUNCTION static C quadratic_degradation_deriv(const PFFracParameters &, const C &c) {
            C imc = 1.0 - c;
            return -2.0 * imc;
        }

        template <typename C>
        UTOPIA_INLINE_FUNCTION static C quadratic_degradation_deriv2(const PFFracParameters &, const C &) {
            return 2.0;
        }

        template <class Strain, class Stress>
        UTOPIA_INLINE_FUNCTION static void compute_stress(const PFFracParameters &params,
                                                          const Scalar &tr,
                                                          const Strain &strain,
                                                          Stress &stress) {
            stress = (2.0 * params.mu * strain) + (params.lambda * tr * (device::identity<Scalar>()));
        }

        Vector &old_solution() { return x_old_; }

        void get_old_solution(Vector &x) const { x = x_old_; }

        void set_old_solution(const Vector &x) {
            x_old_ = x;
            update_history_field(x_old_);
        }

        // E.P Question: What is this doing? Why is it hardcoded?
        void build_irreversility_constraint(Vector &lb) {
            {
                auto d_x_old = const_device_view(x_old_);

                auto lb_view = view_device(lb);
                parallel_for(
                    range_device(lb), UTOPIA_LAMBDA(const SizeType &i) {
                        if (i % (Dim + 1) == 0) {
                            lb_view.set(i, d_x_old.get(i));
                        } else {
                            lb_view.set(i, -9e15);
                        }
                    });
            }
        }

        // E.P Question: What is this doing? Why is it hardcoded?
        void build_irreversility_constraint(Vector &lb, Vector &ub) {
            {
                auto d_x_old = const_device_view(x_old_);

                auto lb_view = view_device(lb);
                auto ub_view = view_device(ub);
                parallel_for(
                    range_device(lb), UTOPIA_LAMBDA(const SizeType &i) {
                        if (i % (Dim + 1) == 0) {  // Is this for the phase field equation
                            lb_view.set(i, d_x_old.get(i));
                            ub_view.set(i, 1.0);
                        } else {
                            lb_view.set(i, -9e15);  // WHIS IS THIS HARD CODED??? Cant we just set nothing there...?
                            ub_view.set(i, 9e15);   // WHAT IS THIS FOR? Is this a constraint for the displacement?
                        }
                    });
            }
        }

        void make_iterate_feasible(const Vector &lb, const Vector &ub, Vector &x) {
            {
                auto d_x_old = view_device(x);

                auto lb_view = const_device_view(lb);
                auto ub_view = const_device_view(ub);
                parallel_for(
                    range_device(lb), UTOPIA_LAMBDA(const SizeType &i) {
                        if (i % (Dim + 1) == 0) {
                            // Scalar li = lb_view.get(i);
                            Scalar ui = ub_view.get(i);
                            auto xi = d_x_old.get(i);
                            // if (li >= xi) {
                            //     d_x_old.set(i, li);
                            // }
                            if (ui <= xi) {
                                d_x_old.set(i, ui);
                            }
                            // else {
                            //     // d_x_old.set(i, (ui <= xi) ? ui : xi);
                            // }
                        }
                    });
            }
        }

        // this 2 functions need to be moved to BC conditions
        void apply_zero_constraints_irreversibiblity(Vector &g) const {
            {
                auto d_x_old = const_device_view(x_old_);

                auto g_view = view_device(g);
                parallel_for(
                    range_device(g), UTOPIA_LAMBDA(const SizeType &i) {
                        if (i % (Dim + 1) == 0) {
                            if (d_x_old.get(i) > params_.crack_set_tol) {
                                g_view.set(i, 0.0);
                            }
                        }
                    });
            }
        }

        // this 2 functions need to be moved to BC conditions
        virtual void apply_zero_constraints_irreversibiblity(Vector &g, const Vector &x) const {
            {
                auto d_x_old = const_device_view(x_old_);
                auto d_x = const_device_view(x);

                auto g_view = view_device(g);
                parallel_for(
                    range_device(g), UTOPIA_LAMBDA(const SizeType &i) {
                        if (i % (Dim + 1) == 0) {
                            if (d_x_old.get(i) > params_.crack_set_tol || d_x.get(i) > params_.crack_set_tol) {
                                g_view.set(i, 0.0);
                            }
                        }
                    });
            }
        }

        // this 2 functions need to be moved to BC conditions
        void add_irr_values_markers(Vector &val, Vector &flg) const {
            {
                auto d_x_old = const_device_view(x_old_);

                auto val_view = view_device(val);
                parallel_for(
                    range_device(val), UTOPIA_LAMBDA(const SizeType &i) {
                        if (i % (Dim + 1) == 0) {
                            if (d_x_old.get(i) > params_.crack_set_tol) {
                                val_view.set(i, d_x_old.get(i));
                            }
                        }
                    });

                auto flg_view = view_device(flg);
                parallel_for(
                    range_device(flg), UTOPIA_LAMBDA(const SizeType &i) {
                        if (i % (Dim + 1) == 0) {
                            if (d_x_old.get(i) > params_.crack_set_tol) {
                                flg_view.set(i, 1.0);
                            }
                        }
                    });
            }
        }

        // this 2 functions need to be moved to BC conditions
        void add_irr_values_markers(Vector &val, Vector &flg, const Vector &x_current) const {
            {
                auto d_x_old = const_device_view(x_current);

                auto val_view = view_device(val);
                parallel_for(
                    range_device(val), UTOPIA_LAMBDA(const SizeType &i) {
                        if (i % (Dim + 1) == 0) {
                            if (d_x_old.get(i) >= params_.crack_set_tol) {
                                val_view.set(i, d_x_old.get(i));
                            }
                        }
                    });

                auto flg_view = view_device(flg);
                parallel_for(
                    range_device(flg), UTOPIA_LAMBDA(const SizeType &i) {
                        if (i % (Dim + 1) == 0) {
                            if (d_x_old.get(i) >= params_.crack_set_tol) {
                                flg_view.set(i, 1.0);
                            }
                        }
                    });
            }
        }

        void add_pf_constraints(const Vector &x_current) const {
            auto *p_this = const_cast<PhaseFieldFracBase<FunctionSpace> *>(this);

            Vector &bc_flgs = p_this->get_eq_constrains_flg();
            Vector &bc_values = p_this->get_eq_constrains_values();
            p_this->space_.apply_constraints(bc_values);
            p_this->space_.build_constraints_markers(bc_flgs);
            p_this->add_irr_values_markers(bc_values, bc_flgs, x_current);
        }

        // we should move this to BC conditions
        // also, will not run efficienetly in parallel
        virtual void apply_zero_constraints_irreversibiblity(Matrix &H) const {
            std::vector<SizeType> indices;
            {
                Read<Vector> r(x_old_);

                Range range_w = range(x_old_);
                for (SizeType i = range_w.begin(); i != range_w.end(); i++) {
                    if (i % (Dim + 1) == 0 && x_old_.get(i) > params_.crack_set_tol) {
                        indices.push_back(i);
                    }
                }
            }

            set_zero_rows(H, indices, 1.);
        }

        virtual void apply_zero_constraints_irreversibiblity(Matrix &H, const Vector &x) const {
            std::vector<SizeType> indices;
            {
                Read<Vector> r(x_old_);
                Read<Vector> r2(x);

                Range range_w = range(x_old_);
                for (SizeType i = range_w.begin(); i != range_w.end(); i++) {
                    if (i % (Dim + 1) == 0 &&
                        (x_old_.get(i) > params_.crack_set_tol || x.get(i) > params_.crack_set_tol)) {
                        indices.push_back(i);
                    }
                }
            }

            set_zero_rows(H, indices, 1.);
        }

        void pressure_field(const Vector &pressure_field) { pressure_field_ = pressure_field; }

        void setup_constant_pressure_field(const Scalar &p_val) {
            if (empty(pressure_field_)) {
                space_.create_vector(pressure_field_);
            }

            pressure_field_.set(p_val);
        }

        Vector &pressure_field() {
            if (empty(pressure_field_)) space_.create_vector(pressure_field_);

            return pressure_field_;
        }

        // void set_dt(const Scalar &dt) { dt_ = dt; }

        void set_time(const Scalar &t) {
            time_ = t;

            if (!neumann_bcs.empty()) {
                if (empty(force_field_)) {
                    space_.create_vector(force_field_);
                }

                force_field_.set(0.0);

                for (auto &&n : neumann_bcs) {
                    n->set_time(t);
                    n->apply(force_field_);
                }
            }
        }

        // Scalar get_dt() const { return dt_; }

        class RestartFiles {
        public:
            std::string last_file;
            int n_files;
        };

        inline static RestartFiles list_checkpoint_folder(const std::string &pattern) {
            glob_t gl;
            glob(pattern.c_str(), GLOB_MARK, NULL, &gl);

            int n_files = gl.gl_pathc;

            printf("n_files (%d):\n", n_files);
            for (int np = 0; np < n_files; np++) {
                printf("%s\n", gl.gl_pathv[np]);
            }

            RestartFiles ret;
            ret.n_files = n_files;

            if (n_files) {
                ret.last_file = gl.gl_pathv[n_files - 1];
            }

            globfree(&gl);
            return ret;
        }

        virtual bool restart_from_checkpoint(Vector &x, SizeType &time_step_counter, Scalar &time) {
            // RestartFiles rf = list_checkpoint_folder(checkpoint_path_ / ("_X_" + "*.dat"));
            // if (!rf.n_files) return false;
            // std::string time_mark = rf.last_file.substr(0, path.size());
            // time_mark = time_mark.substr(0, time_mark.size() - 4);
            // time = atof(time_mark.c_str());
            // time_step_counter = rf.n_files;
            // return space_.read_field(rf.last_file, x);

            Path check_file = checkpoint_path_ / "info.yaml";
            if (!check_file.exists()) return false;

            auto in = open_istream(check_file);
            if (!in) {
                return false;
            }

            Path last_checkpoint;
            in->require("last_checkpoint", last_checkpoint);
            in->require("time", time);
            in->require("step_count", time_step_counter);

            write_count_ = time_step_counter;
            return x.load(last_checkpoint);
        }

        virtual bool write_checkpoint(const Scalar time, const Vector &x) const {
            checkpoint_path_.make_dir();
            Path last_checkpoint =
                checkpoint_path_ / ("X_" + std::to_string(write_count_) + "_" + std::to_string(time) + ".dat");
            x.write(last_checkpoint);

            if (!x.comm().rank()) {
                // Only root process
                Path check_file = checkpoint_path_ / "info.yaml";
                std::ofstream os(check_file.c_str());
                os << "last_checkpoint: " << last_checkpoint << "\n";
                os << "time: " << time << "\n";
                os << "step_count: " << write_count_ << "\n";
                os.close();
            }

            return true;
        }

        virtual void write_to_file(const std::string &output_path, const Vector &x, const Scalar time) {
            if (!x.comm().rank()) {
                utopia::out() << "Saving: " << (output_path + "_X_" + std::to_string(time) + ".vtr") << "\n";
            }

            space_.write(output_path + "_X_" + std::to_string(time) + ".vtr", x);

            if (write_count_ % checkpoint_each_ == 0 && write_count_ != 0) {
                write_checkpoint(time, x);
            }

            write_count_++;
        }

        std::vector<double> WriteParametersToVector() {
            std::vector<double> v{
                params_.Length_x, params_.Length_y, params_.E, params_.length_scale, params_.fracture_toughness};
            return v;
        }

    protected:
        FunctionSpace &space_;
        PFFracParameters params_;
        DiffController<Matrix, Vector> diff_ctrl_;

        bool use_dense_hessian_{false};
        bool check_derivatives_{false};

        Vector x_old_;           // stores old solution  - used for treatment of
                                 // irreversibility constraint
        Vector pressure_field_;  // stores heterogenous pressure field - ideally, this
                                 // vector would have lower size than all 4 variables

        Vector force_field_;

        std::shared_ptr<Vector> local_x_;
        std::shared_ptr<Vector> local_pressure_field_;
        std::shared_ptr<Vector> local_c_old_;

        // Scalar dt_;
        Scalar time_{0};
        std::vector<std::unique_ptr<NeumannBoundaryCondition<FunctionSpace>>> neumann_bcs;

        ptrdiff_t checkpoint_each_{10};
        ptrdiff_t write_count_{0};
        Path checkpoint_path_{"checkpoint"};
    };

}  // namespace utopia

// clean-up macros
#undef UNROLL_FACTOR
#undef U_MIN
#endif
