#ifndef UTOPIA_FLOW_HPP
#define UTOPIA_FLOW_HPP

#include "utopia_FEForwardDeclarations.hpp"
#include "utopia_Integral.hpp"
#include "utopia_Model.hpp"
#include "utopia_UIForcingFunction.hpp"
#include "utopia_UIScalarSampler.hpp"

#include <iostream>

namespace utopia {

    template <typename Scalar>
    class UIScalarFunction final : public UIFunction<Scalar>, public Configurable {
    public:
        void read(Input &in) override {
            type_ = "";
            in.get("type", type_);

            if (type_.empty()) {
                clear();
                return;
            }

            if (type_ == "grid-function" || type_ == "ScalarSampler") {
                auto grid_sampler = std::make_shared<UIScalarSampler<Scalar>>();
                in.get("function", *grid_sampler);

                if (!grid_sampler->empty()) {
                    sampler_ = grid_sampler;
                }

            } else if (type_ == "subdomain-function" || type_ == "SubdomainFunction") {
                auto subdomain_fun = utopia::make_unique<UISubdomainFunction<Scalar>>();

                in.get("function", *subdomain_fun);

                if (!subdomain_fun->good()) {
                    sampler_ = std::make_shared<UIConstantFunction<Scalar>>(1.);
                } else {
                    if (!subdomain_fun->has_default()) {
                        subdomain_fun->set_default(utopia::make_unique<UIConstantFunction<Scalar>>(1.));
                    }

                    sampler_ = std::move(subdomain_fun);
                }
            } else if (type_ == "value") {
                Scalar value = 1.0;
                in.get("function", value);
                sampler_ = std::make_shared<UIConstantFunction<Scalar>>(value);
            }
        }

        UIScalarFunction() { clear(); }

        void clear() {
            type_ = "";
            sampler_ = std::make_shared<UIConstantFunction<Scalar>>(1.);
        }

        inline std::shared_ptr<UIFunction<Scalar>> sampler() { return sampler_; }

        bool set_current_block(const int subdomain_id) override { return sampler_->set_current_block(subdomain_id); }

        inline Scalar eval(const std::vector<Scalar> &x) const override { return sampler_->eval(x); }

        void describe(std::ostream &os) const override {
            if (sampler_) {
                os << type_ << '\n';
                sampler_->describe(os);
            }
        }

    private:
        std::string type_;
        std::shared_ptr<UIFunction<Scalar>> sampler_;
    };

    template <class FunctionSpace, class Matrix, class Vector>
    class Flow final : public Model<Matrix, Vector> {
    public:
        using Scalar = UTOPIA_SCALAR(Vector);
        typedef utopia::Traits<FunctionSpace> TraitsT;
        typedef typename TraitsT::Matrix ElementMatrix;

        Flow(FunctionSpace &space) : space_(space), forcing_function_(space), rescale_(1.0) {}

        inline bool assemble_hessian_and_gradient(const Vector &x, Matrix &hessian, Vector &gradient) override {
            Chrono c;
            c.start();

            auto u = trial(space_);
            auto v = test(space_);

            auto bilinear_form = inner(diffusion_tensor_ * grad(u), ctx_fun(permeability_.sampler()) * grad(v)) * dX;

            utopia::assemble(bilinear_form, hessian);
            forcing_function_.eval(x, gradient);

            // for newton methods
            // gradient -= hessian * x;

            if (rescale_ != 1.0) {
                hessian *= rescale_;
                gradient *= rescale_;
            }

            c.stop();
            std::cout << "Flow assemly time: " << c << std::endl;

            assemble_lower_dimensional_features(x, hessian, gradient);

            const Scalar sum_g = sum(gradient);
            std::cout << "sum_g: " << sum_g << std::endl;

            return true;
        }

        inline bool assemble_lower_dimensional_features(const Vector &x, Matrix &hessian, Vector &gradient) {
            if (lower_dimensional_tags_.empty()) {
                return true;
            }

            Chrono c;
            c.start();

            auto u = trial(space_);
            auto v = test(space_);

            const std::size_t n = lower_dimensional_tags_.size();

            Matrix trace_hessian;
            for (std::size_t i = 0; i < n; ++i) {
                auto side = lower_dimensional_tags_[i];

                auto bilinear_form = surface_integral(inner(
                                                          /* diffusion_tensor_ * */
                                                          grad(u),
                                                          ctx_fun(lower_dimensional_permeability_[i]) * grad(v)),
                                                      side);

                utopia::assemble(bilinear_form, trace_hessian);

                if (rescale_ != 1.0) {
                    trace_hessian *= rescale_;
                }

                hessian += trace_hessian;
            }

            c.stop();
            std::cout << "Flow assemly (sub-dimensional) time: " << c << std::endl;
            return false;
        }

        inline bool is_linear() const override { return true; }

        inline void clear() override {}

        inline void read(Input &in) override {
            read_permeability_tensor(in);
            in.get("permeability-function", permeability_);
            in.get("permeability_function", permeability_);
            in.get("coeff_function", permeability_);
            in.get("forcing-function", forcing_function_);
            in.get("forcing_functions", forcing_function_);

            lower_dimensional_tags_.clear();
            lower_dimensional_tags_.clear();

            std::cout << "lower-dimensional-permeability:\n";

            in.get("lower-dimensional-permeability", [this](Input &in) {
                in.get_all([this](Input &in) {
                    int tag = -1;

                    Scalar value = 1.0;
                    in.get("value", value);
                    in.get("side", tag);

                    std::cout << "side(" << tag << "): " << value << std::endl;

                    if (tag != -1) {
                        auto fun = std::make_shared<UIConstantFunction<Scalar>>(value);

                        lower_dimensional_permeability_.push_back(fun);
                        lower_dimensional_tags_.push_back(tag);
                    } else {
                        std::cerr << "[Error] malformed input for lower-dimensional-permeability!" << std::endl;
                    }
                });
            });
        }

        inline void rescale(const Scalar rescale) { rescale_ = rescale; }

    private:
        FunctionSpace &space_;
        UIForcingFunction<FunctionSpace, Vector> forcing_function_;

        UIScalarFunction<Scalar> permeability_;
        ElementMatrix diffusion_tensor_;

        std::vector<std::shared_ptr<UIFunction<Scalar>>> lower_dimensional_permeability_;
        std::vector<int> lower_dimensional_tags_;
        Scalar rescale_;

        void read_permeability_tensor(Input &in) {
            Scalar constant_permeability = 1.;
            in.get("permeability", constant_permeability);
            in.get("coeff", constant_permeability);

            Scalar permeabilities[3] = {constant_permeability, constant_permeability, constant_permeability};

            in.get("permeability-x", permeabilities[0]);
            in.get("permeability-y", permeabilities[1]);
            in.get("permeability-z", permeabilities[2]);

            const int dim = space_.mesh().spatial_dimension();

            diffusion_tensor_ = identity(dim, dim);

            {
                Write<ElementMatrix> w(diffusion_tensor_);
                for (int i = 0; i < dim; ++i) {
                    diffusion_tensor_.set(i, i, permeabilities[i]);
                }
            }

            std::cout << "global permeabilty tensor: ";

            for (int i = 0; i < dim; ++i) {
                std::cout << permeabilities[i] << " ";
            }

            std::cout << std::endl;
        }
    };

}  // namespace utopia

#endif  // UTOPIA_FLOW_HPP
