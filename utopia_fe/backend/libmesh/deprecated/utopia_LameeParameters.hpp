#ifndef UTOPIA_LAMEE_PARAMETERS_HPP
#define UTOPIA_LAMEE_PARAMETERS_HPP

#include <iostream>
#include <map>
#include <utility>
#include <vector>
#include "utopia_block_var.hpp"

#include "utopia_StressStrainParameters.hpp"

namespace utopia {
    class LameeParameters : public Configurable {
    public:
        LameeParameters(const double default_mu = 1, const double default_lambda = 1)
            : default_mu(default_mu), default_lambda(default_lambda) {}

        void set_mu(const int id, const double mu) { mu_[id] = mu; }

        void set_lambda(const int id, const double lambda) { lambda_[id] = lambda; }

        double mu(const int id) const {
            auto it = mu_.find(id);
            if (it == mu_.end()) {
                return default_mu;
            }

            return it->second;
        }

        double lambda(const int id) const {
            auto it = lambda_.find(id);
            if (it == lambda_.end()) {
                return default_lambda;
            }

            return it->second;
        }

        static BlockVar<double> make_var_from_map(const double default_val, const std::map<int, double> &pairs) {
            std::vector<std::pair<int, double> > vec;
            vec.reserve(pairs.size());

            for (const auto &p : pairs) {
                vec.push_back(p);
            }

            return BlockVar<double>(default_val, vec);
        }

        inline BlockVar<double> var_lambda() const { return make_var_from_map(default_lambda, lambda_); }

        inline BlockVar<double> var_mu() const { return make_var_from_map(default_mu, mu_); }

        void describe(std::ostream &os) const {
            os << "mu: " << std::endl;
            os << default_mu << std::endl;

            for (auto mu : mu_) {
                os << " " << mu.first << " -> " << mu.second << "\n";
            }

            os << "lambda: " << std::endl;
            os << default_lambda << std::endl;

            for (auto lambda : lambda_) {
                os << " " << lambda.first << " -> " << lambda.second << "\n";
            }
        }

        static bool read_parameters(Input &in, ShearModulus<double> &mu, FirstLameParameter<double> &lambda) {
            in.get("mu", mu.value);
            in.get("lambda", lambda.value);

            if (!lambda.valid()) {
                PoissonRatio<double> ni;
                in.get("poisson-ratio", ni.value);

                assert(ni.valid());
                lambda.init(mu, ni);
            }

            return mu.valid() && lambda.valid();
        }

        inline void read(Input &in) override {
            ShearModulus<double> mu;
            FirstLameParameter<double> lambda;

            if (read_parameters(in, mu, lambda)) {
                default_mu = mu.value;
                default_lambda = lambda.value;
            }

            in.get("sub-domains", [this, &mu, &lambda](Input &in) {
                in.get_all([this, &mu, &lambda](Input &array_in) {
                    int block = -1;

                    array_in.get("block", block);

                    if (block == -1) {
                        std::cerr << "[Error] specify block for each sub-domain" << std::endl;
                        return;
                    }

                    if (read_parameters(array_in, mu, lambda)) {
                        mu_[block] = mu.value;
                        lambda_[block] = lambda.value;
                    } else {
                        mu_[block] = default_mu;
                        lambda_[block] = default_lambda;
                    }
                });
            });
        }

        double default_mu, default_lambda;
        std::map<int, double> mu_;
        std::map<int, double> lambda_;
    };
}  // namespace utopia

#endif  // UTOPIA_LAMEE_PARAMETERS_HPP
