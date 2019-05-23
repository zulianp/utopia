#ifndef UTOPIA_LAMEE_PARAMETERS_HPP
#define UTOPIA_LAMEE_PARAMETERS_HPP

#include "utopia_block_var.hpp"
#include <map>
#include <vector>
#include <utility>
#include <iostream>

namespace utopia {
    class LameeParameters : public Configurable {
    public:
        LameeParameters(const double default_mu = 1, const double default_lambda = 1)
        : default_mu(default_mu), default_lambda(default_lambda)
        {}

        void set_mu(const int id, const double mu)
        {
            mu_[id] = mu;
        }

        void set_lambda(const int id, const double lambda)
        {
            lambda_[id] = lambda;
        }

        double mu(const int id) const
        {
            auto it = mu_.find(id);
            if(it == mu_.end()) {
                return default_mu;
            }

            return it->second;
        }

        double lambda(const int id) const
        {
            auto it = lambda_.find(id);
            if(it == lambda_.end()) {
                return default_lambda;
            }

            return it->second;
        }

        static BlockVar<double> make_var_from_map(const double default_val, const std::map<int, double> &pairs)
        {
            std::vector<std::pair<int, double> > vec;
            vec.reserve(pairs.size());

            for(const auto &p : pairs) {
                vec.push_back(p);
            }

            return BlockVar<double>(default_val, vec);
        }

        inline BlockVar<double> var_lambda() const
        {
            return make_var_from_map(default_lambda, lambda_);
        }

        inline BlockVar<double> var_mu() const
        {
            return make_var_from_map(default_mu, mu_);
        }

        void describe(std::ostream &os) const {
            os << "mu: " << std::endl;
            os << default_mu << std::endl;

            for(auto mu : mu_) {
                os << " " << mu.first << " -> " << mu.second << "\n";
            }

            os << "lambda: " << std::endl;
            os << default_lambda << std::endl;

            for(auto lambda : lambda_) {
                os << " " << lambda.first << " -> " << lambda.second << "\n";
            }
        }

        class ShearModulus;
        class PoissonRatio;

        class FirstLameParameter {
        public:
            double value;

            FirstLameParameter() : value(-1) {}

            void init(const ShearModulus &mu, const PoissonRatio &ni)
            {
                value = 2.0 * (mu.value * ni.value)/(1 - 2 * ni.value);
            }

            inline bool valid() const
            {
                return value > 0.;
            }
        };

        class ShearModulus {
        public:
            double value;

            ShearModulus() : value(-1) {}

            inline bool valid() const
            {
                return value > 0.;
            }
        };

        class PoissonRatio {
        public:
            double value;

            PoissonRatio() : value(-1) {}

            inline bool valid() const
            {
                return value > 0.;
            }
        };

        static bool read_parameters(Input &in, ShearModulus &mu, FirstLameParameter &lambda)
        {
            in.get("mu", mu.value);
            in.get("lambda", lambda.value);

            if(!lambda.valid()) {
                PoissonRatio ni;
                in.get("poisson-ratio", ni.value);

                assert(ni.valid());
                lambda.init(mu, ni);
            }

            return mu.valid() && lambda.valid();
        }

        inline void read(Input &in) override
        {
            ShearModulus mu;
            FirstLameParameter lambda;

            if(read_parameters(in, mu, lambda)) {
                default_mu = mu.value;
                default_lambda = lambda.value;
            }

            in.get("sub-domains", [this, &mu, &lambda](Input &in) {
                in.get_all([this, &mu, &lambda](Input &array_in) {
                    int block = -1;

                    array_in.get("block", block);

                    if(block == -1) {
                        std::cerr << "[Error] specify block for each sub-domain" << std::endl;
                        return;
                    }

                    if(read_parameters(array_in, mu, lambda)) {
                        mu_[block]     = mu.value;
                        lambda_[block] = lambda.value;
                    } else {
                        mu_[block]     = default_mu;
                        lambda_[block] = default_lambda;
                    }
                });
            });
        }

        double default_mu, default_lambda;
        std::map<int, double> mu_;
        std::map<int, double> lambda_;
    };
}

#endif //UTOPIA_LAMEE_PARAMETERS_HPP
