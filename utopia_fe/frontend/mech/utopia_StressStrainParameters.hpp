#ifndef UTOPIA_STRESS_STRAIN_PARAMETERS_HPP
#define UTOPIA_STRESS_STRAIN_PARAMETERS_HPP

#include "utopia_Input.hpp"

namespace utopia {

    template <typename T>
    class ShearModulus;

    template <typename T>
    class PoissonRatio;

    template <typename T>
    class YoungModulus;

    template <typename T>
    class StressStrainParameter : public Configurable {
    public:
        virtual ~StressStrainParameter() = default;

        T value;

        StressStrainParameter() : value(-1) {}

        inline const T &get() const {
            assert(valid());
            return value;
        }

        inline T &get() { return value; }

        inline bool valid() const { return value > 0.; }
    };

    template <typename T>
    class FirstLameParameter : public StressStrainParameter<T> {
    public:
        void read(Input &in) override {
            in.get("lambda", this->get());
            in.get("first_lame_parameter", this->get());
        }

        template <typename TSM, typename TPR>
        void init(const ShearModulus<TSM> &mu, const PoissonRatio<TPR> &ni) {
            this->get() = 2.0 * (mu.get() * ni.get()) / (1 - 2 * ni.get());
        }
    };

    template <typename T>
    class ShearModulus : public StressStrainParameter<T> {
    public:
        void read(Input &in) override {
            in.get("mu", this->get());
            in.get("shear_modulus", this->get());
        }

        template <typename TSM, typename TPR>
        void init(const YoungModulus<TSM> &E, const PoissonRatio<TPR> &ni) {
            this->get() = E.get() / (2 * (1 + ni.get()));
        }
    };

    template <typename T>
    class YoungModulus : public StressStrainParameter<T> {
    public:
        void read(Input &in) override {
            in.get("E", this->get());
            in.get("young_modulus", this->get());
        }
    };

    template <typename T>
    class PoissonRatio : public StressStrainParameter<T> {
    public:
        void read(Input &in) override {
            in.get("ni", this->get());
            in.get("poisson_ratio", this->get());
        }
    };

    template <typename FirstLameParameter_,
              typename ShearModulus_ = FirstLameParameter_,
              typename PoissonRatio_ = ShearModulus_,
              typename YoungModulus_ = PoissonRatio_>
    class StressStrainParameters : public Configurable {
    public:
        void read(Input &in) override {
            first_lame_parameter.read(in);
            shear_modulus.read(in);
            poisson_ratio.read(in);
            young_modulus.read(in);

            if (!shear_modulus.valid() && young_modulus.valid() && poisson_ratio.valid()) {
                shear_modulus.init(young_modulus, poisson_ratio);
            }

            if (!first_lame_parameter.valid() && poisson_ratio.valid() && shear_modulus.valid()) {
                first_lame_parameter.init(shear_modulus, poisson_ratio);
            }

            assert(first_lame_parameter.valid());
            assert(shear_modulus.valid());

            bool verbose = false;
            in.get("verbose", verbose);

            if (verbose) {
                utopia::out() << "first_lame_parameter:\t" << first_lame_parameter.get() << '\n';
                utopia::out() << "shear_modulus:\t" << shear_modulus.get() << '\n';
                utopia::out() << "poisson_ratio:\t" << poisson_ratio.get() << '\n';
                utopia::out() << "young_modulus:\t" << young_modulus.get() << '\n';
            }
        }

        FirstLameParameter<FirstLameParameter_> first_lame_parameter;
        ShearModulus<ShearModulus_> shear_modulus;
        PoissonRatio<PoissonRatio_> poisson_ratio;
        YoungModulus<YoungModulus_> young_modulus;
    };

}  // namespace utopia

#endif  // UTOPIA_STRESS_STRAIN_PARAMETERS_HPP
