#ifndef UTOPIA_MATRIX_TRANSFORMER_HPP
#define UTOPIA_MATRIX_TRANSFORMER_HPP

#include <memory>
#include "utopia_Input.hpp"
#include "utopia_SimulationTime.hpp"
#include "utopia_Tracer.hpp"
#include "utopia_Traits.hpp"

// #include "utopia_StabilizeTransport.hpp"

namespace utopia {
    template <class Matrix>
    class MatrixTransformer : public Configurable {
    public:
        using Scalar = typename Traits<Matrix>::Scalar;
        using MatrixTransformerPtr = std::unique_ptr<MatrixTransformer>;
        virtual ~MatrixTransformer() = default;
        virtual void apply(Matrix &mat) = 0;
        virtual void apply_to_matrices(Matrix &mass_matrix, Matrix &op) = 0;
        virtual void set_time(const std::shared_ptr<SimulationTime<Scalar>> &time) = 0;
        void read(Input &) override {}

        class Registry {
        public:
            static Registry &instance() {
                static Registry instance_;
                return instance_;
            }

            template <class Type>
            void register_transformer(const std::string &name) {
                transformers[name] = []() -> MatrixTransformerPtr { return utopia::make_unique<Type>(); };
            }

            MatrixTransformerPtr find_transformer(const std::string &name) const {
                MatrixTransformerPtr ret;
                auto it = transformers.find(name);

                if (it != transformers.end()) {
                    ret = it->second();
                }

                return ret;
            }

        private:
            std::map<std::string, std::function<MatrixTransformerPtr()>> transformers;
        };

        static std::unique_ptr<MatrixTransformer> New(const std::string &type) {
            return Registry::instance().find_transformer(type);
        }

        template <class SubType>
        static void Register(const std::string &type) {
            return Registry::instance().template register_transformer<SubType>(type);
        }
    };

    template <class Vector>
    class PostProcessor {
    public:
        virtual ~PostProcessor() = default;
        virtual void post_process(const Vector &rhs, Vector &x) const = 0;
    };

    // template <class Matrix>
    // void register_transfomers() {
    //     static bool once = false;

    //     if (!once) {
    //         MatrixTransformer<Matrix>::template Register<StabilizeTransport<Matrix>>("StabilizeTransport");
    //         once = true;
    //     }
    // }

}  // namespace utopia

#endif  // UTOPIA_MATRIX_TRANSFORMER_HPP
