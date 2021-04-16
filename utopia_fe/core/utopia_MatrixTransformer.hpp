#ifndef UTOPIA_MATRIX_TRANSFORMER_HPP
#define UTOPIA_MATRIX_TRANSFORMER_HPP

#include "utopia_Input.hpp"
#include "utopia_Traits.hpp"

#include "utopia_StabilizeTransport.hpp"

namespace utopia {
    template <class Matrix>
    class MatrixTransformer : public Configurable {
    public:
        using MatrixTransformerPtr = std::unique_ptr<MatrixTransformer>;
        virtual ~MatrixTransformer() = default;
        virtual void apply(Matrix &mat) const = 0;
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

    template <class Matrix>
    class StabilizeTransport final : public MatrixTransformer<Matrix> {
    public:
        virtual void apply(Matrix &mat) const {
            // using Scalar = typename Traits<Matrix>::Scalar;
            // using SizeType = typename Traits<Matrix>::SizeType;
            // using Vector = typename Traits<Matrix>::Vector;

            // Matrix mat_t = transpose(mat);

            // mat.transform([&](const SizeType i, const SizeType j, const Scalar &value) {
            //     if (i == j) {
            //         return 0.0;
            //     } else {
            //         const Scalar value_t = mat_t.get(i, j);
            //         Scalar max_val = std::max(value, value_t);

            //         if (max_val > 0.0) {
            //             max_val *= -1.0;
            //             return max_val;
            //         } else {
            //             return 0.0;
            //         }
            //     }
            // });

            // Vector diag_elem = -1.0 * sum(mat, 1);
            // mat += diag(diag_elem);

            Matrix stab;
            transport_stabilization(mat, stab);
            mat += stab;
        }
    };

    template <class Matrix>
    void register_transfomers() {
        static bool once = false;

        if (!once) {
            MatrixTransformer<Matrix>::template Register<StabilizeTransport<Matrix>>("StabilizeTransport");
            once = true;
        }
    }

}  // namespace utopia

#endif  // UTOPIA_MATRIX_TRANSFORMER_HPP
