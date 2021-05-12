#ifndef UTOPIA_INTREPID2_L2_NORM_HPP
#define UTOPIA_INTREPID2_L2_NORM_HPP

#include "utopia_Input.hpp"
#include "utopia_Traits.hpp"

#include "utopia_intrepid2_L2ScalarProduct.hpp"

#include <cmath>
#include <vector>

namespace utopia {

    namespace intrepid2 {
        template <class FunctionSpace>
        class L2Norm : public Configurable {
        public:
            using Scalar_t = typename Traits<FunctionSpace>::Scalar;
            using L2ScalarProduct_t = utopia::L2ScalarProduct<Scalar_t>;
            using L2ScalarProductAssembler_t = intrepid2::Assemble<L2ScalarProduct_t>;
            using Intrepid2Field_t = intrepid2::Field<Scalar_t>;

            void read(Input &) override {}

            void compute(const utopia::Field<FunctionSpace> &field, std::vector<Scalar_t> &results) {
                if (!fe_) {
                    fe_ = std::make_shared<intrepid2::FE<Scalar_t>>();
                    create_fe(*field.space(), *fe_, 2);
                }

                int n_components = field.tensor_size();

                assert(n_components >= 1);

                L2ScalarProduct_t op;
                op.n_components = n_components;
                l2_scalar_product_ = std::make_shared<L2ScalarProductAssembler_t>(fe_, op);

                if (!interpid2_field_) {
                    interpid2_field_ = std::make_shared<Intrepid2Field_t>(fe_);
                }

                convert_field(field, *interpid2_field_);

                results.resize(n_components);

                for (int c = 0; c < n_components; ++c) {
                    results[c] = l2_scalar_product_->reduce(interpid2_field_->data(), interpid2_field_->data(), c);
                }

                field.data().comm().sum(n_components, &results[0]);
                for (int c = 0; c < n_components; ++c) {
                    results[c] = std::sqrt(results[c]);
                }
            }

        private:
            std::shared_ptr<intrepid2::FE<Scalar_t>> fe_;
            std::shared_ptr<L2ScalarProductAssembler_t> l2_scalar_product_;
            std::shared_ptr<Intrepid2Field_t> interpid2_field_;
        };

        template <class FunctionSpace>
        auto l2_norm(const utopia::Field<FunctionSpace> &field) -> typename Traits<FunctionSpace>::Scalar {
            using Scalar = typename Traits<FunctionSpace>::Scalar;
            std::vector<Scalar> results;
            L2Norm<FunctionSpace>().compute(field, results);

            assert(results.size() == 1);
            return results[0];
        }

        template <class FunctionSpace>
        void l2_norm(const utopia::Field<FunctionSpace> &field,
                     std::vector<typename Traits<FunctionSpace>::Scalar> &results) {
            L2Norm<FunctionSpace>().compute(field, results);
        }

    }  // namespace intrepid2
}  // namespace utopia

#endif  // UTOPIA_INTREPID2_L2_NORM_HPP
