#ifndef UTOPIA_KOKKOS_L2_NORM_HPP
#define UTOPIA_KOKKOS_L2_NORM_HPP

#include "utopia_Input.hpp"
#include "utopia_Traits.hpp"

#include "utopia_kokkos_L2ScalarProduct.hpp"

#include <cmath>
#include <vector>

namespace utopia {

    namespace kokkos {
        template <class FunctionSpace, class FE>
        class L2Norm : public Configurable {
        public:
            using Scalar_t = typename Traits<FunctionSpace>::Scalar;
            using L2ScalarProduct_t = utopia::kokkos::L2ScalarProduct<FE>;
            using Intrepid2Field_t = kokkos::Field<FE>;

            void read(Input &) override {}

            void compute(const utopia::Field<FunctionSpace> &field, std::vector<Scalar_t> &results) {
                if (!fe_) {
                    fe_ = std::make_shared<FE>();
                    create_fe(*field.space(), *fe_, 2);
                }

                int n_components = field.tensor_size();

                assert(n_components >= 1);

                typename L2ScalarProduct_t::Params op;
                op.n_components = n_components;
                l2_scalar_product_ = std::make_shared<L2ScalarProduct_t>(fe_, op);

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
            std::shared_ptr<FE> fe_;
            std::shared_ptr<L2ScalarProduct_t> l2_scalar_product_;
            std::shared_ptr<Intrepid2Field_t> interpid2_field_;
        };

    }  // namespace kokkos
}  // namespace utopia

#endif  // UTOPIA_KOKKOS_L2_NORM_HPP
