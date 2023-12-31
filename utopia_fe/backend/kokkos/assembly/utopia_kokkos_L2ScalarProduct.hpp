#ifndef UTOPIA_KOKKOS_L2_SCALAR_PRODUCT_HPP
#define UTOPIA_KOKKOS_L2_SCALAR_PRODUCT_HPP

#include "utopia_kokkos_LaplaceOperator.hpp"
#include "utopia_kokkos_SubdomainValue.hpp"

#include "utopia_Views.hpp"

namespace utopia {

    namespace kokkos {

        template <class FE_, typename Fun = typename FE_::Scalar>
        class L2ScalarProduct {
        public:
            using FE = FE_;
            using Scalar = typename Traits<Fun>::Scalar;
            using SizeType = typename FE::SizeType;
            using DynRankView = typename FE::DynRankView;
            using ExecutionSpace = typename FE::ExecutionSpace;
            using Interpolate = typename Field<FE>::Interpolate;

            class Params : public Configurable {
            public:
                using Scalar = typename Traits<Fun>::Scalar;
                // FIXME this type should be generalized to any backend
                using DensityFunction = utopia::kokkos::SubdomainValue<FE>;

                void read(Input &in) override {
                    in.get("density", density);
                    in.get("n_components", n_components);
                    in.get("verbose", verbose);

                    in.get("density_function", [this](Input &in) {
                        density_function = std::make_shared<DensityFunction>(1.0);
                        density_function->read(in);
                    });

                    if (verbose) {
                        utopia::out() << "-----------------------------\n";
                        utopia::out() << "L2ScalarProduct\n";
                        utopia::out() << "density:\t" << density << '\n';
                        utopia::out() << "n_components:\t" << n_components << '\n';
                        utopia::out() << "-----------------------------\n";
                    }
                }

                Params(const Fun &density = Fun(1.0)) : density(density) {}
                UTOPIA_FUNCTION Params(const Params &) = default;

                Fun density;
                int n_components{1};
                std::shared_ptr<DensityFunction> density_function;

                // Testing an printing
                bool verbose{false};
            };

            L2ScalarProduct(const std::shared_ptr<FE> &fe, Params op = Params()) : fe_(fe), op_(std::move(op)) {}

            inline int n_vars() const { return op_.n_components; }
            inline std::string name() const { return "L2ScalarProduct"; }

            class Op {
            public:
                UTOPIA_INLINE_FUNCTION Op(const Fun &density,
                                          const DynRankView &fun,
                                          const DynRankView &x,
                                          const DynRankView &y,
                                          const DynRankView &measure,
                                          const int n_components)
                    : density(density),
                      x(fun, x),
                      y(fun, y),
                      measure(measure),
                      n_components(n_components),
                      n_qp(measure.extent(1)) {}

                UTOPIA_INLINE_FUNCTION Scalar operator()(const int &cell, const int sub_i) const {
                    assert(sub_i < n_components);

                    Scalar ret = 0.0;
                    // Scalar vol = 0.0;

                    for (int qp = 0; qp < n_qp; ++qp) {
                        auto dX = measure(cell, qp);

                        auto x_qp = x(cell, qp, sub_i);
                        auto y_qp = y(cell, qp, sub_i);

                        // assert(device::approxeq(1., x_qp, 1e-8));
                        // assert(device::approxeq(1., y_qp, 1e-8));

                        ret += x_qp * y_qp * density * dX;
                        // vol += dX;
                    }

                    // assert(device::approxeq(vol, ret, 1e-8));

                    assert(ret == ret);
                    return ret;
                }

                inline int dim() const { return n_components; }

                const Fun density;
                const Interpolate x;
                const Interpolate y;
                const DynRankView measure;
                const int n_components;
                const int n_qp;
            };

            inline Op make_op(const DynRankView &x, const DynRankView &y) const {
                return Op(op_.density, fe_->fun(), x, y, fe_->measure(), op_.n_components);
            }

            inline void ensure_buffer() { products = DynRankView("products", fe_->n_cells(), op_.n_components); }

            void apply(const DynRankView &x, const DynRankView &y) {
                ensure_buffer();
                auto op = make_op(x, y);

                auto products = this->products;

                ::Kokkos::parallel_for(
                    "L2ScalarProduct::apply", fe_->cell_range(), UTOPIA_LAMBDA(const int i) {
                        for (int c = 0; c < op.n_components; ++c) {
                            products(i, c) = op(i, c);
                        }
                    });
            }

            Scalar reduce(const DynRankView &x, const DynRankView &y, const int c) const {
                auto op = make_op(x, y);

                assert(c < op_.n_components);

                Scalar ret = 0.0;
                ::Kokkos::parallel_reduce(
                    "L2ScalarProduct::reduce",
                    fe_->cell_range(),
                    UTOPIA_LAMBDA(const int i, Scalar &val) { val += op(i, c); },
                    ret);

                return ret;
            }

            inline Params &user_op() { return op_; }
            inline const Params &user_op() const { return op_; }

            // NVCC_PRIVATE :
            std::shared_ptr<FE> fe_;
            Params op_;
            DynRankView products;
        };
    }  // namespace kokkos
}  // namespace utopia

#endif  // UTOPIA_KOKKOS_L2_SCALAR_PRODUCT_HPP
