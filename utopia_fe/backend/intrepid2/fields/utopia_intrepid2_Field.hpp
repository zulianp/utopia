#ifndef UTOPIA_INTREPID2_FIELD_HPP
#define UTOPIA_INTREPID2_FIELD_HPP

#include "utopia_intrepid2_FE.hpp"
#include "utopia_intrepid2_SubdomainFunction.hpp"

#include <memory>

namespace utopia {
    namespace intrepid2 {
        template <typename Scalar_>
        class Field : public Describable {
        public:
            using Scalar = Scalar_;
            using FE = utopia::intrepid2::FE<Scalar>;
            using DynRankView = typename FE::DynRankView;

            virtual ~Field() = default;

            Field(const std::shared_ptr<FE> &fe) : fe_(fe) {}
            inline DynRankView &data() { return data_; }
            inline const DynRankView &data() const { return data_; }
            inline std::shared_ptr<FE> fe() { return fe_; }
            inline std::shared_ptr<FE> fe() const { return fe_; }

            virtual std::string name() const { return name_; }

            inline void set_tensor_size(const int tensor_size) { tensor_size_ = tensor_size; }
            inline int tensor_size() const { return tensor_size_; }
            virtual void ensure_field() { assert(false && "IMPLEMENT ME in subclass"); }

            inline void set_name(const std::string &name) { name_ = name; }

            virtual void scale(const Scalar &a) {
                auto data = data_.data();
                auto n_elements = data_.size();

                Kokkos::parallel_for(
                    fe_->range(0, n_elements), UTOPIA_LAMBDA(int i) { data[i] *= a; });
            }

            virtual void scale(const SubdomainValue<Scalar> &value) {
                auto tags = fe_->element_tags;
                const int n_values = data_.extent(1);

                auto data = data_;
                ::Kokkos::parallel_for(
                    fe_->cell_range(), UTOPIA_LAMBDA(int cell) {
                        const auto v = value.value(tags(cell));

                        for (int i = 0; i < n_values; ++i) {
                            data(cell, i) *= v;
                        }
                    });
            }

            Scalar dot(const DynRankView &coeff) const {
                assert(data_.extent(0) == coeff.extent(0));
                assert(data_.extent(1) == coeff.extent(1));
                assert(fe_->num_fields() * tensor_size_ == data_.extent(1));

                const int n_values = data_.extent(1);
                auto data = data_;

                Scalar ret = 0.0;
                ::Kokkos::parallel_reduce(
                    fe_->cell_range(),
                    UTOPIA_LAMBDA(int cell, Scalar &val) {
                        Scalar temp = 0.0;
                        for (int d = 0; d < n_values; ++d) {
                            temp = data(cell, d) * coeff(cell, d);
                        }

                        val += temp;
                    },
                    ret);

                return ret;
            }

            void describe(std::ostream &os) const override {
                const int tensor_size = tensor_size_;

                auto data = data_;
                ::Kokkos::parallel_for(
                    fe_->cell_range(), UTOPIA_LAMBDA(int cell) {
                        for (int i = 0; i < tensor_size; ++i) {
                            printf("%g ", data(cell, i));
                        }

                        printf("\n");
                    });
            }

            class Interpolate {
            public:
                UTOPIA_INLINE_FUNCTION Interpolate(const DynRankView &fun, const DynRankView &coefficients)
                    : fun_(fun),
                      coefficients_(coefficients),
                      num_fields_(fun.extent(1)),
                      tensor_size_(coefficients_.extent(1)) {}

                UTOPIA_INLINE_FUNCTION Scalar operator()(const int cell, const int qp, const int var) const {
                    assert(var < tensor_size_);

                    Scalar ret = 0.0;
                    for (int i = 0; i < num_fields_; ++i) {
                        ret += fun_(cell, i, qp) * coefficients_(cell, i * tensor_size_ + var);
                    }

                    return ret;
                }

                UTOPIA_INLINE_FUNCTION Scalar operator()(const int cell, const int qp) const {
                    assert(1 == tensor_size_);

                    Scalar ret = 0.0;
                    for (int i = 0; i < num_fields_; ++i) {
                        ret += fun_(cell, i, qp) * coefficients_(cell, i);
                    }

                    return ret;
                }

                DynRankView fun_;
                DynRankView coefficients_;
                const int tensor_size_;
                const int num_fields_;
            };

            Interpolate interpolate() const {
                assert(fe_);
                return Interpolate(fe_->fun, data_);
            }

            virtual bool is_coefficient() const { return true; }

        private:
            std::shared_ptr<FE> fe_;
            DynRankView data_;
            int tensor_size_{1};
            std::string name_{"UnknownField"};
        };

        template <typename Scalar_>
        class QPField : public Field<Scalar_> {
        public:
            using Super = utopia::intrepid2::Field<Scalar_>;
            using Scalar = Scalar_;
            using FE = utopia::intrepid2::FE<Scalar>;
            using DynRankView = typename FE::DynRankView;
            using Super::scale;

            virtual ~QPField() = default;

            QPField(const std::shared_ptr<FE> &fe) : Super(fe) {}

            void scale(const SubdomainValue<Scalar> &value) override {
                auto tags = this->fe()->element_tags;
                const int tensor_size = this->tensor_size();

                auto data = this->data();
                ::Kokkos::parallel_for(
                    this->fe()->cell_qp_range(), UTOPIA_LAMBDA(int cell, int qp) {
                        const auto v = value.value(tags(cell));

                        for (int i = 0; i < tensor_size; ++i) {
                            data(cell, qp, i) *= v;
                        }
                    });
            }

            void describe(std::ostream &os) const override {
                const int tensor_size = this->tensor_size();

                auto data = this->data();
                ::Kokkos::parallel_for(
                    this->fe()->cell_qp_range(), UTOPIA_LAMBDA(int cell, int qp) {
                        if (qp == 0) {
                            printf("\n");
                        }

                        for (int i = 0; i < tensor_size; ++i) {
                            printf("%g ", data(cell, qp, i));
                        }

                        printf("\n");
                    });
            }

            bool is_coefficient() const override { return false; }
        };
    }  // namespace intrepid2
}  // namespace utopia

#endif  // UTOPIA_INTREPID2_FIELD_HPP
