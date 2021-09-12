#ifndef UTOPIA_KOKKOS_FIELD_HPP
#define UTOPIA_KOKKOS_FIELD_HPP

#include "utopia_Algorithms.hpp"
#include "utopia_MeshElementType.hpp"

#include "utopia_kokkos_FE.hpp"
#include "utopia_kokkos_SubdomainValue.hpp"

#include <memory>

namespace utopia {
    namespace kokkos {
        template <typename FE_>
        class Field : public Describable {
        public:
            using FE = FE_;
            using Scalar = typename FE::Scalar;
            using Function = typename FE::Function;
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

            virtual void scale(const SubdomainValue<FE> &value) {
                auto tags = fe_->element_tags();
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
                assert(fe_->n_shape_functions() * tensor_size_ == data_.extent(1));

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

            void describe(std::ostream &) const override {
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
                UTOPIA_INLINE_FUNCTION Interpolate(const Function &fun, const DynRankView &coefficients)
                    : fun_(fun),
                      coefficients_(coefficients),
                      n_shape_functions_(fun.extent(0)),
                      tensor_size_(coefficients_.extent(1) / n_shape_functions_) {}

                UTOPIA_INLINE_FUNCTION Scalar operator()(const int cell, const int qp, const int var) const {
                    assert(var < tensor_size_);

#ifndef NDEBUG
                    Scalar PU_test = 0.0;
#endif

                    Scalar ret = 0.0;
                    for (int i = 0; i < n_shape_functions_; ++i) {
                        const Scalar c = coefficients_(cell, i * tensor_size_ + var);
                        // assert(device::approxeq(1.0, c, 1e-8));
                        ret += fun_(i, qp) * c;

#ifndef NDEBUG
                        PU_test += fun_(i, qp);
#endif
                    }

                    assert(device::approxeq(1.0, PU_test, 1e-8));

                    return ret;
                }

                UTOPIA_INLINE_FUNCTION Scalar operator()(const int cell, const int qp) const {
                    assert(1 == tensor_size_);

                    Scalar ret = 0.0;
                    for (int i = 0; i < n_shape_functions_; ++i) {
                        ret += fun_(i, qp) * coefficients_(cell, i);
                    }

                    return ret;
                }

                Function fun_;
                DynRankView coefficients_;
                const int n_shape_functions_;
                const int tensor_size_;
            };

            Interpolate interpolate() const {
                assert(fe_);
                return Interpolate(fe_->fun(), data_);
            }

            virtual bool is_coefficient() const { return true; }

            void set_elem_type(MeshElementType type) { elem_type_ = type; }
            inline MeshElementType elem_type() const { return elem_type_; }

        private:
            std::shared_ptr<FE> fe_;
            DynRankView data_;
            int tensor_size_{1};
            std::string name_{"UnknownField"};
            MeshElementType elem_type_{UNDEFINED_TYPE};
        };

        template <typename FE_>
        class QPField : public Field<FE_> {
        public:
            using FE = FE_;
            using Scalar = typename FE::Scalar;
            using Super = utopia::kokkos::Field<FE>;
            using DynRankView = typename FE::DynRankView;
            using Super::scale;

            virtual ~QPField() = default;

            QPField(const std::shared_ptr<FE> &fe) : Super(fe) {}

            void scale(const SubdomainValue<FE> &value) override {
                auto tags = this->fe()->element_tags();
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

            void avg(Field<FE> &field) {
                assert(field.is_coefficient());

                const int tensor_size = this->tensor_size();
                field.data() = DynRankView("avg", this->fe()->n_cells(), tensor_size);
                field.set_tensor_size(tensor_size);

                auto data = this->data();
                auto avg_data = field.data();
                auto measure = this->fe()->measure();

                const int n_quad_points = this->fe()->n_quad_points();

                ::Kokkos::parallel_for(
                    this->fe()->cell_range(), UTOPIA_LAMBDA(int cell) {
                        for (int i = 0; i < tensor_size; ++i) {
                            Scalar m = 0.;
                            for (int qp = 0; qp < n_quad_points; ++qp) {
                                m += measure(cell, qp);
                                avg_data(cell, i) += data(cell, qp, i);
                            }

                            avg_data(cell, i) /= m;
                        }
                    });
            }

            void describe(std::ostream &) const override {
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
    }  // namespace kokkos
}  // namespace utopia

#endif  // UTOPIA_KOKKOS_FIELD_HPP
