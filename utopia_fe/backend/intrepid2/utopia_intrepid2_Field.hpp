#ifndef UTOPIA_INTREPID2_FIELD_HPP
#define UTOPIA_INTREPID2_FIELD_HPP

#include "utopia_intrepid2_FE.hpp"
#include "utopia_intrepid2_SubdomainFunction.hpp"

#include <memory>

namespace utopia {
    namespace intrepid2 {
        template <typename Scalar_>
        class Field {
        public:
            using Scalar = Scalar_;
            using FE = utopia::intrepid2::FE<Scalar>;
            using DynRankView = typename FE::DynRankView;

            virtual ~Field() = default;

            Field(const std::shared_ptr<FE> &fe) : fe_(fe) {}
            inline DynRankView &data() { return data_; }
            inline const DynRankView &data() const { return data_; }
            inline std::shared_ptr<FE> fe() { return fe_; }

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
                const int tensor_size = tensor_size_;

                auto data = data_;
                ::Kokkos::parallel_for(
                    fe_->cell_range(), UTOPIA_LAMBDA(int cell) {
                        const auto v = value.value(tags(cell));

                        for (int i = 0; i < tensor_size; ++i) {
                            data(cell, i) *= v;
                        }
                    });
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

            bool is_coefficient() const override { return false; }
        };
    }  // namespace intrepid2
}  // namespace utopia

#endif  // UTOPIA_INTREPID2_FIELD_HPP
