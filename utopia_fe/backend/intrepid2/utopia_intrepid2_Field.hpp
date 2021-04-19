#ifndef UTOPIA_INTREPID2_FIELD_HPP
#define UTOPIA_INTREPID2_FIELD_HPP

#include "utopia_intrepid2_FE.hpp"

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

            inline void set_n_components(const int n_components) { n_components_ = n_components; }
            inline int n_components() const { return n_components_; }
            virtual void ensure_field() { assert(false && "IMPLEMENT ME in subclass"); }

            inline void set_name(const std::string &name) { name_ = name; }

            void scale(const Scalar &a) {
                auto n_elements = data_.size();

                auto data = data_.data();
                Kokkos::parallel_for(
                    n_elements, UTOPIA_LAMBDA(int i) { data[i] *= a; });
            }

        private:
            std::shared_ptr<FE> fe_;
            DynRankView data_;
            int n_components_{1};
            std::string name_{"UnknownField"};
        };
    }  // namespace intrepid2
}  // namespace utopia

#endif  // UTOPIA_INTREPID2_FIELD_HPP
