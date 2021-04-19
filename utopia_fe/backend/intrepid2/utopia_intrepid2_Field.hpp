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

            Field(const std::shared_ptr<FE> &fe) : fe_(fe) {}
            DynRankView &data() { return data_; }

        private:
            std::shared_ptr<FE> fe_;
            DynRankView data_;
            int n_components_{1};
        };
    }  // namespace intrepid2
}  // namespace utopia

#endif  // UTOPIA_INTREPID2_FIELD_HPP
