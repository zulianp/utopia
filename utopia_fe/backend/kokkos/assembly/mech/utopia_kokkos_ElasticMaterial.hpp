#ifndef UTOPIA_KOKKOS_ELASTIC_MATERIAL_HPP
#define UTOPIA_KOKKOS_ELASTIC_MATERIAL_HPP

namespace utopia {
    namespace kokkos {

        template <class FE_>
        class ElasticMaterial : public FEAssembler<FE_, DefaultView<typename FE_::Scalar>> {
        public:
            using Super = utopia::kokkos::FEAssembler<FE_, DefaultView<typename FE_::Scalar>>;
            using FE = FE_;
            using DynRankView = typename FE::DynRankView;
            using VectorView = typename Super::VectorView;
            using Super::Super;

            virtual ~ElasticMaterial() = default;
            virtual bool principal_stresses(const VectorView &x, VectorView &y) = 0;
        };

    }  // namespace kokkos
}  // namespace utopia

#endif  // UTOPIA_KOKKOS_ELASTIC_MATERIAL_HPP
