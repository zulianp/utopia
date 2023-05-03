#ifndef UTOPIA_KOKKOS_MATERIALFACTORY_IMPL_HPP
#define UTOPIA_KOKKOS_MATERIALFACTORY_IMPL_HPP

#include "utopia_kokkos_MaterialFactory.hpp"

#include "utopia_kokkos_ForcingFunction_new.hpp"
#include "utopia_kokkos_LaplaceOperator_new.hpp"
#include "utopia_kokkos_Mass_new.hpp"

#include "utopia_hyperelasticity_NeoHookeanOgden_2.hpp"
#include "utopia_hyperelasticity_NeoHookeanOgden_3.hpp"

#include "utopia_hyperelasticity_SaintVenantKirchoff_2.hpp"
#include "utopia_hyperelasticity_SaintVenantKirchoff_3.hpp"

#include "utopia_hyperelasticity_Yeoh_2.hpp"
#include "utopia_hyperelasticity_Yeoh_3.hpp"

#include "utopia_kokkos_AutoHyperElasticityNew.hpp"

#include "utopia_kokkos_TransportNew_impl.hpp"

#include <memory>
#include <string>

namespace utopia {
    namespace kokkos {

        template <class FunctionSpace, class FE>
        class MaterialFactory<FunctionSpace, FE>::Impl {
        public:
            using Scalar_t = typename Traits<FunctionSpace>::Scalar;

            using Material_t = utopia::Material<FunctionSpace, FE>;
            // using MaterialPtr_t = std::unique_ptr<Material_t>;
            using MaterialPtr_t = std::unique_ptr<utopia::AbstractMaterial<FunctionSpace>>;

            MaterialPtr_t make(const int ndims, const std::string name) {
                std::string type = name;

                if (has_variant.find(type) != has_variant.end()) {
                    type = create_variant_name(ndims, type);
                }

                auto it = materials.find(type);

                if (it == materials.end()) {
                    assert(false);
                    return nullptr;
                }

                auto mat = it->second();
                return mat;
            }

            void register_rhs_materials() {
                register_material<ForcingFunctionNew<FunctionSpace, FE>>("ForcingFunction");
            }

            void register_materials() {
                register_material<LaplaceOperatorNew<FunctionSpace, FE>>("LaplaceOperator");
                register_material<MassNew<FunctionSpace, FE>>("Mass");
                register_porous_media_materials();
                register_hyperelastic_materials();
            }

            void register_porous_media_materials() {
                //
                register_material<TransportNew<FunctionSpace, FE>>("Transport");
            }

            void register_hyperelastic_materials() {
                register_hyperelastic_material<utopia::kernels::NeoHookeanOgden<Scalar_t, 3>>(3, "NeoHookeanOgden");
                register_hyperelastic_material<utopia::kernels::NeoHookeanOgden<Scalar_t, 2>>(2, "NeoHookeanOgden");

                register_hyperelastic_material<utopia::kernels::SaintVenantKirchoff<Scalar_t, 3>>(
                    3, "SaintVenantKirchoff");
                register_hyperelastic_material<utopia::kernels::SaintVenantKirchoff<Scalar_t, 2>>(
                    2, "SaintVenantKirchoff");

                register_hyperelastic_material<utopia::kernels::Yeoh<Scalar_t, 3>>(3, "Yeoh");
                register_hyperelastic_material<utopia::kernels::Yeoh<Scalar_t, 2>>(2, "Yeoh");
            }

            template <class Material>
            void register_material(std::string name) {
                materials[name] = []() -> MaterialPtr_t { return utopia::make_unique<Material>(); };
            }

            std::string create_variant_name(const int ndims, const std::string &name) {
                return name + std::to_string(ndims);
            }

            template <class Material>
            void register_material_variant(const int ndims, std::string name) {
                register_material<Material>(create_variant_name(ndims, name));
                has_variant.insert(name);
            }

            template <class Kernel>
            void register_hyperelastic_material(const int ndims, std::string name) {
                register_material_variant<utopia::kokkos::AutoHyperElasticityNew<FunctionSpace, FE, Kernel>>(ndims,
                                                                                                             name);
            }

            Impl() { register_materials(); }

            std::map<std::string, std::function<MaterialPtr_t()>> materials;
            std::set<std::string> has_variant;
        };

        template <class FunctionSpace, class FE>
        std::unique_ptr<utopia::AbstractMaterial<FunctionSpace>> MaterialFactory<FunctionSpace, FE>::make(
            const int ndims,
            const std::string &name) {
            static Impl impl;
            return impl.make(ndims, name);
        }

    }  // namespace kokkos
}  // namespace utopia

#endif  // UTOPIA_KOKKOS_MATERIALFACTORY_IMPL_HPP