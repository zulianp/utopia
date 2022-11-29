#ifndef UTOPIA_ASSEMBLY_KOKKOS_GENERATED_IMPL_HPP
#define UTOPIA_ASSEMBLY_KOKKOS_GENERATED_IMPL_HPP

#include "utopia_assembly_kokkos_generated.hpp"

#include "utopia_material_LaplaceOperator.hpp"
#include "utopia_material_LaplaceOperator_AxisAlignedHex8_3_impl.hpp"
#include "utopia_material_LaplaceOperator_AxisAlignedQuad4_2_impl.hpp"
#include "utopia_material_LaplaceOperator_Hex8_3_impl.hpp"
#include "utopia_material_LaplaceOperator_Line2_1_impl.hpp"
#include "utopia_material_LaplaceOperator_Pentatope5_4_impl.hpp"
#include "utopia_material_LaplaceOperator_Quad4_2_impl.hpp"
#include "utopia_material_LaplaceOperator_Tet4_3_impl.hpp"
#include "utopia_material_LaplaceOperator_Tri3_2_impl.hpp"
#include "utopia_material_Mass.hpp"
#include "utopia_material_Mass_AxisAlignedHex8_3_impl.hpp"
#include "utopia_material_Mass_AxisAlignedQuad4_2_impl.hpp"
#include "utopia_material_Mass_Hex8_3_impl.hpp"
#include "utopia_material_Mass_Line2_1_impl.hpp"
#include "utopia_material_Mass_Pentatope5_4_impl.hpp"
#include "utopia_material_Mass_Quad4_2_impl.hpp"
#include "utopia_material_Mass_Tet4_3_impl.hpp"
#include "utopia_material_Mass_Tri3_2_impl.hpp"

namespace utopia {
    namespace kokkos {

        template <class FunctionSpace, class... Args>
        void register_generated_assemblers(AssemblerRegistry<FunctionSpace, kokkos::FE<Args...>> &) {
            m_utopia_warning_once("register_generated_assemblers does not work for this backend!");
        }

        template <class FunctionSpace, class FE>
        void register_generated_assemblers(AssemblerRegistry<FunctionSpace, FE> &registry) {
            registry.template register_assembler_variant<utopia::kokkos::MassTri3<FunctionSpace, FE>>("MassLine2", 1);
            registry.template register_assembler_variant<utopia::kokkos::MassTri3<FunctionSpace, FE>>("MassTri3", 2);
            registry.template register_assembler_variant<utopia::kokkos::MassQuad4<FunctionSpace, FE>>("MassQuad4", 2);
            registry.template register_assembler_variant<utopia::kokkos::MassAxisAlignedQuad4<FunctionSpace, FE>>(
                "MassAxisAlignedQuad4", 2);

            registry.template register_assembler_variant<utopia::kokkos::MassTet4<FunctionSpace, FE>>("MassTet4", 3);
            registry.template register_assembler_variant<utopia::kokkos::MassHex8<FunctionSpace, FE>>("MassHex8", 3);
            registry.template register_assembler_variant<utopia::kokkos::MassAxisAlignedHex8<FunctionSpace, FE>>("MassAxisAlignedHex8",
                                                                                                  3);

            registry.template register_assembler_variant<utopia::kokkos::LaplaceOperatorTri3<FunctionSpace, FE>>(
                "LaplaceOperatorLine2", 1);
            registry.template register_assembler_variant<utopia::kokkos::LaplaceOperatorTri3<FunctionSpace, FE>>("LaplaceOperatorTri3",
                                                                                                  2);
            registry.template register_assembler_variant<utopia::kokkos::LaplaceOperatorQuad4<FunctionSpace, FE>>(
                "LaplaceOperatorQuad4", 2);
            registry.template register_assembler_variant<utopia::kokkos::LaplaceOperatorAxisAlignedQuad4<FunctionSpace, FE>>(
                "LaplaceOperatorAxisAlignedQuad4", 2);

            registry.template register_assembler_variant<utopia::kokkos::LaplaceOperatorTet4<FunctionSpace, FE>>("LaplaceOperatorTet4",
                                                                                                  3);
            registry.template register_assembler_variant<utopia::kokkos::LaplaceOperatorHex8<FunctionSpace, FE>>("LaplaceOperatorHex8",
                                                                                                  3);
            registry.template register_assembler_variant<utopia::kokkos::LaplaceOperatorAxisAlignedHex8<FunctionSpace, FE>>(
                "LaplaceOperatorAxisAlignedHex8", 3);
        }

    }  // namespace kokkos
}  // namespace utopia

#endif  // UTOPIA_ASSEMBLY_KOKKOS_GENERATED_IMPL_HPP
