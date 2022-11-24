#ifndef UTOPIA_STK_INTREPID2_MATERIAL_HPP
#define UTOPIA_STK_INTREPID2_MATERIAL_HPP

#include "utopia_Material.hpp"
#include "utopia_kokkos_Material.hpp"

#include "utopia_intrepid2_FE.hpp"
#include "utopia_stk_FunctionSpace.hpp"
#include "utopia_stk_intrepid2_Assembler.hpp"

namespace utopia {
    using stk_FE_t = utopia::intrepid2::FE<Traits<utopia::stk::FunctionSpace>::Scalar>;

    template <>
    class Material<utopia::stk::FunctionSpace, stk_FE_t>
        : public utopia::kokkos::Material<utopia::stk::FunctionSpace,
                                          stk_FE_t,
                                          utopia::kokkos::FEAssembler<utopia::stk::FunctionSpace, stk_FE_t>> {
    public:
        void initialize(const std::shared_ptr<utopia::stk::FunctionSpace> &space);
        // FIXME
        inline int order() { return 2; }
    };

}  // namespace utopia

#endif  // UTOPIA_STK_INTREPID2_MATERIAL_HPP
