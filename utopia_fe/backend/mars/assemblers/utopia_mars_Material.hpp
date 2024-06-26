#ifndef UTOPIA_MARS_INTREPID2_MATERIAL_HPP
#define UTOPIA_MARS_INTREPID2_MATERIAL_HPP

#include "utopia_Material.hpp"
#include "utopia_kokkos_Material.hpp"

#include "utopia_kokkos_FEAssembler.hpp"
#include "utopia_kokkos_UniformFE.hpp"
#include "utopia_mars_FunctionSpace.hpp"

namespace utopia {
    using mars_FE_t = utopia::kokkos::UniformFE<Traits<utopia::mars::FunctionSpace>::Scalar>;

    template <>
    class Material<utopia::mars::FunctionSpace, mars_FE_t>
        : public utopia::kokkos::Material<utopia::mars::FunctionSpace,
                                          mars_FE_t,
                                          utopia::kokkos::FEAssembler<utopia::mars::FunctionSpace, mars_FE_t>> {
    public:
        void initialize(const std::shared_ptr<utopia::mars::FunctionSpace> &space);
        // FIXME
        // inline int order() { return 2; }
    };

}  // namespace utopia

#endif  // UTOPIA_MARS_INTREPID2_MATERIAL_HPP
