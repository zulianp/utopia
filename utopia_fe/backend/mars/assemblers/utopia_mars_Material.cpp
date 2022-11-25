#include "utopia_mars_Material.hpp"

#include "utopia_kokkos_UniformFE.hpp"
#include "utopia_mars_Discretization.hpp"
#include "utopia_mars_FunctionSpace.hpp"

namespace utopia {

    void Material<utopia::mars::FunctionSpace, mars_FE_t>::initialize(
        const std::shared_ptr<utopia::mars::FunctionSpace> &space) {
        // using Scalar_t = utopia::Traits<utopia::mars::FunctionSpace>::Scalar;
        // using Assembler_t =
        //     utopia::kokkos::FEAssembler<utopia::mars::FunctionSpace, utopia::kokkos::UniformFE<Scalar_t>>;
        // using Discretization_t =
        //     utopia::Discretization<utopia::mars::FunctionSpace, utopia::kokkos::UniformFE<Scalar_t>>;

        // auto fe_ptr = std::make_shared<utopia::kokkos::UniformFE<Scalar_t>>();
        // create_fe(*space, *fe_ptr, order());
        // auto discretization = std::make_shared<Discretization_t>(space, fe_ptr);
        // auto assembler = std::make_shared<Assembler_t>(fe_ptr);
        // assembler->set_discretization(discretization);
        // this->set_assembler(assembler);
        assert(false);
        Utopia::Abort();
    }

}  // namespace utopia