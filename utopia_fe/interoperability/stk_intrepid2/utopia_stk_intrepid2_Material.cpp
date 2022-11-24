#include "utopia_stk_intrepid2_Material.hpp"

#include "utopia_intrepid2_FE.hpp"
#include "utopia_stk_FunctionSpace.hpp"
#include "utopia_stk_intrepid2.hpp"
#include "utopia_stk_intrepid2_Discretization.hpp"

namespace utopia {

    void Material<utopia::stk::FunctionSpace, stk_FE_t>::initialize(
        const std::shared_ptr<utopia::stk::FunctionSpace> &space) {
        using Scalar_t = utopia::Traits<utopia::stk::FunctionSpace>::Scalar;
        using Assembler_t = utopia::kokkos::FEAssembler<utopia::stk::FunctionSpace, utopia::intrepid2::FE<double>>;
        using Discretization_t = utopia::Discretization<utopia::stk::FunctionSpace, utopia::intrepid2::FE<double>>;

        auto fe_ptr = std::make_shared<utopia::intrepid2::FE<Scalar_t>>();
        create_fe(*space, *fe_ptr, order());
        auto discretization = std::make_shared<Discretization_t>(space, fe_ptr);
        auto assembler = std::make_shared<Assembler_t>(fe_ptr);
        assembler->set_discretization(discretization);
        this->set_assembler(assembler);
    }

}  // namespace utopia