#include "utopia_mars_Material.hpp"

#include "utopia_kokkos_UniformFE.hpp"
#include "utopia_mars_Discretization.hpp"
#include "utopia_mars_FunctionSpace.hpp"

namespace utopia {

    void Material<utopia::mars::FunctionSpace, mars_FE_t>::initialize(
        const std::shared_ptr<utopia::mars::FunctionSpace> &space) {
        using Scalar_t = utopia::Traits<utopia::mars::FunctionSpace>::Scalar;
        using FE_t = utopia::kokkos::UniformFE<Scalar_t>;

        using Assembler_t = utopia::kokkos::FEAssembler<utopia::mars::FunctionSpace, FE_t>;
        using Discretization_t = utopia::Discretization<utopia::mars::FunctionSpace, FE_t>;

        std::vector<std::shared_ptr<FE_t>> fes;
        space->handler()->create(fes, order());
        auto discretization = std::make_shared<Discretization_t>(space, fes[0]);
        assert(fes[0]->n_cells() > 0);
        assert(fes[0]->n_shape_functions() > 0);

        auto assembler = std::make_shared<Assembler_t>(fes[0]);
        assembler->set_discretization(discretization);
        this->set_assembler(assembler);
    }

}  // namespace utopia