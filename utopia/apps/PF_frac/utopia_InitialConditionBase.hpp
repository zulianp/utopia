#ifndef UTOPIA_INITIAL_CONDITION_BASE_HPP
#define UTOPIA_INITIAL_CONDITION_BASE_HPP

#include "utopia_Base.hpp"
// FIXME: this file causes nvcc to fail
#ifndef KOKKOS_ENABLE_CUDA
#include "utopia_Core.hpp"

namespace utopia {

    template <class FunctionSpace>
    class InitialCondition : public Configurable {
    public:
        InitialCondition(FunctionSpace &space) : space_(space) {}

        ~InitialCondition() override = default;

        void read(Input &in) override {}

        virtual void init(PetscVector &x) = 0;
        virtual void init(PetscVector &sol_vec, PetscVector & /*press_vec*/) { this->init(sol_vec); }

    protected:
        FunctionSpace &space_;
    };
}  // namespace utopia

#endif

#endif
