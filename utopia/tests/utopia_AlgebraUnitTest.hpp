#ifndef UTOPIA_PARALLEL_TEST_RUNNER_HPP
#define UTOPIA_PARALLEL_TEST_RUNNER_HPP

#include "utopia_Input.hpp"
#include "utopia_Instance.hpp"
#include "utopia_MPI.hpp"
#include "utopia_RunParallelTest.hpp"
#include "utopia_Traits.hpp"

namespace utopia {

    template <class Tensor>
    class AlgebraUnitTest : public UnitTest<typename Traits<Tensor>::Communicator> {
    public:
        void print_backend_info() const override {
            if (Utopia::instance().verbose() && mpi_world_rank() == 0) {
                utopia::out() << "\nBackend: " << Traits<Tensor>::backend_info().get_name() << std::endl;
            }
        }

        ~AlgebraUnitTest() override = default;
    };

}  // namespace utopia

#define UTOPIA_TEST_CLASS(macro_ClassName, macro_Matrix_, macro_Vector_) \
    template <class macro_Matrix_, class macro_Vector_>                  \
    class macro_ClassName final : public UnitTest<typename Traits<macro_Vector_>::Communicator>

#define UTOPIA_EXPOSE_TYPES(macro_Tensor_)        \
    using Traits = utopia::Traits<macro_Tensor_>; \
    using Scalar = typename Traits::Scalar;       \
    using SizeType = typename Traits::SizeType;   \
    using IndexSet = typename Traits::IndexSet;   \
    using Comm = typename Traits::Communicator;   \
    using Layout = typename Traits::Layout;       \
    using MatrixLayout = typename Traits::MatrixLayout;

#endif  // UTOPIA_PARALLEL_TEST_RUNNER_HPP
