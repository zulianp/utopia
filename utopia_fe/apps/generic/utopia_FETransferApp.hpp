#ifndef UTOPIA_FE_TRANSFER_APP_HPP
#define UTOPIA_FE_TRANSFER_APP_HPP

#include "utopia_fe_base.hpp"

#include "utopia_Agglomerate.hpp"
#include "utopia_BlockAgglomerate.hpp"
#include "utopia_ElementWisePseudoInverse.hpp"
#include "utopia_Field.hpp"
#include "utopia_ILU.hpp"
#include "utopia_MPITimeStatistics.hpp"
#include "utopia_petsc_AdditiveCorrectionTransfer.hpp"
#include "utopia_petsc_DILUAlgorithm.hpp"

#include "utopia_ui.hpp"

#include "utopia_fe_Core.hpp"

#include <memory>

namespace utopia {

    template <class FunctionSpace>
    class FETransferApp : public Configurable {
    public:
        // Extract front-end associated objects
        using Vector_t = typename Traits<FunctionSpace>::Vector;
        using Matrix_t = typename Traits<FunctionSpace>::Matrix;
        using Size_t = typename Traits<FunctionSpace>::SizeType;
        using Scalar_t = typename Traits<FunctionSpace>::Scalar;
        using Communicator_t = typename Traits<FunctionSpace>::Communicator;
        using Mesh_t = typename Traits<FunctionSpace>::Mesh;

        // Use specialized compoenents for function space
        using Obstacle_t = utopia::Obstacle<FunctionSpace>;
        using OmniAssembler_t = utopia::OmniAssembler<FunctionSpace>;

        // Use algorithms from utopia algebra
        using QPSolver_t = utopia::QPSolver<Matrix_t, Vector_t>;
        using SemismoothNewton_t = utopia::SemismoothNewton<Matrix_t, Vector_t>;
        using Factorization_t = utopia::Factorization<Matrix_t, Vector_t>;
        using LinearSolver_t = utopia::LinearSolver<Matrix_t, Vector_t>;
        using IterativeSolver_t = utopia::IterativeSolver<Matrix_t, Vector_t>;
        using AlgebraicMultigrid_t = utopia::AlgebraicMultigrid<Matrix_t, Vector_t>;
        using ProjectedGaussSeidel_t = utopia::ProjectedGaussSeidel<Matrix_t, Vector_t>;
        using KSPSolver_t = utopia::KSPSolver<Matrix_t, Vector_t>;

        void read(Input &in) override {
            in.get("from", [this](Input &in) {
                bool read_state = false;
                in.get("read_state", read_state);

                if (read_state) {
                    from_space.read_with_state(in, field);

                    const Scalar_t norm_field = norm2(field.data());
                    std::cout << "norm_field: " << norm_field << std::endl;
                } else {
                    from_space.read(in);
                }
            });

            if (from_space.empty()) {
                return;
            }

            in.get("to", to_space);

            if (to_space.empty()) {
                return;
            }

            if (field.empty()) {
                from_space.create_field(field);
                field.data().set(1.0);
            }

            in.get("transfer", transfer);
            in.get("output_path", output_path);
        }

        void run() {
            if (!transfer.init(make_ref(from_space), make_ref(to_space))) {
                return;
            }

            Vector_t to;

            transfer.apply(field.data(), to);
            to_space.write(output_path, to);
        }

        FETransferApp() {}

        bool is_valid() const { return !from_space.empty() && !to_space.empty() && !field.empty(); }

    private:
        FunctionSpace from_space;
        FunctionSpace to_space;
        Field<FunctionSpace> field;
        FETransfer<FunctionSpace> transfer;
        Path output_path{"./out.e"};
    };

}  // namespace utopia

#endif  // UTOPIA_FE_TRANSFER_APP_HPP