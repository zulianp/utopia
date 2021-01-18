#ifndef UTOPIA_NONLINEAR_ML_BASE_HPP
#define UTOPIA_NONLINEAR_ML_BASE_HPP
#include "utopia_Core.hpp"
#include "utopia_Function.hpp"
#include "utopia_IPRTransfer.hpp"
#include "utopia_Level.hpp"
#include "utopia_MultiLevelBase.hpp"
#include "utopia_MultiLevelEvaluations.hpp"
#include "utopia_NonlinearMultiLevelInterface.hpp"
#include "utopia_SolutionStatus.hpp"
#include "utopia_Transfer.hpp"

#include <algorithm>
#include <vector>

namespace utopia {
#define CHECK_NUM_PRECISION_mode

    /**
     * @brief      Base class for all nonlinear multilevel solvers. \n
     *             Takes care of inializing multilevel hierarchy - calls into
     * assembly routines on each level. \n Different levels are created by
     * interpolation and restriction.\n
     *
     * @tparam     Matrix
     * @tparam     Vector
     */
    template <class Matrix, class Vector>
    class NonlinearMultiLevelBase : public NonlinearMultiLevelInterface<Matrix, Vector>,
                                    public NonLinearSolver<Vector> {
    public:
        using Scalar = typename Traits<Vector>::Scalar;
        using SizeType = typename Traits<Vector>::SizeType;
        using Layout = typename Traits<Vector>::Layout;

        using NonlinearMultiLevelInterface<Matrix, Vector>::set_transfer_operators;

        NonlinearMultiLevelBase(const SizeType &n_levels) : NonlinearMultiLevelInterface<Matrix, Vector>(n_levels) {
            this->n_levels(n_levels);
        }

        ~NonlinearMultiLevelBase() override = default;

        void read(Input &in) override {
            NonlinearMultiLevelInterface<Matrix, Vector>::read(in);
            NonLinearSolver<Vector>::read(in);
        }

        void print_usage(std::ostream &os) const override {
            NonlinearMultiLevelInterface<Matrix, Vector>::print_usage(os);
            NonLinearSolver<Vector>::print_usage(os);
        }

        /**
         * @brief      Writes CSV file with iteration info
         *
         * @param[in]  it_global  The iterator global
         */
        void print_statistics(const SizeType &it_global) override {
            std::string path = "log_output_path";
            auto non_data_path = Utopia::instance().get(path);

            if (!non_data_path.empty()) {
                CSVWriter writer{};
                if (mpi_world_rank() == 0) {
                    if (!writer.file_exists(non_data_path)) {
                        writer.open_file(non_data_path);
                        writer.write_table_row<std::string>({"v_cycles", "time"});
                    } else
                        writer.open_file(non_data_path);

                    writer.write_table_row<Scalar>({Scalar(it_global), this->get_time()});
                    writer.close_file();
                }
            }
        }

        /**
         * @brief      The solve function for nonlinear multilevel solvers.
         *
         * @param[in]  rhs   Function to be minimized.
         * @param      x_h   The initial guess.
         *
         */
        virtual bool solve(Vector &x_h) = 0;
    };

}  // namespace utopia

#endif  // UTOPIA_NONLINEAR_ML_BASE_HPP
