#ifndef UTOPIA_BELOS_SOLVERS_HPP
#define UTOPIA_BELOS_SOLVERS_HPP

#include "utopia_Config.hpp"
#include "utopia_MatrixFreeLinearSolver.hpp"
#include "utopia_PreconditionedSolver.hpp"

#ifdef UTOPIA_ENABLE_TRILINOS_BELOS

namespace utopia {
    /**@ingroup     Linear
     * @brief       Class provides interface to Trilinos Belos solvers \n
     *              For setting up basic parameters, one can use classic Belos
     *              runtime options
     */
    template <typename Matrix, typename Vector, int Backend = Traits<Matrix>::Backend>
    class BelosSolver {};

    template <typename Matrix, typename Vector>
    class BelosSolver<Matrix, Vector, TRILINOS> : public OperatorBasedLinearSolver<Matrix, Vector> {
    public:
        using Scalar = typename utopia::Traits<Vector>::Scalar;
        using SizeType = typename utopia::Traits<Vector>::SizeType;

        using Preconditioner = utopia::Preconditioner<Vector>;
        using Super = utopia::OperatorBasedLinearSolver<Matrix, Vector>;

        struct Param {
            enum Key {
                BLOCK_SIZE = 0,
                MAX_RESTARTS,
                NUM_BLOCKS,
                OUTPUT_FREQ,
                ORTHOGONALIZATION,
                USE_SINGLE_RED,
                NUM_PARAMS
            };

            const std::string name;
            const std::string belos_name;
            const std::string type;
            const std::string description;
            const std::string default_value;
        };

        using ParamKey = typename Param::Key;

        BelosSolver() = delete;
        BelosSolver(const BelosSolver &other);
        BelosSolver(const std::string &solver_type, std::initializer_list<ParamKey> solver_params = {});
        ~BelosSolver() override;

        void update(const std::shared_ptr<const Matrix> &op, const std::shared_ptr<const Matrix> &prec) override;
        void update(const std::shared_ptr<const Matrix> &op) override;
        bool apply(const Vector &rhs, Vector &lhs) override;

        bool solve(const Operator<Vector> &A, const Vector &rhs, Vector &sol) override;
        void update(const Operator<Vector> &A) override;

        void print_usage(std::ostream &os = std::cout) const override;
        void read(Input &in) override;

        void set_preconditioner(const std::shared_ptr<Preconditioner> &precond) override;
        void set_preconditioner(const std::string &direct_pc_type, const std::string &pc_side = "right");

    private:
        const Param params_[ParamKey::NUM_PARAMS] = {
            {"block_size", "Block Size", "int", "Block size used by solver", std::to_string(deafult_block_size_)},
            {"max_restarts",
             "Maximum Restarts",
             "int",
             "Maximum number of restarts the solver is allowed to perform",
             std::to_string(default_max_restarts_)},
            {"num_blocks",
             "Num Blocks",
             "int",
             "Number of blocks allocated for the Krylov basis",
             std::to_string(deafult_num_blocks_)},
            {"output_freq",
             "Output Frequency",
             "int",
             "How often convergence information should be outputted.",
             std::to_string(default_output_frequency_)},
            {"orthogonalization",
             "Orthogonalization",
             "string",
             "The desired orthogonalization: DGKS ,ICGS, and IMGS",
             default_orthogonalization_},
            {"use_single_red",
             "Use Single Reduction",
             "bool",
             "Use single-reduction variant of solver iteration",
             default_use_single_red_ ? "true" : "false"}};

        static constexpr int deafult_block_size_{1};
        int block_size_{deafult_block_size_};

        static constexpr int default_max_restarts_{20};
        int max_restarts_{default_max_restarts_};

        static constexpr int deafult_num_blocks_{300};
        int num_blocks_{deafult_num_blocks_};

        static constexpr int default_output_frequency_{10};
        int output_frequency_{default_output_frequency_};

        static constexpr char default_orthogonalization_[] = "ICGS";
        std::string orthogonalization_{default_orthogonalization_};

        static constexpr bool default_use_single_red_{false};
        bool use_single_red_{default_use_single_red_};

        const std::string solver_type_;
        const std::set<ParamKey> solver_params_;

        void set_preconditioner(const std::shared_ptr<const Matrix> &op);
        void set_programmatic_config();
        void set_problem();

        bool enable_pc_{false};
        bool enable_direct_pc_{false};
        std::string direct_pc_type_;
        std::string pc_side_;

        class Impl;
        std::unique_ptr<Impl> impl_;
    };

}  // namespace utopia

#else  // UTOPIA_ENABLE_TRILINOS_BELOS
#error "Trilinos was not configured with Belos, hence you cannot use the utopia::BelosSolver."
#endif  // UTOPIA_ENABLE_TRILINOS_BELOS
#endif  // UTOPIA_BELOS_SOLVERS_HPP
