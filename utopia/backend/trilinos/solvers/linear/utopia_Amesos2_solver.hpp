#ifndef UTOPIA_AMESOS2_SOLVERS_HPP
#define UTOPIA_AMESOS2_SOLVERS_HPP

#include "utopia_Base.hpp"
#include "utopia_DirectSolver.hpp"

#ifdef UTOPIA_WITH_TRILINOS_AMESOS2
#include "Amesos2_config.h"

#include <variant>

namespace utopia {
    /**@ingroup     Linear
     * @brief       Class provides interface to Trilinos Amesos2 solvers \n
     *              For setting up basic parameters, one can use classic Amesos2
     * runtime options
     */
    template <typename Matrix, typename Vector, int Backend = Traits<Matrix>::Backend>
    class Amesos2Solver {};

    template <typename Matrix, typename Vector>
    class Amesos2Solver<Matrix, Vector, TRILINOS> : public DirectSolver<Matrix, Vector> {
    public:
        class Param {
        public:
            const std::string amesos_name;
            const std::string description;
            const std::string type_str;
            const std::string value_str;

            Param(const char *amesos_name_, const char *descr_, bool value_)
                : amesos_name(amesos_name_),
                  description(descr_),
                  type_str("bool"),
                  value_str(value_ ? "true" : "false"),
                  value(value_) {}

            Param(const char *amesos_name_, const char *descr_, int value_)
                : amesos_name(amesos_name_),
                  description(descr_),
                  type_str("int"),
                  value_str(std::to_string(value_)),
                  value(value_) {}

            Param(const char *amesos_name_, const char *descr_, const char *value_)
                : amesos_name(amesos_name_),
                  description(descr_),
                  type_str("string"),
                  value_str(value_),
                  value(value_) {}

            template <typename T>
            T default_value() const {
                return std::get<T>(value);
            }

        private:
            std::variant<bool, int, std::string> value;
        };

        Amesos2Solver(const Amesos2Solver &other);
        Amesos2Solver() = delete;
        Amesos2Solver(const std::string &solver_type,
                      const std::initializer_list<std::pair<const std::string, Param>> params = {});
        ~Amesos2Solver() override;

        void update(const std::shared_ptr<const Matrix> &op) override;
        bool apply(const Vector &rhs, Vector &lhs) override;

        void print_usage(std::ostream &os) const override;
        void read(Input &is) override;

        int get_nnzLU() const;
        int get_num_preorder() const;
        int get_num_sym_fact() const;
        int get_num_numeric_fact() const;
        int get_num_solve() const;
        bool get_preordering_done() const;
        bool get_sym_factorization_done() const;
        bool get_num_factorization_done() const;

    private:
        const std::string solver_type_;
        const std::map<std::string, Param> params_;

        static constexpr bool default_keep_symbolic_factorization{false};
        bool keep_symbolic_factorization{default_keep_symbolic_factorization};

        class Impl;
        std::unique_ptr<Impl> impl_;

        bool preordering();
        bool num_factorization();
        bool sym_factorization();
    };

}  // namespace utopia

#endif  // UTOPIA_WITH_TRILINOS_AMESOS2
#endif  // UTOPIA_AMESOS2_SOLVERS_HPP