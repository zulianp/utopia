#ifndef UTOPIA_SOLVER_PRINT_Info_HPP
#define UTOPIA_SOLVER_PRINT_Info_HPP

#include <chrono>
#include <iomanip>
#include <limits>
#include "utopia_Traits.hpp"
#include "utopia_Utils.hpp"

namespace utopia {
    /**
     * @brief      The class helping to print-out information about solver:
     * initialization messages, prinitning iteration status, time-stats and
     * exit/convergance messages. It also helps pass solution and informations about
     * solve back into FEM packages.
     */
    class PrintInfo {
        using Scalar = double;

    public:
        static void print_init(const std::string& method, const std::vector<std::string> status_variables) {
            if (mpi_world_rank() == 0) {
                std::cout << "  \n";
                std::cout << std::setw(10) << std::right << std::string(80, '_') << std::setw(15) << std::endl;
                std::cout << std::setw(45) << std::right << " utopia:: " << method << std::setw(15) << std::endl;
                std::cout << std::setw(10) << std::right << std::string(80, '_') << std::setw(15) << std::endl;
                std::cout << std::endl;

                for (const auto& status_variable : status_variables)
                    std::cout << std::setw(27) << std::right << status_variable;
                std::cout << std::endl;

                auto n = status_variables.size();
                for (size_t i = 0; i < n; i++) std::cout << std::setw(27) << std::right << std::string(10, '-');
                std::cout << std::endl;
            }
        }

        static void print_iter_status(const std::vector<Scalar> status_variables) {
            if (mpi_world_rank() == 0) {
                for (auto status_variable : status_variables)
                    std::cout << std::setw(27) << std::right << std::scientific << status_variable;
                std::cout << std::endl;
            }
        }

        template <class T>
        static void print_iter_status(const T& it, const std::vector<Scalar> scalar_vars) {
            if (mpi_world_rank() == 0) {
                std::cout.precision(15);
                std::cout << std::setw(27) << std::right << it;
                for (auto scalar_var : scalar_vars)
                    std::cout << std::setw(27) << std::right << std::scientific << scalar_var;
                std::cout << std::endl;
            }
        }

    private:
        PrintInfo() = default;
    };

}  // namespace utopia
#endif  // UTOPIA_SOLVER_PRINT_Info_HPP
