#ifndef UTOPIA_BACKEND_PRECONDITIONED_SOLVER_HPP
#define UTOPIA_BACKEND_PRECONDITIONED_SOLVER_HPP

#include <array>
#include <string>

namespace utopia {

    class BackendPreconditionedSolver {
    public:
        enum class PreconditionerSide { INVALID = 0, LEFT, RIGHT };
        enum class PreconditionerType { NONE = 0, CHEBYSHEV, ILUT, JACOBI, MULTIGRID, NUM_TYPES };

        BackendPreconditionedSolver() = delete;
        BackendPreconditionedSolver(const BackendPreconditionedSolver &other)
            : backend_pc_types_(other.backend_pc_types_) {}
        BackendPreconditionedSolver(const std::initializer_list<std::pair<PreconditionerType, std::string>> &pc_types) {
            for (const auto &[pc_type, backend_pc_type] : pc_types) {
                std::string backend_pc_type_lc{backend_pc_type};
                std::transform(
                    backend_pc_type.cbegin(), backend_pc_type.cend(), backend_pc_type_lc.begin(), [](char c) {
                        return std::tolower(c);
                    });
                auto pc_type_index = static_cast<size_t>(pc_type);
                backend_pc_types_.at(pc_type_index) = backend_pc_type_lc;
            }
        }

        void set_preconditioner(const std::string &pc_type, const std::string &pc_side) {
            // map backend-independent name string to enum value
            static const char *pc_type_names[] = {"none", "chebyshev", "ilut", "jacobi", "multigrid"};
            static_assert(std::size(pc_type_names) == NUM_PRECONDITIONER_TYPES);

            PreconditionerType pc_type_{PreconditionerType::NONE};
            if (!pc_type.empty()) {
                auto pc_type_lc{pc_type};
                std::transform(
                    pc_type.cbegin(), pc_type.cend(), pc_type_lc.begin(), [](char c) { return std::tolower(c); });
                for (size_t i = 0; i < std::size(pc_type_names) && pc_type_ == PreconditionerType::NONE; i++) {
                    if (pc_type_lc.compare(pc_type_names[i]) == 0) {
                        pc_type_ = static_cast<PreconditionerType>(i);
                    }
                }
                // should not happen: string for preconditioner type not found
                assert(pc_type_ != PreconditionerType::NONE);
            }

            PreconditionerSide pc_side_{PreconditionerSide::INVALID};
            if (pc_type_ != PreconditionerType::NONE) {
                auto pc_side_lc{pc_side};
                std::transform(
                    pc_side.cbegin(), pc_side.cend(), pc_side_lc.begin(), [](char c) { return std::tolower(c); });
                if (pc_side_lc.compare("left") == 0) {
                    pc_side_ = PreconditionerSide::LEFT;
                } else if (pc_side_lc.compare("right") == 0) {
                    pc_side_ = PreconditionerSide::RIGHT;
                } else {
                    // should not happen: string should be either 'left' or 'right'
                    assert(false);
                }
            }

            // call to backend-specialized method
            set_preconditioner(pc_type_, pc_side_);
        }

        std::string get_backend_preconditioner_name(PreconditionerType pc_type) const {
            const auto &pc_type_name = backend_pc_types_.at(static_cast<size_t>(pc_type));
            assert(pc_type_name.empty() == false);
            return pc_type_name;
        }

        virtual void set_preconditioner(PreconditionerType pc_type, PreconditionerSide pc_side) = 0;

    private:
        static constexpr auto NUM_PRECONDITIONER_TYPES = static_cast<size_t>(PreconditionerType::NUM_TYPES);

        std::array<std::string, NUM_PRECONDITIONER_TYPES> backend_pc_types_;
    };

}  // namespace utopia

#endif  // UTOPIA_BACKEND_PRECONDITIONED_SOLVER_HPP