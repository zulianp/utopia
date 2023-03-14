#ifndef UTOPIA_BACKEND_PRECONDITIONED_SOLVER_HPP
#define UTOPIA_BACKEND_PRECONDITIONED_SOLVER_HPP

#include <array>
#include <cassert>
#include <map>
#include <string>

#define UTOPIA_PCNONE "none"
#define UTOPIA_PCASM "asm"
#define UTOPIA_PCBJACOBI "bjacobi"
#define UTOPIA_PCCHEBYSHEV "chebyshev"
#define UTOPIA_PCHYPRE "hypre"
#define UTOPIA_PCILU "ilu"
#define UTOPIA_PCILUT "ilut"
#define UTOPIA_PCJACOBI "jacobi"
#define UTOPIA_PCLU "lu"
#define UTOPIA_PCMULTIGRID "multigrid"
#define UTOPIA_PCREDUNDANT "redundant"
#define UTOPIA_PCSOR "sor"

namespace utopia {

    class BackendPreconditionedSolver {
    public:
        enum class PreconditionerSide { INVALID = 0, LEFT, RIGHT };
        enum class PreconditionerType {
            NONE = 0,
            ASM,
            BJACOBI,
            CHEBYSHEV,
            HYPRE,
            ILU,
            ILUT,
            JACOBI,
            LU,
            MULTIGRID,
            REDUNDANT,
            SOR,
            NUM_TYPES
        };

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
            assert(utopia_pc_types.size() == NUM_PRECONDITIONER_TYPES);

            PreconditionerType pc_type_{PreconditionerType::NONE};
            if (!pc_type.empty()) {
                auto pc_type_lc{pc_type};
                std::transform(
                    pc_type.cbegin(), pc_type.cend(), pc_type_lc.begin(), [](char c) { return std::tolower(c); });
                auto pc_type_it = utopia_pc_types.find(pc_type_lc);
                // should not happen: string for preconditioner type not found
                assert(pc_type_it != utopia_pc_types.end());
                pc_type_ = pc_type_it->second;
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

        std::string get_preconditioner_name(PreconditionerType pc_type) const {
            auto utopia_pc_type_it = std::find_if(utopia_pc_types.begin(),
                                                  utopia_pc_types.end(),
                                                  [pc_type](const auto p) { return p.second == pc_type; });
            assert(utopia_pc_type_it != utopia_pc_types.end());
            return utopia_pc_type_it->first;
        }

        virtual std::string get_preconditioner_name() const = 0;
        virtual void set_preconditioner(PreconditionerType pc_type, PreconditionerSide pc_side) = 0;

    private:
        // map backend-independent name string to enum value
        const std::map<std::string, PreconditionerType> utopia_pc_types = {
            {UTOPIA_PCNONE, PreconditionerType::NONE},
            {UTOPIA_PCASM, PreconditionerType::ASM},
            {UTOPIA_PCBJACOBI, PreconditionerType::BJACOBI},
            {UTOPIA_PCCHEBYSHEV, PreconditionerType::CHEBYSHEV},
            {UTOPIA_PCHYPRE, PreconditionerType::HYPRE},
            {UTOPIA_PCILU, PreconditionerType::ILU},
            {UTOPIA_PCILUT, PreconditionerType::ILUT},
            {UTOPIA_PCJACOBI, PreconditionerType::JACOBI},
            {UTOPIA_PCLU, PreconditionerType::LU},
            {UTOPIA_PCMULTIGRID, PreconditionerType::MULTIGRID},
            {UTOPIA_PCREDUNDANT, PreconditionerType::REDUNDANT},
            {UTOPIA_PCSOR, PreconditionerType::SOR}};

        static constexpr auto NUM_PRECONDITIONER_TYPES = static_cast<size_t>(PreconditionerType::NUM_TYPES);

        std::array<std::string, NUM_PRECONDITIONER_TYPES> backend_pc_types_;
    };

}  // namespace utopia

#endif  // UTOPIA_BACKEND_PRECONDITIONED_SOLVER_HPP