#ifndef UTOPIA_BLOCK_PRECONDITIONER_HPP
#define UTOPIA_BLOCK_PRECONDITIONER_HPP

#include "utopia_Base.hpp"
#include "utopia_Preconditioner.hpp"
#include "utopia_Traits.hpp"

#include <vector>

namespace utopia {

    template <class Vector>
    class BlockPreconditioner : public Preconditioner<Vector> {
    public:
        using SizeType = typename utopia::Traits<Vector>::SizeType;
        virtual ~BlockPreconditioner() = default;

        inline void set_index_set(const std::vector<SizeType> &index_set) { index_set_ = index_set; }

        inline const std::vector<SizeType> &index_set() const { return index_set_; }

    private:
        std::vector<SizeType> index_set_;
    };

    template <class Matrix, class Vector>
    class BlockSolvePreconditioner : public BlockPreconditioner<Vector> {
    public:
        using SizeType = typename utopia::Traits<Vector>::SizeType;
        typedef utopia::LinearSolver<Matrix, Vector> LinearSolverT;

        virtual ~BlockSolvePreconditioner() = default;

        inline void set_index_set(const std::vector<SizeType> &index_set) { index_set_ = index_set; }

        inline const std::vector<SizeType> &index_set() const { return index_set_; }

    private:
        std::vector<SizeType> index_set_;
        std::shared_ptr<LinearSolver<Matrix, Vector> > linear_solver_;
    };
}  // namespace utopia

#endif  // UTOPIA_BLOCK_PRECONDITIONER_HPP
