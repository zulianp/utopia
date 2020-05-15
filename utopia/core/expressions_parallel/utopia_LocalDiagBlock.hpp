/*
 * @Author: alenakopanicakova
 * @Date:   2016-05-02
 * @Last Modified by:   alenakopanicakova
 * @Last Modified time: 2016-07-29
 */
#ifndef UTOPIA_UTOPIA_LOCAL_DIAG_BLOCK_HPP
#define UTOPIA_UTOPIA_LOCAL_DIAG_BLOCK_HPP

#include <string>
#include "utopia_Expression.hpp"
#include "utopia_StoreAs.hpp"

namespace utopia {
    template <class Expr>
    class LocalDiagBlock : public Expression<LocalDiagBlock<Expr> > {
    public:
        LocalDiagBlock(const Expr &expr) : expr_(expr) {}

        inline const Expr &expr() const { return expr_; }

        inline std::string get_class() const override { return "LocalDiagBlock<" + expr_.get_class() + ">"; }

    private:
        const Expr &expr_;
    };

    /*!
     * @ingroup  parallel_expressions
     * @brief  Extracts square diagonal blocks from parallel distributed matrix and assigns them to the local matrices.
     * \n The ordering of blocks corresponds to the ordering of the  MPI_rank.
     *
     */
    template <class Derived>
    LocalDiagBlock<Derived> local_diag_block(const Expression<Derived> &expr) {
        return LocalDiagBlock<Derived>(expr.derived());
    }
}  // namespace utopia

#endif  // UTOPIA_UTOPIA_LOCAL_DIAG_BLOCK_HPP
