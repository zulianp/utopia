#ifndef UTOPIA_TRILINOS_MAT_CHOP_HPP
#define UTOPIA_TRILINOS_MAT_CHOP_HPP

#include "utopia_MatChop.hpp"
#include "utopia_Tpetra_Matrix.hpp"
#include "utopia_Traits.hpp"

namespace utopia {

    template <>
    class Chop<TpetraMatrix, TRILINOS> {
        using Scalar = typename utopia::Traits<TpetraMatrix>::Scalar;

    public:
        static void apply(TpetraMatrix &mat, const Scalar &eps);
    };

    template <>
    class ChopSmallerThan<TpetraMatrix, TRILINOS> {
        using Scalar = typename utopia::Traits<TpetraMatrix>::Scalar;

    public:
        static void apply(TpetraMatrix &mat, const Scalar &eps);
    };

    template <>
    class ChopGreaterThan<TpetraMatrix, TRILINOS> {
        using Scalar = typename utopia::Traits<TpetraMatrix>::Scalar;

    public:
        static void apply(TpetraMatrix &mat, const Scalar &eps);
    };

}  // namespace utopia

#endif  // UTOPIA_TRILINOS_MAT_CHOP_HPP
