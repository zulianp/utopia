#ifndef UTOPIA_TRILINOS_MAX_ROW_NNZ_HPP
#define UTOPIA_TRILINOS_MAX_ROW_NNZ_HPP

#include "utopia_Traits.hpp"

#include "utopia_MaxRowNNZ.hpp"

#include "utopia_Tpetra_Matrix.hpp"

namespace utopia {

    template <>
    class MaxRowNNZ<TpetraMatrix, TRILINOS> {
    public:
        using SizeType = Traits<TpetraMatrix>::SizeType;
        using LocalSizeType = Traits<TpetraMatrix>::LocalSizeType;
        using Scalar = Traits<TpetraMatrix>::Scalar;

        static SizeType apply(const TpetraMatrix &in);
    };

}  // namespace utopia

#endif  // UTOPIA_TRILINOS_MAX_ROW_NNZ_HPP
