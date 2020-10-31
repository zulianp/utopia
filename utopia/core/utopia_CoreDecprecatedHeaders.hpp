#ifndef UTOPIA_CORE_DECPRECATED_HEADERS_HPP
#define UTOPIA_CORE_DECPRECATED_HEADERS_HPP

#include "utopia_Base.hpp"

#ifdef UTOPIA_DEPRECATED_API
#include "utopia_Each.hpp"
#include "utopia_Eval_Factory.hpp"
#include "utopia_Factory.hpp"
#include "utopia_ParallelEach.hpp"

namespace utopia {

    /**    @defgroup factory Factory
     *      @brief  Factory methods allow for creating basic tensor in an easy way
     */

    /** \addtogroup Global
     * @brief Manipulation with objects on global adress space.
     * @ingroup factory
     *  @{
     */

    /// Returns identity matrix  \f$ I^{row \times cols}  \f$.
    UTOPIA_DEPRECATED_MSG("Use the object variant with 1st argument as the communicator.")
    inline Factory<Identity, 2> identity(const Size::SizeType rows, const Size::SizeType cols) {
        return Factory<Identity, 2>(Size({rows, cols}));
    }
    /// Returns identity matrix  \f$ I^{size_0 \times size_1}  \f$.
    UTOPIA_DEPRECATED_MSG("Use the object variant with 1st argument as the communicator.")
    inline Factory<Identity, 2> identity(const Size &size) { return Factory<Identity, 2>(size); }

    /// Returns identity matrix  \f$ I^{row \times cols}  \f$.
    UTOPIA_DEPRECATED_MSG("Use the object variant with 1st argument as the communicator.")
    inline Factory<DenseIdentity, 2> dense_identity(const Size::SizeType rows, const Size::SizeType cols) {
        return Factory<DenseIdentity, 2>(Size({rows, cols}));
    }
    /// Returns denDensese_identity matrix  \f$ I^{size_0 \times size_1}  \f$.
    UTOPIA_DEPRECATED_MSG("Use the object variant with 1st argument as the communicator.")
    inline Factory<DenseIdentity, 2> dense_identity(const Size &size) { return Factory<DenseIdentity, 2>(size); }

    /// Returns identity matrix  \f$ I^{row \times cols}  \f$.
    inline constexpr SymbolicTensor<Identity, 2> identity() { return {}; }

    /// Returns global \f$ 0^{rows \times rows}  \f$.
    UTOPIA_DEPRECATED_MSG("Use the object variant with 1st argument as the communicator.")
    inline Factory<Zeros, 2> zeros(const Size::SizeType rows, const Size::SizeType cols) {
        return Factory<Zeros, 2>(Size({rows, cols}));
    }

    /// Returns global \f$ 0^{n \times 1}  \f$.
    UTOPIA_DEPRECATED_MSG("Use the object variant with 1st argument as the communicator.")
    inline Factory<Zeros, 1> zeros(const Size::SizeType n) { return Factory<Zeros, 1>(Size({n})); }

    /// Returns global \f$ 0^{size \times size}  \f$.
    UTOPIA_DEPRECATED_MSG("Use the object variant with 1st argument as the communicator.")
    inline Factory<Zeros, utopia::DYNAMIC> zeros(const Size &size) { return Factory<Zeros, utopia::DYNAMIC>(size); }

    /// nnz_x_row_or_col depends if your using a row-major or col-major sparse storage

    template <typename T>
    UTOPIA_DEPRECATED_MSG("Use the object variant with 1st argument as the communicator.")
    inline Factory<NNZ<T>, 2> sparse(const Size::SizeType rows, const Size::SizeType cols, T nnz_x_row_or_col) {
        return Factory<NNZ<T>, 2>(Size({rows, cols}), NNZ<T>(nnz_x_row_or_col));
    }

    template <typename SizeType>
    UTOPIA_DEPRECATED_MSG("Use the object variant with 1st argument as the communicator.")
    inline Factory<NNZXRow<SizeType>, 2> sparse(const Size &gs,
                                                const std::vector<SizeType> &d_nnz,
                                                const std::vector<SizeType> &o_nnz) {
        return Factory<NNZXRow<SizeType>, 2>(gs, NNZXRow<SizeType>(d_nnz, o_nnz));
    }

    ///  Returns global matrix \f$ value * 1^{rows \times cols}  \f$.
    template <typename T>
    UTOPIA_DEPRECATED_MSG("Use the object variant with 1st argument as the communicator.")
    inline Factory<Values<T>, 2> values(const Size::SizeType rows, const Size::SizeType cols, T value) {
        return Factory<Values<T>, 2>(Size({rows, cols}), Values<T>(value));
    }

    ///  Returns global vector \f$ value *1I^{rows \times 1}  \f$.

    template <typename T>
    UTOPIA_DEPRECATED_MSG("Use the object variant with 1st argument as the communicator.")
    inline Factory<Values<T>, 1> values(const Size::SizeType rows, T value) {
        return Factory<Values<T>, 1>(Size({rows, 1}), Values<T>(value));
    }
    /** @}*/

    /** \addtogroup Local
     * @brief Manipulation with objects on local adress space.
     * @ingroup factory
     *  @{
     */

    /// Returns local identity matrix  \f$ I^{row \times cols}  \f$ i.e. each processors owns local identity matrix.
    UTOPIA_DEPRECATED_MSG("Use the object variant with 1st argument as the communicator.")
    inline Factory<LocalIdentity, 2> local_identity(const Size::SizeType rows, const Size::SizeType cols) {
        return Factory<LocalIdentity, 2>(Size({rows, cols}));
    }

    /// Returns local identity matrix  \f$ I^{size \times size}  \f$ i.e. each processors owns local identity matrix.
    UTOPIA_DEPRECATED_MSG("Use the object variant with 1st argument as the communicator.")
    inline Factory<LocalIdentity, 2> local_identity(const Size &size) { return Factory<LocalIdentity, 2>(size); }

    /// Returns local identity matrix  \f$ I^{size \times size}  \f$ i.e. each processors owns local identity matrix.
    UTOPIA_DEPRECATED_MSG("Use the object variant with 1st argument as the communicator.")
    inline Factory<LocalDenseIdentity, 2> local_dense_identity(const Size &size) {
        return Factory<LocalDenseIdentity, 2>(size);
    }

    ///  Returns local zero vector \f$ 0^{n \times 1}  \f$.
    UTOPIA_DEPRECATED_MSG("Use the object variant with 1st argument as the communicator.")
    inline Factory<LocalZeros, 1> local_zeros(const Size::SizeType n) { return Factory<LocalZeros, 1>(Size({n})); }

    ///  Returns local zero matrices \f$ 0^{size \times size}  \f$ i.e. each processors owns local zero matrix.
    UTOPIA_DEPRECATED_MSG("Use the object variant with 1st argument as the communicator.")
    inline Factory<LocalZeros, utopia::DYNAMIC> local_zeros(const Size &size) {
        return Factory<LocalZeros, utopia::DYNAMIC>(size);
    }

    ///  Returns global matricx \f$ ?^{size_0 \times size_1}  \f$.
    UTOPIA_DEPRECATED_MSG("Use the object variant with 1st argument as the communicator.")
    inline Factory<Resize, utopia::DYNAMIC> dense(const Size &size) { return Factory<Resize, utopia::DYNAMIC>(size); }

    ///  Returns global matricx \f$ ?^{size \times size}  \f$.
    UTOPIA_DEPRECATED_MSG("Use the object variant with 1st argument as the communicator.")
    inline Factory<Resize, 2> dense(const Size::SizeType rows, const Size::SizeType cols) {
        return Factory<Resize, 2>(Size({rows, cols}));
    }

    ///  Returns local matrices \f$ value * 1^{rows \times cols}  \f$ i.e. each processors owns local matrix.

    template <typename T>
    UTOPIA_DEPRECATED_MSG("Use the object variant with 1st argument as the communicator.")
    inline Factory<LocalValues<T>, 2> local_values(const Size::SizeType rows, const Size::SizeType cols, T value) {
        return Factory<LocalValues<T>, 2>(Size({rows, cols}), LocalValues<T>(value));
    }

    ///  Returns local vector \f$ value * 1^{rows \times 1}  \f$ i.e. each processors owns local vector.

    template <typename T>
    UTOPIA_DEPRECATED_MSG("Use the object variant with 1st argument as the communicator.")
    inline Factory<LocalValues<T>, 1> local_values(const Size::SizeType rows, T value) {
        return Factory<LocalValues<T>, 1>(Size({rows, 1}), LocalValues<T>(value));
    }

    /// nnz_x_row_or_col depends if your using a row-major or col-major sparse storage

    template <typename T>
    UTOPIA_DEPRECATED_MSG("Use the object variant with 1st argument as the communicator.")
    inline Factory<LocalNNZ<T>, 2> local_sparse(const Size::SizeType rows,
                                                const Size::SizeType cols,
                                                T nnz_x_row_or_col) {
        return Factory<LocalNNZ<T>, 2>(Size({rows, cols}), LocalNNZ<T>(nnz_x_row_or_col));
    }

    template <typename T>
    UTOPIA_DEPRECATED_MSG("Use the object variant with 1st argument as the communicator.")
    inline Factory<LocalNNZ<T>, 2> local_sparse(const Size &s, T nnz_x_row_or_col) {
        return Factory<LocalNNZ<T>, 2>(s, LocalNNZ<T>(nnz_x_row_or_col));
    }

    template <typename T, class... Args>
    UTOPIA_DEPRECATED_MSG("Use the object variant with 1st argument as the communicator.")
    inline auto local_sparse(const Size::SizeType rows, const Size::SizeType cols, T nnz_x_row_or_col, Args &&... opts)
        -> Build<Factory<LocalNNZ<T>, 2>, decltype(options(opts...))> {
        return Build<Factory<LocalNNZ<T>, 2>, decltype(options(opts...))>(local_sparse(rows, cols, nnz_x_row_or_col),
                                                                          options(opts...));
    }

    template <class T, class Traits, int Backend>
    class Eval<Construct<Tensor<T, 2>, Diag<Factory<Zeros, 1>>>, Traits, Backend> {
    public:
        using Expr = utopia::Construct<Tensor<T, 2>, Diag<Factory<Zeros, 1>>>;

        inline static void apply(const Expr &expr) {
            UTOPIA_TRACE_BEGIN(expr);

            auto &l = Eval<Tensor<T, 2>, Traits>::apply(expr.left());
            l.zeros(matrix_layout(layout(expr.right())));

            UTOPIA_TRACE_END(expr);
        }
    };

    /** @}*/
}  // namespace utopia

#endif  // UTOPIA_DEPRECATED_API
#endif  // UTOPIA_CORE_DECPRECATED_HEADERS_HPP
