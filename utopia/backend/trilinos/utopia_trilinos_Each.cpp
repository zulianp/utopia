#include "utopia_trilinos_Each_impl.hpp"

#include <Trilinos_version.h>

#if (TRILINOS_MAJOR_MINOR_VERSION >= 130100 && UTOPIA_REMOVE_TRILINOS_DEPRECATED_CODE)
#include <Tpetra_Access.hpp>
#endif

namespace utopia {

    void TpetraVectorEach::apply_read(const TpetraVector &v, std::function<void(const Scalar &)> &fun) {
        const auto &impl = raw_type(v);

#if (TRILINOS_MAJOR_MINOR_VERSION >= 130100 && UTOPIA_REMOVE_TRILINOS_DEPRECATED_CODE)
        auto view = impl->getLocalViewHost(Tpetra::Access::ReadWrite);
#else
        auto view = impl->getLocalViewHost();
#endif

        const auto r = range(v);

        For<>::apply(0, r.extent(), [&view, &fun](const std::size_t i) { fun(view(i, 0)); });
    }

    void TpetraMatrixEach::apply_read(const TpetraMatrix &mat, std::function<void(const Scalar &)> &fun) {
        // FIXME make performance version
        apply_read(mat, [&fun](const SizeType i, const SizeType j, const Scalar val) {
            UTOPIA_UNUSED(i);
            UTOPIA_UNUSED(j);

            fun(val);
        });
    }

    // template class Each<TpetraMatrix, 2, FillType::SPARSE>;
    template class Each<TpetraVector, 1, FillType::DENSE>;
    template class Each<TpetraVector, 1, FillType::DELEGATE>;
}  // namespace utopia
