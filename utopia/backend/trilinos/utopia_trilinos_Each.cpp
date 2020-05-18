#include "utopia_trilinos_Each_impl.hpp"

namespace utopia {

    void TpetraVectorEach::apply_read(const TpetraVector &v, std::function<void(const Scalar &)> &fun) {
        const auto &impl = raw_type(v);
        auto view = impl->getLocalViewHost();

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

    template class Each<TpetraMatrix, 2, FillType::SPARSE>;
    template class Each<TpetraVector, 1, FillType::DENSE>;
    template class Each<TpetraVector, 1, FillType::DELEGATE>;
}  // namespace utopia
