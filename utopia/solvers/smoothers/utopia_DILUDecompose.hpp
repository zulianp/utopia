#ifndef UTOPIA_DILU_DECOMPOSE_HPP
#define UTOPIA_DILU_DECOMPOSE_HPP

#include "utopia_CRSMatrix.hpp"
#include "utopia_ILUDecompose.hpp"
#include "utopia_Traits.hpp"

#include <vector>

namespace utopia {

    template <class Matrix, class Vector>
    class DILUAlgorithm final : public ILUAlgorithm<Matrix, Vector> {
    public:
        using Scalar = typename Traits<Matrix>::Scalar;

        static bool decompose(const Matrix &in, std::vector<Scalar> &d);

        bool update(const Matrix &mat) override;
        void apply(const Vector &b, Vector &x) override;
        void read(Input &) override;

    private:
        std::shared_ptr<const Matrix> mat_;
        std::vector<Scalar> d_;
        std::vector<SizeType> diag_idx_;
        std::vector<SizeType> temp_;
    };

}  // namespace utopia

#endif  // UTOPIA_DILU_DECOMPOSE_HPP
