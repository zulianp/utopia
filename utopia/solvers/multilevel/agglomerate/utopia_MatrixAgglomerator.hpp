#ifndef UTOPIA_MATRIX_AGGLOMERATOR_HPP
#define UTOPIA_MATRIX_AGGLOMERATOR_HPP

#include "utopia_Clonable.hpp"
#include "utopia_Input.hpp"
#include "utopia_Traits.hpp"
#include "utopia_Transfer.hpp"

#include <memory>

namespace utopia {

    template <class Matrix>
    class MatrixAgglomerator : public Clonable, public Configurable {
    public:
        using Transfer = utopia::Transfer<Matrix, typename Traits<Matrix>::Vector>;

        virtual ~MatrixAgglomerator() = default;
        virtual std::shared_ptr<Transfer> create_transfer(const Matrix &in) = 0;
        virtual std::shared_ptr<Transfer> create_truncated_transfer(const Matrix &in) { return create_transfer(in); }
        virtual void read(Input &) override {}
        MatrixAgglomerator *clone() const override = 0;
    };

}  // namespace utopia

#endif  // UTOPIA_MATRIX_AGGLOMERATOR_HPP
