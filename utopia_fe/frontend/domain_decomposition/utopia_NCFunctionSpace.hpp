#ifndef UTOPIA_NC_FUNCTION_SPACE_HPP
#define UTOPIA_NC_FUNCTION_SPACE_HPP

#include "utopia_Traits.hpp"
#include "utopia_Transfer.hpp"

#include <memory>

namespace utopia {

    template <class FunctionSpace>
    class NCFunctionSpace {
    public:
        using Matrix_t = typename Traits<FunctionSpace>::Matrix;
        using Vector_t = typename Traits<FunctionSpace>::Vector;
        using Scalar_t = typename Traits<FunctionSpace>::Scalar;
        using Size_t = typename Traits<FunctionSpace>::SizeType;
        using Transfer_t = utopia::Transfer<Matrix_t, Vector_t>;

        std::shared_ptr<FunctionSpace> space() const { return space_; }

        NCFunctionSpace(const std::shared_ptr<FunctionSpace> &space) : space_(space) {}

    private:
        std::shared_ptr<FunctionSpace> space_;
        std::shared_ptr<Transfer_t> projector_;
    };

    template <class FunctionSpace>
    class Traits<NCFunctionSpace<FunctionSpace>> : public Traits<FunctionSpace> {};

}  // namespace utopia

#endif  // UTOPIA_NC_FUNCTION_SPACE_HPP