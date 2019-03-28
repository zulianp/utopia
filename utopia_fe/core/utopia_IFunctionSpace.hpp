#ifndef UTOPIA_I_FUNCTION_SPACE_HPP
#define UTOPIA_I_FUNCTION_SPACE_HPP

#include <memory>

namespace utopia {

    template<class FunctionSpaceT>
    class FunctionSpacePolicyChecker {
    public:
        using ElemT    = utopia::Traits<FunctionSpaceT>::Elem;
        using ScalarT  = utopia::Traits<FunctionSpaceT>::Scalar;
        using IntegerT = utopia::Traits<FunctionSpaceT>::Integer;

        static void dofs(const FunctionSpaceT &space, std::vector<Integer> &indices)
        {
            space.dofs(indices);
        }

        static void dofs(const FunctionSpaceT &space, const Integer var_num, std::vector<Integer> &indices)
        {
            space.dofs(var_num, indices);
        }

        // Integer n_variables() const = 0;
        // Integer n_dofs() const = 0;
        // Integer n_local_dofs() const = 0;

        // Integer n_elements() const = 0;
        // Elem & elem(const Integer id) = 0;
        // const Elem & elem(const Integer id) const = 0;
    };
}

#endif //UTOPIA_I_FUNCTION_SPACE_HPP
