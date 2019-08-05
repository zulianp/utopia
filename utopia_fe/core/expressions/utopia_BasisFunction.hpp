#ifndef UTOPIA_BASIS_FUNCTION_HPP
#define UTOPIA_BASIS_FUNCTION_HPP

#include "utopia_Expression.hpp"
#include "utopia_FEExpression.hpp"
#include "utopia_FEForwardDeclarations.hpp"

#include <memory>

namespace utopia {
    template<class Derived, class FunctionSpaceT>
    class BasisFunction : public Expression<Derived>, public FEExpression {
    public:
        BasisFunction(const FunctionSpaceT &space)
        : space_ptr_(std::make_shared<FunctionSpaceT>(space))
        {}

        BasisFunction(std::shared_ptr<FunctionSpaceT> &space_ptr)
        : space_ptr_(space_ptr)
        {}

        BasisFunction(FunctionSpaceT &&space)
        : space_ptr_(std::make_shared<FunctionSpaceT>(std::move(space)))
        {}

        inline std::shared_ptr<FunctionSpaceT> space_ptr() const
        {
            return space_ptr_;
        }

    private:
        std::shared_ptr<FunctionSpaceT> space_ptr_;
    };


    template<template<class> class Function, class SubSpace>
    class BasisFunction<Function< ProductFunctionSpace<SubSpace> >, ProductFunctionSpace<SubSpace> >
    : public Expression< Function< ProductFunctionSpace<SubSpace> > >, public FEExpression {
    public:
        typedef utopia::ProductFunctionSpace<SubSpace> FunctionSpaceT;

        BasisFunction(const FunctionSpaceT &space)
        : space_ptr_(std::make_shared<FunctionSpaceT>(space))
        {}

        BasisFunction(std::shared_ptr<FunctionSpaceT> &space_ptr)
        : space_ptr_(space_ptr)
        {}

        BasisFunction(FunctionSpaceT &&space)
        : space_ptr_(std::make_shared<FunctionSpaceT>(std::move(space)))
        {}

        inline std::shared_ptr<FunctionSpaceT> space_ptr() const
        {
            return space_ptr_;
        }

        inline Function<SubSpace> operator[](const unsigned int sub_space_index) const
        {
            return Function<SubSpace>(space_ptr()->subspace_ptr(sub_space_index));
        }
        std::size_t codim() const
        {
            return space_ptr()->n_subspaces();
        }

    private:
        std::shared_ptr<FunctionSpaceT> space_ptr_;
    };
}

#endif //UTOPIA_BASIS_FUNCTION_HPP
