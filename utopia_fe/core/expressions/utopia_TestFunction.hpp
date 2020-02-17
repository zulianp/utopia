#ifndef UTOPIA_TEST_FUNCTION_HPP
#define UTOPIA_TEST_FUNCTION_HPP

#include "utopia_BasisFunction.hpp"
#include "utopia_FunctionSpace.hpp"
#include "utopia_FunctionalTraits.hpp"
#include <memory>

namespace utopia {
    template<class FunctionSpaceT>
    class TestFunction : public BasisFunction< TestFunction<FunctionSpaceT>, FunctionSpaceT > {
    public:
        static const int Order = utopia::FormTraits<FunctionSpaceT>::Order;
        typedef typename utopia::FormTraits<FunctionSpaceT>::Scalar Scalar;
        typedef typename utopia::FormTraits<FunctionSpaceT>::Implementation Implementation;

        typedef utopia::BasisFunction< TestFunction<FunctionSpaceT>, FunctionSpaceT > super;
        using super::super;

        inline std::string get_class() const override { return "TestFunction"; }
    };

    template<class Derived>
    inline TestFunction<Derived> test(FunctionSpace<Derived> &&space)
    {
        return TestFunction<Derived>(std::move<Derived>(space.derived()));
    }

    template<class Derived>
    inline TestFunction<Derived> test(const FunctionSpace<Derived> &space)
    {
        return TestFunction<Derived>(space.derived());
    }

    template<class Derived>
    inline TestFunction<Derived> test(const std::shared_ptr< FunctionSpace<Derived> > &space)
    {
        return TestFunction<Derived>(std::static_pointer_cast<Derived>(space));
    }

    template<class FunctionSpaceT>
    class Traits< TestFunction<FunctionSpaceT> > : public Traits<FunctionSpaceT> {};


    template<class Space, class AssemblyContext>
    class FunctionalTraits<TestFunction<Space>, AssemblyContext> : public FunctionalTraits<Space, AssemblyContext> {
    public:
        inline static int type(const TestFunction<Space> &expr, const AssemblyContext &ctx)  { return FunctionalTraits<Space, AssemblyContext>::type(*expr.space_ptr(), ctx);  }
        inline static int order(const TestFunction<Space> &expr, const AssemblyContext &ctx) { return FunctionalTraits<Space, AssemblyContext>::order(*expr.space_ptr(), ctx); }
    };
}

#endif //UTOPIA_TEST_FUNCTION_HPP
