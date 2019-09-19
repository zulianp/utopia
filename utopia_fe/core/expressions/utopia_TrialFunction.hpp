#ifndef UTOPIA_TRIAL_FUNCTION_HPP
#define UTOPIA_TRIAL_FUNCTION_HPP

#include "utopia_BasisFunction.hpp"
#include "utopia_FunctionSpace.hpp"
#include "utopia_Traverse.hpp"
#include "utopia_Traits.hpp"
#include "utopia_FunctionalTraits.hpp"

#include <memory>

namespace utopia {
    template<class FunctionSpaceT>
    class TrialFunction : public BasisFunction< TrialFunction<FunctionSpaceT>, FunctionSpaceT > {
    public:
        static const int Order = utopia::FormTraits<FunctionSpaceT>::Order;
        typedef typename utopia::FormTraits<FunctionSpaceT>::Scalar Scalar;
        typedef typename utopia::FormTraits<FunctionSpaceT>::Implementation Implementation;

        typedef utopia::BasisFunction< TrialFunction<FunctionSpaceT>, FunctionSpaceT > super;
        using super::super;

        inline std::string get_class() const override { return "TrialFunction"; }
    };

    template<class Derived>
    TrialFunction<Derived> trial(FunctionSpace<Derived> &&space)
    {
        return TrialFunction<Derived>(std::move<Derived>(space.derived()));
    }

    template<class Derived>
    TrialFunction<Derived> trial(const FunctionSpace<Derived> &space)
    {
        return TrialFunction<Derived>(space.derived());
    }

    template<class Derived>
    inline TrialFunction<Derived> trial(const std::shared_ptr< FunctionSpace<Derived> > &space)
    {
        return TrialFunction<Derived>(std::static_pointer_cast<Derived>(space));
    }

    template<class FunctionSpaceT>
    class Traits< TrialFunction<FunctionSpaceT> > : public Traits<FunctionSpaceT> {};


    template<class Space, class AssemblyContext>
    class FunctionalTraits<TrialFunction<Space>, AssemblyContext> : public FunctionalTraits<Space, AssemblyContext> {
    public:
        inline static int type(const TrialFunction<Space> &expr, const AssemblyContext &ctx)  { return FunctionalTraits<Space, AssemblyContext>::type(*expr.space_ptr(), ctx);  }
        inline static int order(const TrialFunction<Space> &expr, const AssemblyContext &ctx) { return FunctionalTraits<Space, AssemblyContext>::order(*expr.space_ptr(), ctx); }
    };
}

#endif //UTOPIA_TRIAL_FUNCTION_HPP
