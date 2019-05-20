#ifndef UTOPIA_INTEGRATOR_HPP
#define UTOPIA_INTEGRATOR_HPP

#include "utopia_Expression.hpp"
#include "utopia_Traits.hpp"
#include "utopia_FunctionalTraits.hpp"
#include "utopia_Traverse.hpp"
#include "utopia_Range.hpp"
#include "utopia_FEBackend.hpp"

namespace utopia {

    template<class FunctionSpace>
    class Integrator {
    public:
        static const int BACKEND_TAG = Traits<FunctionSpace>::Backend;

        virtual ~Integrator() {}

        Integrator(FunctionSpace &space)
        : space_(space)
        {}

        inline const FunctionSpace &space() const
        {
            return space_;
        }

        virtual int order() const = 0;
        virtual void init_context(AssemblyContext<BACKEND_TAG> &ctx) = 0;
        
        
        virtual const FunctionSpace &test()  const { return space(); }
        
        virtual Range test_range(const AssemblyContext<BACKEND_TAG> &ctx) const 
        {
            return FEBackend<BACKEND_TAG>::range(test(), ctx);
        }

    private:
        FunctionSpace &space_;
    };

    template<class FunctionSpace>
    class LinearIntegrator : public Integrator<FunctionSpace>,
                             public Expression< LinearIntegrator<FunctionSpace> > {
    public:
        enum {
            StoreAs = UTOPIA_BY_REFERENCE
        };

        using ElementVector = typename Traits<FunctionSpace>::ElementVector;
        using Type = ElementVector;

        static const int BACKEND_TAG = Traits<FunctionSpace>::Backend;

        virtual ~LinearIntegrator() {}
        
        virtual bool assemble(
            const Range &test_r,
            const AssemblyContext<BACKEND_TAG> &ctx,
            ElementVector &result) = 0;

        bool assemble(
            const AssemblyContext<BACKEND_TAG> &ctx,
            ElementVector &result)
        {
            auto test_r  = test_range(ctx);
            return assemble(test_r, ctx, result);
        }
    };

    template<class FunctionSpace>
    class Traits< LinearIntegrator<FunctionSpace> > : public Traits<FunctionSpace> {
    public:
        static const int FILL_TYPE = utopia::FillType::DENSE;
        static const int Order = 1;
    };

    template<class FunctionSpace, class AssemblyContext>
    class FunctionalTraits<LinearIntegrator<FunctionSpace>, AssemblyContext>  {
    public:
        inline static int type(const LinearIntegrator<FunctionSpace> &expr,  const AssemblyContext &ctx) { return FunctionalTraits<FunctionSpace, AssemblyContext>::type(expr.space(), ctx);  }
        inline static int order(const LinearIntegrator<FunctionSpace> &expr, const AssemblyContext &ctx) { return expr.order(); }
    };

    ////////////////////////////////////////////////////////////////////////////////////

    template<class FunctionSpace>
    class BilinearIntegrator : public Integrator<FunctionSpace>,
                               public Expression< BilinearIntegrator<FunctionSpace> > {
    public:
        enum {
            StoreAs = UTOPIA_BY_REFERENCE
        };

        using ElementMatrix = typename Traits<FunctionSpace>::ElementMatrix;
        using Type = ElementMatrix;

        static const int BACKEND_TAG = Traits<FunctionSpace>::Backend;

        virtual ~BilinearIntegrator() {}
        
        virtual bool assemble(
            const Range &trial_r,
            const Range &test_r,
            const AssemblyContext<BACKEND_TAG> &ctx,
            ElementMatrix &result) = 0;

        bool assemble(
            const AssemblyContext<BACKEND_TAG> &ctx,
            ElementMatrix &result)
        {
            auto trial_r = trial_range(ctx);
            auto test_r  = test_range(ctx);

            return assemble(trial_r, test_r, ctx, result);
        }

        virtual const FunctionSpace &trial() const override { return this->space(); }

        virtual Range trial_range(const AssemblyContext<BACKEND_TAG> &ctx) const 
        {
            return FEBackend<BACKEND_TAG>::range(trial(), ctx);
        }
    };

    template<class FunctionSpace>
    class Traits< BilinearIntegrator<FunctionSpace> > : public Traits<FunctionSpace> {
    public:
        static const int FILL_TYPE = utopia::FillType::DENSE;
        static const int Order = 2;
    };

    template<class FunctionSpace, class AssemblyContext>
    class FunctionalTraits<BilinearIntegrator<FunctionSpace>, AssemblyContext>  {
    public:
        inline static int type(const BilinearIntegrator<FunctionSpace> &expr,  const AssemblyContext &ctx) { return FunctionalTraits<FunctionSpace, AssemblyContext>::type(expr.space(), ctx);  }
        inline static int order(const BilinearIntegrator<FunctionSpace> &expr, const AssemblyContext &ctx) { return expr.order(); }
    };

    //////////////////////////////////////////////////////////////////////

    template<class FunctionSpace, class Visitor>
    inline static int traverse(LinearIntegrator<FunctionSpace> &expr, Visitor &visitor)
    {
        return visitor.visit(expr);
    }

    template<class FunctionSpace, class Visitor>
    inline static int traverse(BilinearIntegrator<FunctionSpace> &expr, Visitor &visitor)
    {
        return visitor.visit(expr);
    }

}

#endif //UTOPIA_INTEGRATOR_HPP
