#ifndef UTOPIA_INTEGRATOR_HPP
#define UTOPIA_INTEGRATOR_HPP

#include "utopia_Expression.hpp"
#include "utopia_Traits.hpp"
#include "utopia_FunctionalTraits.hpp"
#include "utopia_Traverse.hpp"
#include "utopia_Range.hpp"
#include "utopia_FEBackend.hpp"
#include "utopia_FindSpace.hpp"
#include "utopia_make_unique.hpp"

namespace utopia {

    template<class FunctionSpace>
    class Integrator {
    public:
        using AssemblyValues = typename Traits<FunctionSpace>::AssemblyValues;
        static const int BACKEND_TAG = Traits<FunctionSpace>::Backend;

        virtual ~Integrator() {}

        virtual Range test_range(const AssemblyContext<BACKEND_TAG> &ctx) const 
        {
            return FEBackend<BACKEND_TAG>::range(test(), ctx);
        }

        virtual int order(const AssemblyValues &values) const = 0;
        virtual void init_values(AssemblyValues &values) const = 0;
        virtual const FunctionSpace &test() const = 0; 
        virtual std::shared_ptr<FunctionSpace> test_ptr() const = 0;
        virtual bool is_surface() const = 0;
    };

    template<class FunctionSpace>
    class LinearIntegrator : public virtual Integrator<FunctionSpace>,
                             public Expression< LinearIntegrator<FunctionSpace> > {
    public:
        enum {
            StoreAs = UTOPIA_BY_REFERENCE
        };

        typedef utopia::Traits<FunctionSpace> TraitsT;
        typedef typename TraitsT::Vector ElementVector;

        using Type = ElementVector;

        static const int BACKEND_TAG = Traits<FunctionSpace>::Backend;

        virtual ~LinearIntegrator() {}
        
        virtual bool assemble(
            const Range &test_r,
            AssemblyContext<BACKEND_TAG> &ctx,
            ElementVector &result) = 0;

        bool assemble(
            AssemblyContext<BACKEND_TAG> &ctx,
            ElementVector &result)
        {
            auto test_r  = this->test_range(ctx);
            ctx.init_tensor(result, false);
            return assemble(test_r, ctx, result);
        }

        virtual std::string getClass() const override { return "LinearIntegrator"; }
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
        inline static int order(const LinearIntegrator<FunctionSpace> &expr, const AssemblyContext &ctx) { return expr.order(ctx); }
    };

    ////////////////////////////////////////////////////////////////////////////////////

    template<class FunctionSpace>
    class BilinearIntegrator : public virtual Integrator<FunctionSpace>,
                               public Expression< BilinearIntegrator<FunctionSpace> > {
    public:
        enum {
            StoreAs = UTOPIA_BY_REFERENCE
        };

        typedef utopia::Traits<FunctionSpace> TraitsT;
        typedef typename TraitsT::Matrix ElementMatrix;
        typedef typename TraitsT::Vector ElementVector;

        using Type = ElementMatrix;

        static const int BACKEND_TAG = Traits<FunctionSpace>::Backend;

        virtual ~BilinearIntegrator() {}
        
        virtual bool assemble(
            const Range &trial_r,
            const Range &test_r,
            AssemblyContext<BACKEND_TAG> &ctx,
            ElementMatrix &result) = 0;

        bool assemble(
            AssemblyContext<BACKEND_TAG> &ctx,
            ElementMatrix &result)
        {
            auto trial_r = this->trial_range(ctx);
            auto test_r  = this->test_range(ctx);
            ctx.init_tensor(result, false);
            return assemble(trial_r, test_r, ctx, result);
        }

        //symmetric by default
        virtual const FunctionSpace &trial() const { return this->test(); }
        virtual std::shared_ptr<FunctionSpace> trial_ptr() const { return this->test_ptr(); }

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
        inline static int type(const BilinearIntegrator<FunctionSpace> &expr,  const AssemblyContext &ctx) { return FunctionalTraits<FunctionSpace, AssemblyContext>::type(expr.test(), ctx);  }
        inline static int order(const BilinearIntegrator<FunctionSpace> &expr, const AssemblyContext &ctx) { return expr.order(ctx); }
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

    template<class FunctionSpace, class Visitor>
    inline static int traverse(const LinearIntegrator<FunctionSpace> &expr, Visitor &visitor)
    {
        return visitor.visit(expr);
    }

    template<class FunctionSpace, class Visitor>
    inline static int traverse(const BilinearIntegrator<FunctionSpace> &expr, Visitor &visitor)
    {
        return visitor.visit(expr);
    }

    //temporary fix
    template<class Space>
    std::shared_ptr<Space> find_any_space(const LinearIntegrator<Space> &expr)
    {
        return expr.test_ptr();
    }

    template<class Space>
    std::shared_ptr<Space> find_any_space(const LinearIntegrator<ProductFunctionSpace<Space>> &expr)
    {
        return expr.test_ptr()->subspace_ptr(0);
    }

    //temporary fix
    template<class Space>
    std::shared_ptr<Space> find_any_space(const BilinearIntegrator<Space> &expr)
    {
        return expr.test_ptr();
    }

    template<class Space>
    std::shared_ptr<Space> find_any_space(const BilinearIntegrator<ProductFunctionSpace<Space>> &expr)
    {
        return expr.test_ptr()->subspace_ptr(0);
    }

}

#endif //UTOPIA_INTEGRATOR_HPP
