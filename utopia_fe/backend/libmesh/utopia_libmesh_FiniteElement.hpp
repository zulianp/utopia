#ifndef UTOPIA_LIBMESH_FINITE_ELEMENT_HPP
#define UTOPIA_LIBMESH_FINITE_ELEMENT_HPP

#include "utopia_FiniteElement.hpp"
#include "utopia_libmesh_FunctionSpace.hpp"
#include "utopia_libmesh_AssemblyContext.hpp"

namespace utopia {
    template<class Space>
    class FiniteElement<Space, LIBMESH_TAG> {
    public:
        using SizeType = std::size_t;// typename Traits<Space>::SizeType;
        
        enum {
            StoreAs = UTOPIA_BY_REFERENCE
        };

        void set(const SizeType &i) 
        {
            ctx_.set_current_element(i);
            ctx_.set_has_assembled(false);
        }

        template<class Expr>
        void init(const Expr &expr)
        {
            ctx_.init(expr);
        }

        template<class Expr>
        void reinit(const Expr &expr)
        {
            ctx_.reinit(expr);
        }


        AssemblyContext<LIBMESH_TAG> &ctx()
        {
            return ctx_;
        }

        const AssemblyContext<LIBMESH_TAG> &ctx() const
        {
            return ctx_;
        }

        inline const Space &space() const
        {
            return space_;
        }

        inline Space &space()
        {
            return space_;
        }

        FiniteElement(Space &space) : space_(space)
        {}

        inline static std::string get_class() 
        {
            return "FiniteElement";
        }

    private:
        Space &space_;
        AssemblyContext<LIBMESH_TAG> ctx_;


        FiniteElement(const FiniteElement &other)
        : space_(other.space_)
        {}
    };

}

#endif //UTOPIA_LIBMESH_FINITE_ELEMENT_HPP
