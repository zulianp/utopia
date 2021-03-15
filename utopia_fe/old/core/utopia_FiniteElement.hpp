#ifndef UTOPIA_FINITE_ELEMENT_HPP
#define UTOPIA_FINITE_ELEMENT_HPP

#include <vector>

#include "utopia_FormTraits.hpp"
#include "utopia_TrialFunction.hpp"

namespace utopia {
    template <class Space, int Backend = Traits<Space>::Backend>
    class FiniteElement {
    public:
        using SizeType = std::size_t;  // typename Traits<Space>::SizeType;

        enum { StoreAs = UTOPIA_BY_REFERENCE };

        void set(const SizeType &i) {
            // DO SOMETHING
        }

        FiniteElement(Space &space) : space_(space) {}

        inline const Space &space() const { return space_; }

        inline Space &space() { return space_; }

    private:
        Space &space_;

        FiniteElement(const FiniteElement &other) : space_(other.space_) {}
    };

    template <class Space>
    class Traits<FiniteElement<Space>> : public Traits<Space> {};

    template <class Space>
    class TrialFunction<FiniteElement<Space>>
        : public BasisFunction<TrialFunction<FiniteElement<Space>>, FiniteElement<Space>> {
    public:
        static const int Order = utopia::FormTraits<Space>::Order;
        typedef typename utopia::FormTraits<Space>::Scalar Scalar;
        typedef typename utopia::FormTraits<Space>::Implementation Implementation;

        typedef utopia::BasisFunction<TrialFunction<FiniteElement<Space>>, FiniteElement<Space>> Super;
        using Super::Super;

        inline std::string get_class() const override { return "TrialFunction"; }
    };

    template <class Space>
    TrialFunction<FiniteElement<Space>> trial(FiniteElement<Space> &space) {
        return TrialFunction<FiniteElement<Space>>(make_ref(space));
    }

    template <class Space, class AssemblyContext>
    class FunctionalTraits<FiniteElement<Space>, AssemblyContext> {
    public:
        inline static int type(const FiniteElement<Space> &expr, const AssemblyContext &ctx) {
            return FunctionalTraits<Space, AssemblyContext>::type(expr.space(), ctx);
        }
        inline static int order(const FiniteElement<Space> &expr, const AssemblyContext &ctx) {
            return FunctionalTraits<Space, AssemblyContext>::order(expr.space(), ctx);
        }
    };

}  // namespace utopia

#endif  // UTOPIA_FINITE_ELEMENT_HPP
