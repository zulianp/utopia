#ifndef UTOPIA_FE_FORWARD_DECLARATIONS_HPP
#define UTOPIA_FE_FORWARD_DECLARATIONS_HPP

#include "utopia_ForwardDeclarations.hpp"

namespace utopia {

    class FEExpression;
    class Mesh;

    template<class> class Coefficient;
    template<class> class Curl;
    template<class> class DifferentialOperator;
    template<class> class Divergence;
    template<class> class FormExpressions;
    template<class> class FunctionSpace;
    template<class> class Gradient;
    template<class> class Integral;
    template<class> class TestFunction;
    template<class> class TrialFunction;
    template<class> class BlockVar;
    template<class> class Dual;

    template<class, class> class Interpolate;
    template<class, int> class TimeDerivative;
    template<class> class ProductFunctionSpace;

    template<int> class AssemblyContext;

    template<class ...> class Equations;
    template<class ...> class FEConstraints;
    template<class, class> class ContextFunction;
    template<typename,int> class ConstantCoefficient;

    template<class Expr, class Out_, int Order_, class Fun>
    class Filter;

}

#endif //UTOPIA_FE_FORWARD_DECLARATIONS_HPP

