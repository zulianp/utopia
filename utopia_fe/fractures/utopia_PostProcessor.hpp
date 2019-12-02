#ifndef UTOPIA_POSTPROCESSOR_HPP
#define UTOPIA_POSTPROCESSOR_HPP

#include "utopia_Input.hpp"

namespace utopia {

    template<class FunctionSpace, class Vector>
    class PostProcessor : public Configurable {
    public:
        virtual ~PostProcessor() {}
        virtual void apply(FunctionSpace &V, const Vector &sol) = 0;
        virtual void apply(FunctionSpace &V, const Vector &sol,  const Vector &other) = 0;
        virtual void describe(std::ostream &os = std::cout) const = 0;
        virtual void export_values() const = 0;
    };

}

#endif //UTOPIA_POSTPROCESSOR_HPP
