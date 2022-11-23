#ifndef UTOPIA_LINEPOST_PROCESSOR_HPP
#define UTOPIA_LINEPOST_PROCESSOR_HPP

#include <memory>
#include "utopia_PostProcessor.hpp"

namespace utopia {

    template <class FunctionSpace, class Vector>
    class LinePostProcessor final : public PostProcessor<FunctionSpace, Vector> {
    public:
        LinePostProcessor();
        ~LinePostProcessor();
        void apply(FunctionSpace &V, const Vector &sol) override;
        void apply(FunctionSpace &V, const Vector &sol, const Vector &other) override;
        void describe(std::ostream &os = std::cout) const override;
        void export_values() const override;
        void read(Input &in) override;

    private:
        class Impl;
        std::unique_ptr<Impl> impl_;
    };

}  // namespace utopia

#endif  // UTOPIA_LINEPOST_PROCESSOR_HPP
