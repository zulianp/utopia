#ifndef UTOPIA_MEASURE_RESIDUAL_HPP
#define UTOPIA_MEASURE_RESIDUAL_HPP

#include "utopia_Traits.hpp"

#include "utopia_Norm.hpp"

#include "utopia_Eval_Reduce.hpp"

#include <memory>

namespace utopia {
    template <class Vector>
    class MeasureResidual : public Configurable {
    public:
        using Scalar = typename Traits<Vector>::Scalar;

        virtual ~MeasureResidual() = default;
        void read(Input &) override {}

        virtual Scalar measure(const Vector &r) const { return norm2(r); }
    };

    // TODO a front-end version
    template <class Vector>
    class MeasureResidualComponents {};

    template <class Vector>
    class MeasureResidualFactory {
    public:
        static std::shared_ptr<MeasureResidual<Vector>> make(const std::string &) {
            return std::make_shared<MeasureResidual<Vector>>();
        }
    };

}  // namespace utopia

#endif  // UTOPIA_MEASURE_RESIDUAL_HPP
