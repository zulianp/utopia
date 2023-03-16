#ifndef UTOPIA_QUADRATURE_BASED_ASSEMBLER_HPP
#define UTOPIA_QUADRATURE_BASED_ASSEMBLER_HPP

#include <cassert>
#include <memory>

namespace utopia {

    class QMortarBuilder;

    class QuadratureBasedAssembler {
    public:
        QuadratureBasedAssembler();
        virtual ~QuadratureBasedAssembler();

        inline const QMortarBuilder &get_q_builder() const {
            assert(q_builder);
            return *q_builder;
        }

        void set_q_builder(const std::shared_ptr<QMortarBuilder> &q_builder);

    protected:
        std::shared_ptr<QMortarBuilder> q_builder;
    };

}  // namespace utopia

#endif  // UTOPIA_QUADRATURE_BASED_ASSEMBLER_HPP
