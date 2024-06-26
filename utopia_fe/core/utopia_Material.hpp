#ifndef UTOPIA_MATERIAL_HPP
#define UTOPIA_MATERIAL_HPP

#include "utopia_Function.hpp"
#include "utopia_Operator.hpp"
#include "utopia_Traits.hpp"

namespace utopia {

    template <class FunctionSpace>
    class AbstractMaterial
        : public utopia::Function<typename Traits<FunctionSpace>::Matrix, typename Traits<FunctionSpace>::Vector>,
          public Operator<typename Traits<FunctionSpace>::Vector>,
          public Describable {
    public:
        virtual ~AbstractMaterial() = default;
        virtual void initialize(const std::shared_ptr<FunctionSpace> &space) = 0;

        virtual std::string name() const = 0;

        virtual bool has_hessian() const = 0;
        virtual bool has_gradient() const = 0;
        virtual bool has_value() const = 0;
        virtual bool is_operator() const = 0;
        virtual bool is_linear() const { return false; }
        bool is_hessian_constant() const override { return this->is_linear(); }
        virtual void must_apply_constraints(const bool) {}

        void describe(std::ostream &os = std::cout) const override {
            os << "AbstractMaterial::describe(...) Implement me!!!" << std::endl;
        }
    };

    template <class FunctionSpace, class FE>
    class Material {};

}  // namespace utopia

#endif  // UTOPIA_MATERIAL_HPP
