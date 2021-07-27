#ifndef UTOPIA_STK_FE_ASSEMBLER_HPP
#define UTOPIA_STK_FE_ASSEMBLER_HPP

#include "utopia_stk_FunctionSpace.hpp"

namespace utopia {

    template <>
    class FEAssembler<utopia::stk::FunctionSpace> : public Configurable {
    public:
        using FunctionSpace = utopia::stk::FunctionSpace;
        using Matrix = Traits<FunctionSpace>::Matrix;
        using Vector = Traits<FunctionSpace>::Vector;
        using Environment = utopia::Environment<FunctionSpace>;

        virtual ~FEAssembler() = default;
        virtual void clear() {}
        virtual std::string name() const = 0;
        virtual bool assemble(const Vector &x, Matrix &hessian, Vector &gradient) = 0;

        inline void read(Input &in) override {
            if (env_) {
                std::string space_name;
                in.get("space", space_name);

                if (!space_name.empty()) {
                    assert(!space_);

                    space_ = env_->find_space(space_name);
                }
            }
        }

        inline void set_environment(const std::shared_ptr<Environment> &env) { env_ = env; }
        inline const std::shared_ptr<Environment> &environment() const { return env_; }

        inline void set_space(const std::shared_ptr<FunctionSpace> &space) { space_ = space; }
        inline std::shared_ptr<FunctionSpace> space() { return space_; }

    private:
        std::shared_ptr<Environment> env_;
        std::shared_ptr<FunctionSpace> space_;
    };

}  // namespace utopia

#endif  // UTOPIA_STK_FE_ASSEMBLER_HPP