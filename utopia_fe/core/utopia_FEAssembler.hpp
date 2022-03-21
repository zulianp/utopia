#ifndef UTOPIA_FE_ASSEMBLER_HPP
#define UTOPIA_FE_ASSEMBLER_HPP

#include "utopia_Input.hpp"
#include "utopia_Traits.hpp"

#include "utopia_SimulationTime.hpp"
#include "utopia_fe_Environment.hpp"

namespace utopia {

    template <class FunctionSpace_>
    class FEAssembler : public Configurable {
    public:
        using FunctionSpace = FunctionSpace_;
        using Matrix = typename Traits<FunctionSpace>::Matrix;
        using Vector = typename Traits<FunctionSpace>::Vector;
        using Environment = utopia::Environment<FunctionSpace>;
        using Scalar = typename Traits<Vector>::Scalar;
        using SimulationTime = utopia::SimulationTime<Scalar>;

        virtual ~FEAssembler() = default;
        virtual void clear() {}
        virtual std::string name() const = 0;
        virtual bool assemble(const Vector &x, Matrix &hessian, Vector &gradient) = 0;
        virtual bool assemble(const Vector &x, Matrix &jacobian) = 0;
        virtual bool assemble(const Vector &x, Vector &fun) = 0;
        virtual bool apply(const Vector &x, Vector &hessian_times_x) = 0;
        virtual bool is_operator() const = 0;

        // For linear only
        virtual bool assemble(Matrix &jacobian) = 0;
        virtual bool assemble(Vector &fun) = 0;

        inline void read(Input &in) override {
            if (environment() && !space()) {
                std::string space_name;
                in.get("space", space_name);

                if (!space_name.empty()) {
                    set_space(environment()->find_space(space_name));
                }
            }
        }

        virtual void set_environment(const std::shared_ptr<Environment> &env) = 0;
        virtual std::shared_ptr<Environment> environment() const = 0;

        virtual void set_space(const std::shared_ptr<FunctionSpace> &space) = 0;
        virtual std::shared_ptr<FunctionSpace> space() const = 0;

        virtual bool is_linear() const = 0;

        virtual void set_time(const std::shared_ptr<SimulationTime> &time) = 0;

        // private:
        //     std::shared_ptr<Environment> env_;
        //     std::shared_ptr<FunctionSpace> space_;
    };

}  // namespace utopia

#endif  // UTOPIA_FE_ASSEMBLER_HPP