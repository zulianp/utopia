#ifndef UTOPIA_UI_EXTENDED_FUNCTION_HPP
#define UTOPIA_UI_EXTENDED_FUNCTION_HPP

#include "utopia_ExtendedFunction.hpp"
#include "utopia_Bratu.hpp"
#include "utopia_MinSurf.hpp"
#include "utopia_Poisson.hpp"

#include "utopia_ui.hpp"
#include <string>

namespace utopia {
    template<class FunctionSpace, class Matrix, class Vector>
    class UIExtendedFunction final : public ExtendedFunction<Matrix, Vector>, public Serializable {
    public:
        UIExtendedFunction(FunctionSpace &V) : V_(V) {}

        ~UIExtendedFunction() {}

        void read(InputStream &is) override {
            std::string name;
            is.read("name", name);
            initialize(name);
        }

        void initialize(const std::string &name)
        {
            if(name == "bratu") {
                fun_ = std::make_shared<Bratu<decltype(V_), USparseMatrix, UVector>>(V_);
            } else if(name == "min-surf") {
                fun_ = std::make_shared<MinSurf<decltype(V_), USparseMatrix, UVector>>(V_);
            } else
            // if(name == "poisson") {
            // 	fun_ = std::make_shared<Poisson<decltype(V_), USparseMatrix, UVector>>(V_);
            // } else
            {
                assert(false);
            }
        }

        inline bool value(const Vector &x, typename Vector::Scalar &energy) const override
        {
            return fun_->value(x, energy);
        }

        inline bool gradient(const Vector &x, Vector &gradient) const override
        {
              return fun_->gradient(x, gradient);
        }

        inline bool hessian(const Vector &x, Matrix &hessian) const override
        {
            return fun_->hessian(x, hessian);
        }

    private:
        FunctionSpace &V_;
        std::shared_ptr<ExtendedFunction<Matrix, Vector>> fun_;
    };
}


#endif //UTOPIA_UI_EXTENDED_FUNCTION_HPP
