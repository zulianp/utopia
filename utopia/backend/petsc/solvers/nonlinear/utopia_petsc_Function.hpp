#ifndef UTOPIA_PETSC_FUNCTION_HPP
#define UTOPIA_PETSC_FUNCTION_HPP

#include "utopia_Function.hpp"
#include "utopia_petsc_Types.hpp"

#include <memory>

namespace utopia {
    template<class Matrix, class Vector>
    class Function<Matrix, Vector, PETSC> {
    public:
        DEF_UTOPIA_SCALAR(Matrix);
        
        virtual ~Function() { }
        
        virtual bool value(const Vector &/*point*/, Scalar &/*value*/) const = 0;
        virtual bool gradient(const Vector &/*point*/, Vector &/*result*/) const = 0;
        virtual bool hessian(const Vector &x, Matrix &H) const = 0;
        
        virtual bool hessian(const Vector &/*point*/, Matrix &/*result*/, Matrix &/*preconditioner*/) const
        {
            return false;
        }
        
        virtual bool has_preconditioner() const {
            return false;
        }
        
        virtual bool update(const Vector &/*point*/) { return true; };
        
        class Data {
        public:
            Matrix H;
            Matrix H_pre;
            Vector g;
        };
        
        inline std::shared_ptr<Data> data() {
            if(!data_) {
                data_ = std::make_shared<Data>();
            }
            
            return data_;
        }
        
    private:
        std::shared_ptr<Data> data_;
    };
}


#endif  //UTOPIA_PETSC_FUNCTION_HPP
