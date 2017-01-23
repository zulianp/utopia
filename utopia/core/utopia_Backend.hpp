//
// Created by Patrick Zulian on 18/05/15.
//

#ifndef utopia_utopia_BACKEND_HPP
#define utopia_utopia_BACKEND_HPP

#include <vector>
#include <memory>
#include <iostream>

#include "utopia_Base.hpp"

namespace utopia {

    template<typename Scalar, int BackendType>
    class Backend {
    };

//    template<typename T>
//    class Backend<T, CUSTOM> {
//    public:
//        typedef std::vector<T> Vector;
//        typedef Vector Matrix;
//
//        static Backend &Instance() {
//            static Backend instance;
//            return instance;
//        }
//
//        template<class Fun>
//        void apply(const Vector &left, const Vector &right, Vector &result, const Fun &fun)
//        {
//            using std::move;
//            assert(left.size() == right.size());
//            result.resize(left.size());
//
//            for(typename Vector::size_type i = 0; i < left.size(); ++i) {
//                result[i] = fun(left[i], right[i]);
//            }
//        }
//
//        template<class Operation>
//        Vector apply(const Vector &left, const Vector &right, const Operation &)
//        {
//            using std::move;
//
//            Vector result;
//            apply(left, right, result, Operation::template Fun<T>());
//            return move(result);
//        }
//
//        template<class Operation>
//        inline T reduce(const Vector &v, const Operation &)
//        {
//            auto fun = Operation::template Fun<T>();
//            T result = 0;
//            for(typename Vector::size_type i = 0; i < v.size(); ++i) {
//                result = fun(result, v[i]);
//            }
//
//            return result;
//        }
//
//        template<class Operation>
//        Vector apply(const T scalar, const Vector &v, const Operation &)
//        {
//            using std::move;
//            typename Vector::size_type size = v.size();
//            std::vector<T> result(size);
//
//            auto fun = Operation::template Fun<T>();
//
//            for (typename Vector::size_type i = 0; i < size; ++i) {
//                result[i] = fun(scalar, v[i]);
//            }
//
//            return move(result);
//        }
//
//        template<class Operation>
//        Vector apply(const Vector &v, const T scalar, const Operation &op)
//        {
//            return apply(scalar, v, op);
//        }
//
//        Vector axpy(const T scaleFactor, const Vector &left, const Vector &right) {
//            std::cout << "axpy" << std::endl;
//            using std::move;
//            assert(left.size() == right.size());
//            typename Vector::size_type size = right.size();
//            std::vector<T> result(size);
//
//            for (typename Vector::size_type i = 0; i < size; ++i) {
//                result[i] = scaleFactor * left[i] + right[i];
//            }
//
//            return move(result);
//        }
//
//
//        void assign(Vector &left, Vector &&right)
//        {
//            left = right;
//        }
//    private:
//        Backend() { }
//    };
}
#endif //utopia_utopia_BACKEND_HPP
