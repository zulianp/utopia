//
// Created by Patrick Zulian on 15/05/15.
//

#ifndef utopia_utopia_BASE_HPP
#define utopia_utopia_BASE_HPP

#include "utopia_Config.hpp"

#include <assert.h>
#include <sstream>



//namespace std {
//    template<typename _Tp>
//    constexpr _Tp &&
//    forward(typename std::remove_reference<_Tp>::type &__t) noexcept { return static_cast<_Tp &&>(__t); }

//    template<typename _Tp>
//    constexpr _Tp &&
//    forward(typename std::remove_reference<_Tp>::type &&__t) noexcept {
//        static_assert(!std::is_lvalue_reference<_Tp>::value, "template argument substituting _Tp is an lvalue reference type");
//        return static_cast<_Tp &&>(__t);
//    }
//}

namespace utopia {
    typedef long int SizeType;
    static const SizeType INVALID_INDEX = -1;
    static const int DYNAMIC = -1;
    static const int UTOPIA_MAX_TENSOR_ORDER = 2;

    template<class Iter>
    void disp(const Iter &begin, const Iter &end, std::ostream &os) {
        for (Iter it = begin; it != end; ++it) {
            os << *it << " ";
        }
        os << "\n";
    }



    inline std::string tree_format(const std::string &clazz) {
        std::stringstream ss;

        std::string indent = "";
        std::string space = "|\t";


        for (auto c: clazz) {
            switch (c) {
                case '<' : {
                    indent += space;
                    ss << "\n" << indent;
                    break;
                }

                case '>' : {

                    indent.resize(indent.size() - space.size());
                    break;
                }

                case ',' : {
                    ss << "\n" << indent;
                    break;
                }

                case ' ': {
                    break;
                }

                default: {
                    ss << c;
                    break;
                }
            }
        }
        return ss.str();
    }


    //@decprecated use tree_format
   inline std::string treeFormat(const std::string &clazz) {
        return tree_format(clazz);
    }


    template<class TensorType>
    class is_parallel { public: enum { value = false }; };
}

#define CONST_DERIVED_CRT(Derived) inline const Derived &derived() const { return static_cast<const Derived &>(*this); }
#define DERIVED_CRT(Derived) inline Derived &derived() { return static_cast<Derived &>(*this); }

#define ALL_DERIVED_CRT(Derived)  \
CONST_DERIVED_CRT(Derived)        \
DERIVED_CRT(Derived)

#define DEF_UTOPIA_SCALAR(Tensor) typedef typename utopia::Traits<Tensor>::Scalar Scalar;
#define UTOPIA_SCALAR(Tensor) typename utopia::Traits<Tensor>::Scalar
#define UTOPIA_SIZE_TYPE(Tensor) typename utopia::Traits<Tensor>::SizeType
#define UTOPIA_MATRIX(Tensor) typename utopia::Traits<Tensor>::Matrix
#define UTOPIA_VECTOR(Tensor) typename utopia::Traits<Tensor>::Vector





#endif //utopia_utopia_BASE_HPP
