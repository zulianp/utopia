//////////////////////////////////////////////////////////////////////////////////
// Copyright (c) 2015, Patrick Zulian - patrick.zulian@usi.ch                   //
// Institute of Computational Science - USI Università della Svizzera Italiana  //
//                                                                              //
// Redistribution and use in source and binary forms, with or without           //
// modification, are permitted provided that the following conditions are met:  //
//                                                                              //
// 1. Redistributions of source code must retain the above copyright notice,    //
//    this list of conditions and the following disclaimer.                     //
//                                                                              //
// 2. Redistributions in binary form must reproduce the above copyright notice, //
//    this list of conditions and the following disclaimer in the documentation //
//    and/or other materials provided with the distribution.                    //
//                                                                              //
// 3. Neither the name of the copyright holder nor the names of its contributors//
//    may be used to endorse or promote products derived from this software     //
//    without specific prior written permission.                                //
//                                                                              //
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"  //
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE    //
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE   //
// ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE    //
// LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR          //
// CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF         //
// SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS     //
// INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN      //
// CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)      //
// ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE   //
// POSSIBILITY OF SUCH DAMAGE.                                                  //
//////////////////////////////////////////////////////////////////////////////////

// Created by Patrick Zulian on 15/05/15.

#ifndef utopia_utopia_BASE_HPP
#define utopia_utopia_BASE_HPP

#include "utopia_Config.hpp"

#include <cassert>
#include <iomanip>
#include <iostream>
#include <sstream>

// namespace std {
//    template<typename _Tp>
//    constexpr _Tp &&
//    forward(typename std::remove_reference<_Tp>::type &__t) noexcept { return static_cast<_Tp &&>(__t); }

//    template<typename _Tp>
//    constexpr _Tp &&
//    forward(typename std::remove_reference<_Tp>::type &&__t) noexcept {
//        static_assert(!std::is_lvalue_reference<_Tp>::value, "template argument substituting _Tp is an lvalue
//        reference type"); return static_cast<_Tp &&>(__t);
//    }
//}

namespace utopia {
    using SizeType = long;
    static const SizeType INVALID_INDEX = -1;
    static const int DYNAMIC = -1;
    static const int UTOPIA_MAX_TENSOR_ORDER = 2;

    template <class Iter>
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

        for (auto c : clazz) {
            switch (c) {
                case '<': {
                    indent += space;
                    ss << "\n" << indent;
                    break;
                }

                case '>': {
                    indent.resize(indent.size() - space.size());
                    break;
                }

                case ',': {
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
    inline std::string treeFormat(const std::string &clazz) { return tree_format(clazz); }

    template <class TensorType>
    class is_parallel {
    public:
        enum { value = false };
    };

    template <typename T>
    void assert_equal(const T &left,
                      const T &right,
                      const std::string &filename,
                      const int line,
                      const std::string &expr_string) {
        if (left != right) {
            std::cerr << "assertion failure: " << expr_string << " with " << left << " != " << right << std::endl;
            std::cerr << "at " << filename << ":" << line << std::endl;
            abort();
        }
    }
}  // namespace utopia

#define CONST_DERIVED_CRT(Derived) \
    inline const Derived &derived() const { return static_cast<const Derived &>(*this); }
#define DERIVED_CRT(Derived) \
    inline Derived &derived() { return static_cast<Derived &>(*this); }

#define ALL_DERIVED_CRT(Derived) \
    CONST_DERIVED_CRT(Derived)   \
    DERIVED_CRT(Derived)

#define DEF_UTOPIA_SCALAR(Tensor) typedef typename utopia::Traits<Tensor>::Scalar Scalar
#define UTOPIA_SCALAR(Tensor) typename utopia::Traits<Tensor>::Scalar
#define UTOPIA_SIZE_TYPE(Tensor) typename utopia::Traits<Tensor>::SizeType
#define UTOPIA_MATRIX(Tensor) typename utopia::Traits<Tensor>::Matrix
#define UTOPIA_VECTOR(Tensor) typename utopia::Traits<Tensor>::Vector
#define UTOPIA_UNUSED(macro_var_) (void)macro_var_
#define utopia_assert(macro_expr__) assert((macro_expr__))

#ifndef NDEBUG
#define utopia_test_assert(macro_expr_) assert((macro_expr_))
#define utopia_assert_equal(macro_left_, macro_right_) \
    utopia::assert_equal(macro_left_, macro_right_, __FILE__, __LINE__, #macro_left_ " == " #macro_right_)
#else
namespace utopia {
    void test_check_assertion(const bool expr,
                              const std::string &filename,
                              const int line,
                              const std::string &expr_string);
}  // namespace utopia
#define utopia_test_assert(macro_expr_) utopia::test_check_assertion(macro_expr_, __FILE__, __LINE__, #macro_expr_)
#define utopia_assert_equal(...)
#endif  // NDEBUG

#endif  // utopia_utopia_BASE_HPP
