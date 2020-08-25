#ifndef UTOPIA_ML_EVAL_HPP
#define UTOPIA_ML_EVAL_HPP

#include "utopia_Core.hpp"
#include "utopia_ExtendedFunction.hpp"
#include "utopia_Function.hpp"
#include "utopia_LevelMemory.hpp"

namespace utopia {
enum LocalSolveType { PRE_SMOOTHING = 1, POST_SMOOTHING = 2, COARSE_SOLVE = 0 };

enum MultiLevelCoherence {
  FIRST_ORDER = 1,
  FIRST_ORDER_DF = 3,
  FIRST_ORDER_MGOPT = 4,
  FIRST_ORDER_MULTIPLICATIVE_DF = 6,
  FIRST_ORDER_ADDITIVE_MULTIPLICATIVE_DF = 7,

  SECOND_ORDER = 2,
  SECOND_ORDER_DF = 5,

  GALERKIN = 0
};

template <MultiLevelCoherence T, MultiLevelCoherence U>
struct is_same : std::false_type {};

template <MultiLevelCoherence T>
struct is_same<T, T> : std::true_type {};

template <MultiLevelCoherence T, MultiLevelCoherence... Rest>
struct is_any : std::false_type {};

template <MultiLevelCoherence T, MultiLevelCoherence First>
struct is_any<T, First> : is_same<T, First> {};

template <MultiLevelCoherence T, MultiLevelCoherence First,
          MultiLevelCoherence... Rest>
struct is_any<T, First, Rest...>
    : std::integral_constant<bool, is_same<T, First>::value ||
                                       is_any<T, Rest...>::value> {};

// helper function, should be implemented in c++ 14
template <bool B, class T = void>
using enable_if_t = typename std::enable_if<B, T>::type;

template <class Matrix, class Vector, MultiLevelCoherence CONSISTENCY_TYPE>
class MultilevelDerivEval {};

template <MultiLevelCoherence T>
struct is_first_order {
  static const bool value = false;
};

template <>
struct is_first_order<FIRST_ORDER> {
  static const bool value = true;
};

template <>
struct is_first_order<FIRST_ORDER_MGOPT> {
  static const bool value = true;
};

template <>
struct is_first_order<FIRST_ORDER_DF> {
  static const bool value = true;
};

template <>
struct is_first_order<FIRST_ORDER_MULTIPLICATIVE_DF> {
  static const bool value = true;
};

template <>
struct is_first_order<FIRST_ORDER_ADDITIVE_MULTIPLICATIVE_DF> {
  static const bool value = true;
};

template <MultiLevelCoherence T>
struct is_derivative_free {
  static const bool value = false;
};

template <>
struct is_derivative_free<FIRST_ORDER_DF> {
  static const bool value = true;
};

template <>
struct is_derivative_free<FIRST_ORDER_MULTIPLICATIVE_DF> {
  static const bool value = true;
};

template <>
struct is_derivative_free<FIRST_ORDER_ADDITIVE_MULTIPLICATIVE_DF> {
  static const bool value = true;
};

template <>
struct is_derivative_free<SECOND_ORDER_DF> {
  static const bool value = true;
};

}  // namespace utopia

#endif  // UTOPIA_ML_EVAL_HPP