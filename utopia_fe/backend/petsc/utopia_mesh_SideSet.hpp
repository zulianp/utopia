#ifndef UTOPIA_PETSC_SIDE_SETS_HPP
#define UTOPIA_PETSC_SIDE_SETS_HPP

#include <algorithm>
#include <string>
#include "utopia_ArrayView.hpp"

namespace utopia {
    namespace mesh {

        class SideSet {
        public:
            using BoundaryIdType = int;

            class Cube {
            public:
                inline static constexpr BoundaryIdType left_id() { return 1; }
                inline static constexpr BoundaryIdType right_id() { return 2; }
                inline static constexpr BoundaryIdType bottom_id() { return 3; }
                inline static constexpr BoundaryIdType top_id() { return 4; }
                inline static constexpr BoundaryIdType front_id() { return 5; }
                inline static constexpr BoundaryIdType back_id() { return 6; }
                inline static constexpr BoundaryIdType invalid_id() { return -1; }

                inline static std::string left() { return "left"; }
                inline static std::string right() { return "right"; }
                inline static std::string bottom() { return "bottom"; }
                inline static std::string top() { return "top"; }
                inline static std::string front() { return "front"; }
                inline static std::string back() { return "back"; }

                inline static BoundaryIdType id_from_user_space_name(const std::string &str_in) {
                    return convert(str_in);
                }

                inline static std::string name_from_id(const BoundaryIdType id) {
                    switch (id) {
                        case left_id():
                            return left();
                        case right_id():
                            return right();
                        case bottom_id():
                            return bottom();
                        case top_id():
                            return top();
                        case front_id():
                            return front();
                        case back_id():
                            return back();
                        default: {
                            return "Undefined";
                        }
                    }
                }

                inline static BoundaryIdType convert(const std::string &str_in) {
                    if (str_in.empty()) {
                        return invalid_id();
                    }

                    std::string str;
                    str.resize(str_in.size(), '0');

                    std::transform(
                        str_in.begin(), str_in.end(), str.begin(), [](unsigned char c) { return std::tolower(c); });

                    if (str == left()) {
                        return left_id();
                    }

                    if (str == right()) {
                        return right_id();
                    }

                    if (str == bottom()) {
                        return bottom_id();
                    }

                    if (str == top()) {
                        return top_id();
                    }

                    if (str == front()) {
                        return front_id();
                    }

                    if (str == back()) {
                        return back_id();
                    }

                    return invalid_id();
                }
            };
        };

        template <int Dim>
        class SideSets {};

        template <>
        class SideSets<1> {
        public:
            static const int n_sides = 2;

            using Sides = utopia::ArrayView<SideSet::BoundaryIdType, n_sides>;

            static constexpr Sides ids_ = {SideSet::Cube::left_id(), SideSet::Cube::right_id()};

            UTOPIA_INLINE_FUNCTION static constexpr const Sides &sides() { return ids_; }

            template <class TensorIndex, class Dims>
            UTOPIA_INLINE_FUNCTION static constexpr bool on_side(const SideSet::BoundaryIdType &boundary_id,
                                                                 const TensorIndex &tensor_index,
                                                                 const Dims &dims) {
                switch (boundary_id) {
                    case SideSet::Cube::left_id(): {
                        return tensor_index[0] == 0;
                    }

                    case SideSet::Cube::right_id(): {
                        return tensor_index[0] == (dims[0] - 1);
                    }

                    default: {
                        return false;
                    }
                }
            }
        };

        template <>
        class SideSets<2> {
        public:
            static const int n_sides = 4;

            using Sides = utopia::ArrayView<SideSet::BoundaryIdType, n_sides>;

            static constexpr Sides ids_ = {SideSet::Cube::left_id(),
                                           SideSet::Cube::right_id(),
                                           SideSet::Cube::bottom_id(),
                                           SideSet::Cube::top_id()};

            UTOPIA_INLINE_FUNCTION static constexpr const Sides &sides() { return ids_; }

            template <class TensorIndex, class Dims>
            UTOPIA_INLINE_FUNCTION static constexpr bool on_side(const SideSet::BoundaryIdType &boundary_id,
                                                                 const TensorIndex &tensor_index,
                                                                 const Dims &dims) {
                switch (boundary_id) {
                    case SideSet::Cube::left_id(): {
                        return tensor_index[0] == 0;
                    }

                    case SideSet::Cube::right_id(): {
                        return tensor_index[0] == (dims[0] - 1);
                    }

                    case SideSet::Cube::bottom_id(): {
                        return tensor_index[1] == 0;
                    }

                    case SideSet::Cube::top_id(): {
                        return tensor_index[1] == (dims[1] - 1);
                    }

                    default: {
                        return false;
                    }
                }
            }
        };

        template <>
        class SideSets<3> {
        public:
            static const int n_sides = 6;

            using Sides = utopia::ArrayView<SideSet::BoundaryIdType, n_sides>;

            static constexpr Sides ids_ = {SideSet::Cube::left_id(),
                                           SideSet::Cube::right_id(),
                                           SideSet::Cube::bottom_id(),
                                           SideSet::Cube::top_id(),
                                           SideSet::Cube::front_id(),
                                           SideSet::Cube::back_id()};

            UTOPIA_INLINE_FUNCTION static constexpr const Sides &sides() { return ids_; }

            template <class TensorIndex, class Dims>
            UTOPIA_INLINE_FUNCTION static constexpr bool on_side(const SideSet::BoundaryIdType &boundary_id,
                                                                 const TensorIndex &tensor_index,
                                                                 const Dims &dims) {
                switch (boundary_id) {
                    case SideSet::Cube::left_id(): {
                        return tensor_index[0] == 0;
                    }

                    case SideSet::Cube::right_id(): {
                        return tensor_index[0] == (dims[0] - 1);
                    }

                    case SideSet::Cube::bottom_id(): {
                        return tensor_index[1] == 0;
                    }

                    case SideSet::Cube::top_id(): {
                        return tensor_index[1] == (dims[1] - 1);
                    }

                    case SideSet::Cube::back_id(): {
                        return tensor_index[2] == 0;
                    }

                    case SideSet::Cube::front_id(): {
                        return tensor_index[2] == (dims[2] - 1);
                    }

                    default: {
                        return false;
                    }
                }
            }
        };

        template <>
        class SideSets<DYNAMIC_SIZE> {
        public:
            template <class TensorIndex, class Dims>
            UTOPIA_INLINE_FUNCTION static bool on_side(const SideSet::BoundaryIdType &boundary_id,
                                                       const TensorIndex &tensor_index,
                                                       const Dims &dims) {
                return SideSets<3>::on_side(boundary_id, tensor_index, dims);
            }
        };

    }  // namespace mesh
}  // namespace utopia

#endif  // UTOPIA_PETSC_SIDE_SETS_HPP
