#ifndef UTOPIA_SIDE_SETS_HPP
#define UTOPIA_SIDE_SETS_HPP

#include "utopia_ArrayView.hpp"
#include <algorithm>
#include <string>

namespace utopia {

    class SideSet {
    public:
        using BoundaryIdType = int;
        inline static constexpr BoundaryIdType left()    { return 1; }
        inline static constexpr BoundaryIdType right()   { return 2; }
        inline static constexpr BoundaryIdType bottom()  { return 3; }
        inline static constexpr BoundaryIdType top()     { return 4; }
        inline static constexpr BoundaryIdType front()   { return 5; }
        inline static constexpr BoundaryIdType back()    { return 6; }
        inline static constexpr BoundaryIdType invalid() { return 0; }

        inline static BoundaryIdType from_name(const std::string &str_in)
        {
            if(str_in.empty()) {
                return invalid();
            }

            std::string str;
            str.resize(str_in.size(), '0');

            std::transform(str_in.begin(), str_in.end(), str.begin(),
                [](unsigned char c){ return std::tolower(c);
            });

            if(str == "left") {
                return left();
            }

            if(str == "right") {
                return right();
            }

            if(str == "bottom") {
                return bottom();
            }

            if(str == "top") {
                return top();
            }

            if(str == "front") {
                return front();
            }

            if(str == "back") {
                return back();
            }

            return invalid();
        }
    };

    template<int Dim>
    class SideSets {};

    template<>
    class SideSets<1> {
    public:
        static const int n_sides = 2;

        using Sides = utopia::ArrayView<SideSet::BoundaryIdType, n_sides>;

        static constexpr Sides ids_ = {
            SideSet::left(), SideSet::right()
        };

        UTOPIA_INLINE_FUNCTION static constexpr const Sides &sides()
        {
            return ids_;
        }

        template<class TensorIndex, class Dims>
        UTOPIA_INLINE_FUNCTION static constexpr bool on_side(
            const SideSet::BoundaryIdType &boundary_id,
            const TensorIndex &tensor_index,
            const Dims &dims)
        {
            switch(boundary_id) {
                case SideSet::left():
                {
                    return tensor_index[0] == 0;
                }

                case SideSet::right():
                {
                    return tensor_index[0] == (dims[0] - 1);
                }

                default:
                {
                    return false;
                }
            }
        }
    };

    template<>
    class SideSets<2> {
    public:
        static const int n_sides = 4;

        using Sides = utopia::ArrayView<SideSet::BoundaryIdType, n_sides>;

        static constexpr Sides ids_ =  {
            SideSet::left(),
            SideSet::right(),
            SideSet::bottom(),
            SideSet::top()
        };

        UTOPIA_INLINE_FUNCTION static constexpr const Sides &sides()
        {
            return ids_;
        }

        template<class TensorIndex, class Dims>
        UTOPIA_INLINE_FUNCTION static constexpr bool on_side(
            const SideSet::BoundaryIdType &boundary_id,
            const TensorIndex &tensor_index,
            const Dims &dims)
        {
            switch(boundary_id) {
                case SideSet::left():
                {
                    return tensor_index[0] == 0;
                }

                case SideSet::right():
                {
                    return tensor_index[0] == (dims[0] - 1);
                }

                case SideSet::bottom():
                {
                    return tensor_index[1] == 0;
                }

                case SideSet::top():
                {
                    return tensor_index[1] == (dims[1] - 1);
                }

                default:
                {
                    return false;
                }
            }
        }
    };

    template<>
    class SideSets<3> {
    public:
        static const int n_sides = 6;

        using Sides = utopia::ArrayView<SideSet::BoundaryIdType, n_sides>;

        static constexpr Sides ids_ = {
            SideSet::left(),
            SideSet::right(),
            SideSet::bottom(),
            SideSet::top(),
            SideSet::front(),
            SideSet::back()
        };

        UTOPIA_INLINE_FUNCTION static constexpr const Sides &sides()
        {
            return ids_;
        }

        template<class TensorIndex, class Dims>
        UTOPIA_INLINE_FUNCTION static constexpr bool on_side(
            const SideSet::BoundaryIdType &boundary_id,
            const TensorIndex &tensor_index,
            const Dims &dims)
        {
            switch(boundary_id) {
                case SideSet::left():
                {
                    return tensor_index[0] == 0;
                }

                case SideSet::right():
                {
                    return tensor_index[0] == (dims[0] - 1);
                }

                case SideSet::bottom():
                {
                    return tensor_index[1] == 0;
                }

                case SideSet::top():
                {
                    return tensor_index[1] == (dims[1] - 1);
                }

                case SideSet::back():
                {
                    return tensor_index[2] == 0;
                }

                case SideSet::front():
                {
                    return tensor_index[2] == (dims[2] - 1);
                }

                default:
                {
                    return false;
                }
            }
        }
    };

    template<>
    class SideSets<DYNAMIC_SIZE> {
    public:
        template<class TensorIndex, class Dims>
        UTOPIA_INLINE_FUNCTION static bool on_side(
            const SideSet::BoundaryIdType &boundary_id,
            const TensorIndex &tensor_index,
            const Dims &dims)
        {
            return SideSets<3>::on_side(boundary_id, tensor_index, dims);
        }

    };

}

#endif //UTOPIA_SIDE_SETS_HPP
