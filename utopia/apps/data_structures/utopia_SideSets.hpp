#ifndef UTOPIA_SIDE_SETS_HPP
#define UTOPIA_SIDE_SETS_HPP

#include "utopia_ArrayView.hpp"

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
    };

}

#endif //UTOPIA_SIDE_SETS_HPP
