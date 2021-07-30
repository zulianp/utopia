#ifndef UTOPIA_AABB_HPP
#define UTOPIA_AABB_HPP

namespace utopia {

    template <class Point_>
    class AABB {
    public:
        using Point = Point_;

        Point min;
        Point max;
    };

}  // namespace utopia

#endif  // UTOPIA_AABB_HPP
