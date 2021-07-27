#ifndef UTOPIA_MOONOLITH_HOUSEHOLDER_REFLECTION_HPP
#define UTOPIA_MOONOLITH_HOUSEHOLDER_REFLECTION_HPP

namespace utopia {

    template <class Matrix, class Vector, int Dim>
    class HouseholderReflectionForContact {
    public:
        static void build(const Vector &is_contact, const Vector &normal, Matrix &trafo);
    };

}  // namespace utopia

#endif  // UTOPIA_MOONOLITH_HOUSEHOLDER_REFLECTION_HPP
