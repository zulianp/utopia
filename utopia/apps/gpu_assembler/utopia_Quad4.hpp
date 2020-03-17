#ifndef UTOPIA_QUAD_4_HPP
#define UTOPIA_QUAD_4_HPP

namespace utopia {

    template<class Scalar_, int PhysicalDim>
    class Quad4 {
    public:
        using Scalar = Scalar_;
        static const int Dim = PhysicalDim;
        static const int NNodes = 4;
        static const int NFunctions = 4;
        using Point = utopia::StaticVector<Scalar, Dim>;
        using GradValue = utopia::StaticVector<Scalar, Dim>;
        using STGradX   = utopia::StaticVector<Scalar, Dim-1>;
        using FunValue  = Scalar;
        using MemType   = utopia::Varying<>;
        using Side      = utopia::Node1<Scalar, PhysicalDim>;

        //IMPLEMENT ME

    };
}

#endif //UTOPIA_QUAD_4_HPP
