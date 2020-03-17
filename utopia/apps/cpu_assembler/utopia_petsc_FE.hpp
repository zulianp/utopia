#ifndef UTOPIA_PETSC_FE_HPP
#define UTOPIA_PETSC_FE_HPP

#include "utopia_petsc_ForwardDeclarations.hpp"
#include "utopia_Box.hpp"
#include "utopia_ArrayView.hpp"
#include "utopia_VectorView.hpp"
#include "utopia_Quadrature.hpp"
#include "utopia_UniformQuad4.hpp"
#include "utopia_UniformHex8.hpp"
#include "utopia_petsc_Vector.hpp"

namespace utopia {

    class PetscUniformQuad4 : public UniformQuad4<PetscScalar> {};
    class PetscUniformHex8  : public UniformHex8<PetscScalar>  {};

    template<typename View, std::size_t N>
    class Accessor<std::array<VectorView<View>, N>> {
    public:
        using Matrix = std::array<VectorView<View>, N>;
        using T = typename Traits<View>::Scalar;

        static const T &get(const Matrix &t, const std::size_t &i, const std::size_t &j)
        {
            return t[i][j];
        }

        static void set(Matrix &t, const std::size_t &i, const std::size_t &j, const T &val)
        {
            t[i][j] = val;
        }
    };


    //FIXME these quadrature defs are redundant

    template<>
    class Quadrature<PetscUniformQuad4, 2, 2> {
    public:
        using Scalar   = PetscUniformQuad4::Scalar;
        // using SizeType = PetscUniformQuad4::SizeType;
        using Point    = PetscUniformQuad4::Point;
        using ViewDevice = Quadrature;
        using ViewHost   = Quadrature;

        static const int Order   = 2;
        static const int Dim     = 2;
        static const int NPoints = 6;

        inline static constexpr int n_points()
        {
            return NPoints;
        }

        inline static constexpr int dim()
        {
            return Dim;
        }

        void init()
        {
            Quad4Quadrature<Scalar, Order, Dim, NPoints>::get(points_, weights_);
        }

        template<class Point>
        inline void point(const int qp_idx, Point &p) const
        {
            p[0] = points_[qp_idx][0];
            p[1] = points_[qp_idx][1];
        }

        inline Point point(const int qp_idx) const
        {
            Point p;
            point(qp_idx, p);
            return p;
        }

        inline const Scalar &weight(const int qp_idx) const
        {
            return weights_[qp_idx];
        }

        Quadrature()
        {
            init();
        }

        Quadrature(const Quadrature &other)
        {
            for(int i = 0; i < NPoints; ++i) {
                points_[i][0] = other.points_[i][0];
                points_[i][1] = other.points_[i][1];
                weights_[i] = other.weights_[i];
            }
        }

        inline const ViewDevice &view_device() const
        {
            return *this;
        }

        inline const ViewHost &view_host() const
        {
            return *this;
        }

    private:
        std::array<Point, NPoints> points_;
        std::array<Scalar, NPoints> weights_;
    };

    template<>
    class Quadrature<PetscUniformHex8, 2, 3> {
    public:
        using Scalar   = PetscUniformHex8::Scalar;
        // using SizeType = PetscUniformHex8::SizeType;
        using Point    = PetscUniformHex8::Point;
        using ViewDevice = Quadrature;
        using ViewHost   = Quadrature;

        static const int Order   = 2;
        static const int Dim     = 3;
        static const int NPoints = 27;

        inline static constexpr int n_points()
        {
            return NPoints;
        }

        inline static constexpr int dim()
        {
            return Dim;
        }

        void init()
        {
            Hex8Quadrature<Scalar, Order, Dim, NPoints>::get(points_, weights_);
        }

        template<class Point>
        inline void point(const int qp_idx, Point &p) const
        {
            p[0] = points_[qp_idx][0];
            p[1] = points_[qp_idx][1];
            p[2] = points_[qp_idx][2];
        }

        inline Point point(const int qp_idx) const
        {
            Point p;
            point(qp_idx, p);
            return p;
        }

        inline const Scalar &weight(const int qp_idx) const
        {
            return weights_[qp_idx];
        }

        Quadrature()
        {
            init();
        }

        Quadrature(const Quadrature &other)
        {
            for(int i = 0; i < NPoints; ++i) {
                points_[i][0] = other.points_[i][0];
                points_[i][1] = other.points_[i][1];
                points_[i][2] = other.points_[i][2];
                weights_[i] = other.weights_[i];
            }
        }

        inline const ViewDevice &view_device() const
        {
            return *this;
        }

        inline const ViewHost &view_host() const
        {
            return *this;
        }

    private:
        std::array<Point, NPoints> points_;
        std::array<Scalar, NPoints> weights_;
    };


    template<typename Scalar_, int PhysicalDim>
    class Quadrature<Edge2<Scalar_, PhysicalDim>, 2, 1> {
    public:
        using Scalar   = Scalar_;
        // using SizeType = PetscUniformHex8::SizeType;
        using Point      = utopia::StaticVector<Scalar, 1>;
        using ViewDevice = Quadrature;
        using ViewHost   = Quadrature;

        static const int Order   = 2;
        static const int Dim     = 1;
        static const int NPoints = 12;

        inline static constexpr int n_points()
        {
            return NPoints;
        }

        inline static constexpr int dim()
        {
            return Dim;
        }

        void init()
        {
            utopia::Quadrature<Scalar, 6, 1>::get(points_, weights_);
        }

        template<class Point>
        inline void point(const int qp_idx, Point &p) const
        {
            p[0] = points_[qp_idx][0];
        }

        inline Point point(const int qp_idx) const
        {
            Point p;
            point(qp_idx, p);
            return p;
        }

        inline const Scalar &weight(const int qp_idx) const
        {
            return weights_[qp_idx];
        }

        Quadrature()
        {
            init();
        }

        Quadrature(const Quadrature &other)
        {
            for(int i = 0; i < NPoints; ++i) {
                points_[i][0] = other.points_[i][0];
                weights_[i] = other.weights_[i];
            }
        }

        inline const ViewDevice &view_device() const
        {
            return *this;
        }

        inline const ViewHost &view_host() const
        {
            return *this;
        }

    private:
        std::array<Point, NPoints> points_;
        std::array<Scalar, NPoints> weights_;
    };

}

#endif //UTOPIA_PETSC_FE_HPP
