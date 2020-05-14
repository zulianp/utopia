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
#include "utopia_Quad4.hpp"
#include "utopia_Tri3.hpp"

namespace utopia {

    using PetscUniformQuad4 = utopia::UniformQuad4<PetscScalar>;
    using PetscUniformHex8  = utopia::UniformHex8<PetscScalar>;

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

    template<int Order_, int NPoints_, int Dim_>
    class QuadratureBase {
    public:
        static const int Order   = Order_;
        static const int Dim     = Dim_;
        static const int NPoints = NPoints_;

        using Scalar   = PetscScalar;
        using Point    = utopia::StaticVector<PetscScalar, Dim_>;

        inline static constexpr int n_points()
        {
            return NPoints;
        }

        inline static constexpr int dim()
        {
            return Dim;
        }

        template<class Point>
        inline void point(const int qp_idx, Point &p) const
        {
            for(int d = 0; d < Dim; ++d) {
                p[d] = points_[qp_idx][d];
            }
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

        void copy(const QuadratureBase &other)
        {
            for(int i = 0; i < NPoints; ++i) {
                points_[i].copy(other.points_[i]);
                weights_[i] = other.weights_[i];
            }
        }

        QuadratureBase() = default;

        QuadratureBase(const QuadratureBase &other)
        {
            copy(other);
        }

        void describe(std::ostream &os = std::cout) const
        {
            for(auto w : weights_) {
                os << w << " ";
            }

            os << std::endl;
        }


    protected:
        std::array<Point, NPoints> points_;
        std::array<Scalar, NPoints> weights_{};
    };

    //FIXME these quadrature defs are redundant

    template<>
    class Quadrature<PetscUniformQuad4, 2, 2> : public QuadratureBase<2, 6, 2> {
    public:
        using Super      = utopia::QuadratureBase<2, 6, 2>;
        using Scalar     = Super::Scalar;
        using Point      = Super::Point;
        using ViewDevice = Quadrature;
        using ViewHost   = Quadrature;

        Quadrature()
        {
            init();
        }

        Quadrature(const Quadrature &other)

            = default;

        void init()
        {
            Quad4Quadrature<Scalar, Super::Order, Super::Dim, Super::NPoints>::get(this->points_, this->weights_);
        }

        inline const ViewDevice &view_device() const
        {
            return *this;
        }

        inline const ViewHost &view_host() const
        {
            return *this;
        }
    };

    template<>
    class Quadrature<PetscUniformQuad4, 0, 2> : public QuadratureBase<0, 4, 2> {
    public:
        using Super      = utopia::QuadratureBase<0, 4, 2>;
        using Scalar     = Super::Scalar;
        using Point      = Super::Point;
        using ViewDevice = Quadrature;
        using ViewHost   = Quadrature;

        Quadrature()
        {
            init();
        }

        Quadrature(const Quadrature &other)

            = default;

        void init()
        {
            this->points_[0][0] = 0.211324865405187;
            this->points_[0][1] = 0.211324865405187;
            this->points_[1][0] = 0.788675134594813;
            this->points_[1][1] = 0.211324865405187;
            this->points_[2][0] = 0.211324865405187;
            this->points_[2][1] = 0.788675134594813;
            this->points_[3][0] = 0.788675134594813;
            this->points_[3][1] = 0.788675134594813;

            this->weights_ = {
                0.25,
                0.25,
                0.25,
                0.25
            };
        }

        inline const ViewDevice &view_device() const
        {
            return *this;
        }

        inline const ViewHost &view_host() const
        {
            return *this;
        }
    };

    template<typename Scalar_, int PhysicalDim, int Order>
    class Quadrature<Quad4<Scalar_, PhysicalDim>, Order, 2> : public Quadrature<PetscUniformQuad4, Order, 2>  {
    public:
        using ViewDevice = utopia::Quadrature<PetscUniformQuad4, Order, 2>;
        using ViewHost   = utopia::Quadrature<PetscUniformQuad4, Order, 2>;
    };

    template<>
    class Quadrature<PetscUniformHex8, 0, 3> : public QuadratureBase<0, 6, 3>{
    public:
        using Super      = utopia::QuadratureBase<0, 6, 3>;
        using Scalar     = Super::Scalar;
        using Point      = Super::Point;
        using ViewDevice = Quadrature;
        using ViewHost   = Quadrature;

        void init()
        {
            this->points_[0][0] = 0.0;
            this->points_[0][1] = 0.5;
            this->points_[0][2] = 0.5;
            this->points_[1][0] = 0.5;
            this->points_[1][1] = 0.0;
            this->points_[1][2] = 0.5;
            this->points_[2][0] = 0.5;
            this->points_[2][1] = 0.5;
            this->points_[2][2] = 0.0;
            this->points_[3][0] = 0.5;
            this->points_[3][1] = 0.5;
            this->points_[3][2] = 1.0;
            this->points_[4][0] = 0.5;
            this->points_[4][1] = 1.0;
            this->points_[4][2] = 0.5;
            this->points_[5][0] = 1.0;
            this->points_[5][1] = 0.5;
            this->points_[5][2] = 0.5;

            this->weights_ = {
                0.16666666666666666666666666666667,
                0.16666666666666666666666666666667,
                0.16666666666666666666666666666667,
                0.16666666666666666666666666666667,
                0.16666666666666666666666666666667,
                0.16666666666666666666666666666667
            };
        }

        Quadrature()
        {
            init();
        }

        Quadrature(const Quadrature &other)

            = default;

        inline const ViewDevice &view_device() const
        {
            return *this;
        }

        inline const ViewHost &view_host() const
        {
            return *this;
        }
    };

    template<>
    class Quadrature<PetscUniformHex8, 2, 3> : public QuadratureBase<2, 27, 3>{
    public:
        using Super      = utopia::QuadratureBase<2, 27, 3>;
        using Scalar     = Super::Scalar;
        using Point      = Super::Point;
        using ViewDevice = Quadrature;
        using ViewHost   = Quadrature;

        void init()
        {
            Hex8Quadrature<Scalar, Super::Order, Super::Dim, Super::NPoints>::get(this->points_, this->weights_);
        }

        Quadrature()
        {
            init();
        }

        Quadrature(const Quadrature &other)

            = default;

        inline const ViewDevice &view_device() const
        {
            return *this;
        }

        inline const ViewHost &view_host() const
        {
            return *this;
        }
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

        void describe(std::ostream &os = std::cout) const
        {
            for(auto w : weights_) {
                os << w << " ";
            }

            os << std::endl;
        }

    private:
        std::array<Point, NPoints> points_;
        std::array<Scalar, NPoints> weights_;
    };

    template<typename Scalar_, int PhysicalDim>
    class Quadrature<Tri3<Scalar_, PhysicalDim>, 2, 2> : public QuadratureBase<2, 6, 2> {
    public:
        using Super      = utopia::QuadratureBase<2, 6, 2>;
        using Scalar     = Super::Scalar;
        using Point      = Super::Point;
        using ViewDevice = Quadrature;
        using ViewHost   = Quadrature;

        void init()
        {
            this->points_[0][0] = 0.5;
            this->points_[0][1] = 0.5;
            this->points_[1][0] = 0.5;
            this->points_[1][1] = 0.0;
            this->points_[2][0] = 0.0;
            this->points_[2][1] = 0.5;
            this->points_[3][0] = 1.0/6.0;
            this->points_[3][1] = 1.0/6.0;
            this->points_[4][0] = 1.0/6.0;
            this->points_[4][1] = 2.0/3.0;
            this->points_[5][0] = 2.0/3.0;
            this->points_[5][1] = 1.0/6.0;

            this->weights_ = { 1.0/30.0, 1.0/30.0, 1.0/30.0, 0.3, 0.3, 0.3 };
        }

        Quadrature()
        {
            init();
        }

        Quadrature(const Quadrature &other)
        : Super(other)
        {}

        inline const ViewDevice &view_device() const
        {
            return *this;
        }

        inline const ViewHost &view_host() const
        {
            return *this;
        }
    };

    template<typename Scalar_, int PhysicalDim>
    class Quadrature<Tri3<Scalar_, PhysicalDim>, 0, 2> : public QuadratureBase<0, 1, 2> {
    public:
        using Super      = utopia::QuadratureBase<0, 1, 2>;
        using Scalar     = Super::Scalar;
        using Point      = Super::Point;
        using ViewDevice = Quadrature;
        using ViewHost   = Quadrature;

        void init()
        {
            this->points_[0][0] = 0.;
            this->points_[0][1] = 0.;

            this->weights_ = { 1.0 };
        }

        Quadrature()
        {
            init();
        }

        Quadrature(const Quadrature &other)
        : Super(other)
        {}

        inline const ViewDevice &view_device() const
        {
            return *this;
        }

        inline const ViewHost &view_host() const
        {
            return *this;
        }
    };

}

#endif //UTOPIA_PETSC_FE_HPP
