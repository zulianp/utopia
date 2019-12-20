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

    template<int Dim_>
    class PetscElem {
    public:
        static const int Dim = Dim_;
        using SizeType = PetscInt;
        using Scalar   = PetscScalar;
        using Point    = utopia::StaticVector<Scalar, Dim>;
        using Grad     = utopia::StaticVector<Scalar, Dim>;
        using NodeIndex = utopia::ArrayView<const SizeType>;

        virtual ~PetscElem() {}
        PetscElem()
        {}

        inline void set(const NodeIndex &nodes)
        {
            nodes_ = nodes;
        }

        inline const NodeIndex &nodes() const
        {
            return nodes_;
        }


        inline NodeIndex &nodes()
        {
            return nodes_;
        }

        inline SizeType node_id(const SizeType k) const
        {
            return nodes_[k];
        }

        inline const SizeType &node(const SizeType &i) const
        {
            return nodes_[i];
        }

        inline SizeType n_nodes() const
        {
            return nodes_.size();
        }

        inline constexpr static SizeType dim() { return Dim; }
        virtual Scalar fun(const SizeType &i, const Point &p) const = 0;

        inline void idx(const SizeType &idx)
        {
            idx_ = idx;
        }

        inline SizeType idx() const
        {
            assert(idx_ >= 0);
            return idx_;
        }

    private:
        NodeIndex nodes_;
        SizeType idx_;
    };

    class PetscUniformQuad4 final : public PetscElem<2> {
    public:
        using Super    = utopia::PetscElem<2>;
        using SizeType = Super::SizeType;
        using Scalar   = Super::Scalar;
        using Point    = Super::Point;
        using Grad     = Super::Grad;
        using MemType  = Uniform<>;
        static const int Dim = 2;
        static const int NNodes = 4;

        inline Scalar fun(const SizeType &i, const Point &p) const
        {
            return impl_.fun(i, p);
        }

        inline void node(const SizeType &i, Point &p) const
        {
            return impl_.node(i, p);
        }

        inline void point(const Point &in, Point &out) const
        {
            impl_.point(in, out);
        }

        inline void centroid(Point &out) const
        {
            impl_.centroid(out);
        }

        inline void grad(const int i, const Point &p, Grad &g) const
        {
           impl_.grad(i, p, g);
        }

        inline constexpr static bool is_affine()
        {
            return UniformQuad4<Scalar>::is_affine();
        }

        inline constexpr static Scalar reference_measure()
        {
            return UniformQuad4<Scalar>::reference_measure();
        }

        inline Scalar measure() const
        {
            return impl_.measure();
        }

        inline constexpr static int n_nodes()
        {
            return UniformQuad4<Scalar>::n_nodes();
        }

        inline void set(
            const StaticVector2<Scalar> &translation,
            const StaticVector2<Scalar> &h)
        {
            impl_.set(translation, h);
        }

        void describe(std::ostream &os = std::cout) const
        {
            auto &t = impl_.translation();
            os << t(0) << " " << t(1) << "\n";
        }

    private:
        UniformQuad4<Scalar> impl_;
    };

    class PetscUniformHex8 final : public PetscElem<3> {
    public:
        using Super    = utopia::PetscElem<3>;
        using SizeType = Super::SizeType;
        using Scalar   = Super::Scalar;
        using Point    = Super::Point;
        using Grad     = Super::Grad;
        using MemType  = Uniform<>;
        static const int Dim = 3;
        static const int NNodes = 8;

        inline Scalar fun(const SizeType &i, const Point &p) const
        {
            return impl_.fun(i, p);
        }

        inline void node(const SizeType &i, Point &p) const
        {
            return impl_.node(i, p);
        }

        inline void point(const Point &in, Point &out) const
        {
            impl_.point(in, out);
        }

        inline void centroid(Point &out) const
        {
            impl_.centroid(out);
        }

        inline void grad(const int i, const Point &p, Grad &g) const
        {
           impl_.grad(i, p, g);
        }

        inline constexpr static bool is_affine()
        {
            return UniformHex8<Scalar>::is_affine();
        }

        inline constexpr static Scalar reference_measure()
        {
            return UniformHex8<Scalar>::reference_measure();
        }

        inline Scalar measure() const
        {
            return impl_.measure();
        }

        inline constexpr static int n_nodes()
        {
            return UniformHex8<Scalar>::n_nodes();
        }

        inline void set(
            const StaticVector3<Scalar> &translation,
            const StaticVector3<Scalar> &h)
        {
            impl_.set(translation, h);
        }

        void describe(std::ostream &os = std::cout) const
        {
            auto &t = impl_.translation();
            os << t(0) << " " << t(1) << " " << t(2) << "\n";
        }

    private:
        UniformHex8<Scalar> impl_;
    };

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

    template<>
    class Quadrature<PetscUniformQuad4, 2, 2> {
    public:
        using Scalar   = PetscUniformQuad4::Scalar;
        using SizeType = PetscUniformQuad4::SizeType;
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
        using SizeType = PetscUniformHex8::SizeType;
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

}

#endif //UTOPIA_PETSC_FE_HPP
