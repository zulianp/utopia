#ifndef UTOPIA_LIBMESH_TRANSFORM_HPP
#define UTOPIA_LIBMESH_TRANSFORM_HPP

#include "moonolith_transform.hpp"
#include "libmesh/elem.h"
#include "libmesh/dense_matrix.h"

namespace utopia {

    class Transform {
    public:
        virtual ~Transform() {}
        virtual void transform_to_reference(const libMesh::Point &world, libMesh::Point &ref) const = 0;
        virtual void apply(const libMesh::Point &ref, libMesh::Point &world) const = 0;
        virtual void jacobian(const libMesh::Point &ref, libMesh::DenseMatrix<libMesh::Real> &result) const
        {
            (void) ref;
            (void) result;
            assert(false && "implement me");
        }
    };

    class Transform1 : 
        public Transform,
        public moonolith::Transform<double, 1, 1>,
        public moonolith::Transform<double, 1, 2> {
    public:
        Transform1(const libMesh::Elem &elem)
        : elem_(elem)
        {}

        void transform_to_reference(const libMesh::Point &world, libMesh::Point &ref) const override;
        void apply(const libMesh::Point &ref, libMesh::Point &world) const override;

        inline bool apply(const moonolith::Vector<double, 1> &in, moonolith::Vector<double, 1> &out) override
        {
            libMesh::Point ref(in.x), world;
            apply(ref, world);
            out.x = world(0);
            //FIXME
            return true;
        }

        inline bool apply_inverse(const moonolith::Vector<double, 1> &in, moonolith::Vector<double, 1> &out) override
        {
            libMesh::Point ref, world(in.x);
            transform_to_reference(ref, world);
            out.x = ref(0);
            //FIXME
            return true;
        }

        inline bool apply(const moonolith::Vector<double, 1> &in, moonolith::Vector<double, 2> &out) override
        {
            libMesh::Point ref(in.x), world;
            apply(ref, world);
            out.x = world(0);
            out.y = world(1);
            //FIXME
            return true;
        }

        inline bool apply_inverse(const moonolith::Vector<double, 2> &in, moonolith::Vector<double, 1> &out) override
        {
            libMesh::Point ref, world(in.x, in.y);
            
            transform_to_reference(ref, world);
            out.x = ref(0);
            //FIXME
            return true;
        }

    private:
        const libMesh::Elem &elem_;
    };


    class Transform2 :
        public Transform,
        public moonolith::Transform<double, 2, 2>,
        public moonolith::Transform<double, 2, 3>  {
    public:
        // Transform2(const libMesh::DenseMatrix<libMesh::Real> &polygon)
        // : polygon_(polygon)
        // {
        //  assert(polygon_.m() == 3 || polygon_.m() == 4 && "must be either a triangle or a quad");
        // }

        Transform2(const libMesh::Elem &elem)
        : elem_(elem)
        {}

        void transform_to_reference(const libMesh::Point &world, libMesh::Point &ref) const override;
        void apply(const libMesh::Point &ref, libMesh::Point &world) const override;

        static void compute_affine_transformation(const libMesh::Elem * elem, libMesh::DenseMatrix<libMesh::Real> &A_inv);

        inline bool apply(const moonolith::Vector<double, 2> &in, moonolith::Vector<double, 2> &out) override
        {
            libMesh::Point ref(in.x, in.y), world;
            apply(ref, world);
            out.x = world(0);
            out.y = world(1);
           
            //FIXME
            return true;
        }

        inline bool apply_inverse(const moonolith::Vector<double, 2> &in, moonolith::Vector<double, 2> &out) override
        {
            libMesh::Point ref, world(in.x, in.y);
            transform_to_reference(ref, world);
            out.x = ref(0);
            out.y = ref(1);
            //FIXME
            return true;
        }

        inline bool apply(const moonolith::Vector<double, 2> &in, moonolith::Vector<double, 3> &out) override
        {
            libMesh::Point ref(in.x, in.y), world;
            apply(ref, world);
            out.x = world(0);
            out.y = world(1);
            out.z = world(2);
            //FIXME
            return true;
        }

        inline bool apply_inverse(const moonolith::Vector<double, 3> &in, moonolith::Vector<double, 2> &out) override
        {
            libMesh::Point ref, world(in.x, in.y, in.z);
            
            transform_to_reference(ref, world);
            out.x = ref(0);
            out.y = ref(1);
            //FIXME
            return true;
        }

    private:

        const libMesh::Elem &elem_;

    };

    class Transform3 :
        public Transform,
        public moonolith::Transform<double, 3, 3> {
    public:

        Transform3(const libMesh::Elem &elem)
        : elem_(elem)
        { }

        void transform_to_reference(const libMesh::Point &world, libMesh::Point &refm) const override;
        void apply(const libMesh::Point &ref, libMesh::Point &world) const override;


        inline bool apply(const moonolith::Vector<double, 3> &in, moonolith::Vector<double, 3> &out) override
        {
            libMesh::Point ref(in.x, in.y, in.z), world;
            apply(ref, world);
            out.x = world(0);
            out.y = world(1);
            out.z = world(2);
           
            //FIXME
            return true;
        }

        inline bool apply_inverse(const moonolith::Vector<double, 3> &in, moonolith::Vector<double, 3> &out) override
        {
            libMesh::Point ref, world(in.x, in.y, in.z);
            
            transform_to_reference(ref, world);
            out.x = ref(0);
            out.y = ref(1);
            out.z = ref(2);
           
            //FIXME
            return true;
        }

    private:
        // const Polyhedron &polyhedron_;
        const libMesh::Elem &elem_;

    };


    class AffineTransform2 :
        public Transform,
        public moonolith::Transform<double, 2, 2> {
    public:
        AffineTransform2(const libMesh::Elem &elem)
        {
            compute_affine_transformation(elem, A_inv_ , A_inv_m_b_);
        }

        AffineTransform2(const libMesh::DenseMatrix<libMesh::Real> &A_inv,
                         const libMesh::DenseVector<libMesh::Real> A_inv_m_b)
        : A_inv_(A_inv), A_inv_m_b_(A_inv_m_b)
        {}

        void transform_to_reference(const libMesh::Point &world, libMesh::Point &ref) const override;
        void apply(const libMesh::Point &ref, libMesh::Point &world) const override { assert(false && "implement me"); }

        AffineTransform2() {}


        libMesh::DenseMatrix<libMesh::Real> &A_inv()
        {
            return A_inv_;
        }

        libMesh::DenseVector<libMesh::Real> &A_inv_m_b()
        {
            return A_inv_m_b_;
        }

        inline bool apply(const moonolith::Vector<double, 2> &in, moonolith::Vector<double, 2> &out) override
        {
            libMesh::Point ref(in.x, in.y), world;
            apply(ref, world);
            out.x = world(0);
            out.y = world(1);
            //FIXME
            return true;
        }

        inline bool apply_inverse(const moonolith::Vector<double, 2> &in, moonolith::Vector<double, 2> &out) override
        {
            libMesh::Point ref, world(in.x, in.y);
            
            transform_to_reference(ref, world);
            out.x = ref(0);
            out.y = ref(1);
            //FIXME
            return true;
        }

    private:
        libMesh::DenseMatrix<libMesh::Real> A_inv_;
        libMesh::DenseVector<libMesh::Real> A_inv_m_b_;



        static void compute_affine_transformation(const libMesh::Elem &elem,
                                                  libMesh::DenseMatrix<libMesh::Real> &A_inv,
                                                  libMesh::DenseVector<libMesh::Real> &A_inv_m_b);
    };


    class AffineTransform3 :
        public Transform,
        public moonolith::Transform<double, 3, 3> {
    public:
        AffineTransform3(const libMesh::Elem &elem)
        {
            compute_affine_transformation(elem, A_inv_ , A_inv_m_b_);
        }

        void transform_to_reference(const libMesh::Point &world, libMesh::Point &ref) const override;
        void apply(const libMesh::Point &ref, libMesh::Point &world) const override { assert(false && "implement me"); }

        AffineTransform3(const libMesh::DenseMatrix<libMesh::Real> &A_inv,
                         const libMesh::DenseVector<libMesh::Real> A_inv_m_b)
        : A_inv_(A_inv), A_inv_m_b_(A_inv_m_b)
        {}

        AffineTransform3()
        {}

        libMesh::DenseMatrix<libMesh::Real> &A_inv()
        {
            return A_inv_;
        }

        libMesh::DenseVector<libMesh::Real> &A_inv_m_b()
        {
            return A_inv_m_b_;
        }


        inline bool apply(const moonolith::Vector<double, 3> &in, moonolith::Vector<double, 3> &out) override
        {
            libMesh::Point ref(in.x, in.y, in.z), world;
            apply(ref, world);
            out.x = world(0);
            out.y = world(1);
            out.z = world(2);
           
            //FIXME
            return true;
        }

        inline bool apply_inverse(const moonolith::Vector<double, 3> &in, moonolith::Vector<double, 3> &out) override
        {
            libMesh::Point ref, world(in.x, in.y, in.z);
            
            transform_to_reference(ref, world);
            out.x = ref(0);
            out.y = ref(1);
            out.z = ref(2);
           
            //FIXME
            return true;
        }

    private:
        libMesh::DenseMatrix<libMesh::Real> A_inv_;
        libMesh::DenseVector<libMesh::Real> A_inv_m_b_;

        static void compute_affine_transformation(const libMesh::Elem &elem,
                                                  libMesh::DenseMatrix<libMesh::Real> &A_inv,
                                                  libMesh::DenseVector<libMesh::Real> &A_inv_m_b);
    };

    class SideAffineTransform3 :
        public Transform,
        public moonolith::Transform<double, 2, 3>,
        public moonolith::Transform<double, 3, 3>  {
    public:
        inline SideAffineTransform3(const libMesh::Elem &elem, const int side)
        : a_trafo_()
        {
            compute_affine_transformation(elem, side, a_trafo_.A_inv(), a_trafo_.A_inv_m_b());
        }

        inline void transform_to_reference(const libMesh::Point &world, libMesh::Point &ref) const override
        {
            a_trafo_.transform_to_reference(world, ref);
            assert( std::abs(ref(2)) < 1e-8 );
        }

        void apply(const libMesh::Point &ref, libMesh::Point &world) const override { assert(false && "implement me"); }


        inline bool apply(const moonolith::Vector<double, 2> &in, moonolith::Vector<double, 3> &out) override
        {
            libMesh::Point ref(in.x, in.y), world;
            apply(ref, world);
            out.x = world(0);
            out.y = world(1);
            out.z = world(2);
            //FIXME
            return true;
        }

        inline bool apply_inverse(const moonolith::Vector<double, 3> &in, moonolith::Vector<double, 2> &out) override
        {
            libMesh::Point ref, world(in.x, in.y, in.z);
            
            transform_to_reference(ref, world);
            out.x = ref(0);
            out.y = ref(1);
            //FIXME
            return true;
        }

        inline bool apply(const moonolith::Vector<double, 3> &in, moonolith::Vector<double, 3> &out) override
        {
            libMesh::Point ref(in.x, in.y, in.z), world;
            apply(ref, world);
            out.x = world(0);
            out.y = world(1);
            out.z = world(2);
           
            //FIXME
            return true;
        }

        inline bool apply_inverse(const moonolith::Vector<double, 3> &in, moonolith::Vector<double, 3> &out) override
        {
            libMesh::Point ref, world(in.x, in.y, in.z);
            
            transform_to_reference(ref, world);
            out.x = ref(0);
            out.y = ref(1);
            out.z = ref(2);
           
            //FIXME
            return true;
        }

    private:
        AffineTransform3 a_trafo_;

        static void compute_affine_transformation(const libMesh::Elem &elem,
                                                  const int side,
                                                  libMesh::DenseMatrix<libMesh::Real> &A_inv,
                                                  libMesh::DenseVector<libMesh::Real> &A_inv_m_b);
    };

    class SideAffineTransform2 :
        public Transform,
        public moonolith::Transform<double, 1, 2>,
        public moonolith::Transform<double, 2, 2> {
    public:
        inline SideAffineTransform2(const libMesh::Elem &elem, const int side)
        : a_trafo_()
        {
            compute_affine_transformation(elem, side, a_trafo_.A_inv(), a_trafo_.A_inv_m_b());
        }

        inline void transform_to_reference(const libMesh::Point &world, libMesh::Point &ref) const override
        {
            a_trafo_.transform_to_reference(world, ref);

            //reference segment is (-1, 1)
            ref(0) *= 2.;
            ref(0) -= 1.;

            assert( std::abs(ref(1)) < 1e-8 );
            assert( std::abs(ref(2)) < 1e-8 );
        }

        void apply(const libMesh::Point &ref, libMesh::Point &world) const override { assert(false && "implement me"); }

        inline bool apply(const moonolith::Vector<double, 1> &in, moonolith::Vector<double, 2> &out) override
        {
            libMesh::Point ref(in.x), world;
            apply(ref, world);
            out.x = world(0);
            out.y = world(1);
            //FIXME
            return true;
        }

        inline bool apply_inverse(const moonolith::Vector<double, 2> &in, moonolith::Vector<double, 1> &out) override
        {
            libMesh::Point ref, world(in.x, in.y);
            
            transform_to_reference(ref, world);
            out.x = ref(0);
            //FIXME
            return true;
        }

        inline bool apply(const moonolith::Vector<double, 2> &in, moonolith::Vector<double, 2> &out) override
        {
            libMesh::Point ref(in.x, in.y), world;
            apply(ref, world);
            out.x = world(0);
            out.y = world(1);
           
            //FIXME
            return true;
        }

        inline bool apply_inverse(const moonolith::Vector<double, 2> &in, moonolith::Vector<double, 2> &out) override
        {
            libMesh::Point ref, world(in.x, in.y);
            transform_to_reference(ref, world);
            out.x = ref(0);
            out.y = ref(1);
            //FIXME
            return true;
        }

    private:
        AffineTransform2 a_trafo_;

        static void compute_affine_transformation(const libMesh::Elem &elem,
                                                  const int side,
                                                  libMesh::DenseMatrix<libMesh::Real> &A_inv,
                                                  libMesh::DenseVector<libMesh::Real> &A_inv_m_b);
    };
}


#endif //UTOPIA_LIBMESH_TRANSFORM_HPP
