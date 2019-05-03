#ifndef UTOPIA_LIBMESH_TOMOONOLITH_CONVERTIONS_HPP
#define UTOPIA_LIBMESH_TOMOONOLITH_CONVERTIONS_HPP

#include "MortarAssemble.hpp"

#include "moonolith_polygon.hpp"
#include "moonolith_line.hpp"
#include "moonolith_polyhedron.hpp"
#include "moonolith_polygon.hpp"
#include "moonolith_line.hpp"
#include "moonolith_map_quadrature.hpp"
#include "moonolith_affine_transform.hpp"
#include "moonolith_elem_shape.hpp"
#include "moonolith_vector.hpp"
#include "moonolith_elem_triangle.hpp"
#include "moonolith_elem_quad.hpp"
#include "moonolith_elem_segment.hpp"

#include "libmesh/point.h"
#include "libmesh/elem.h"

namespace utopia {

        template<typename T, int Dim>
        inline void make(const libMesh::Point &p, moonolith::Vector<T, Dim> &q)
        {
            for(int i = 0; i < Dim; ++i) {
                q[i] = p(i);
            }
        }

        template<typename T, int Dim>
        inline void make(const libMesh::Elem &elem, moonolith::Line<T, Dim> &poly)
        {   
            make(elem.node_ref(0), poly.p0);
            make(elem.node_ref(1), poly.p1);
        }


        template<typename T, int Dim>
        inline void make_non_affine(const libMesh::Elem &elem, moonolith::Storage<moonolith::Vector<T, Dim>> &poly_line)
        {   
            const int n_nodes = elem.n_nodes();

            poly_line.resize(n_nodes);

            if(n_nodes == 2) {
                //P1
                make(elem.node_ref(0), poly_line[0]);
                make(elem.node_ref(1), poly_line[1]);
            } else if(n_nodes == 3) {
                //P2
                make(elem.node_ref(0), poly_line[0]);
                make(elem.node_ref(2), poly_line[1]);
                make(elem.node_ref(1), poly_line[2]);
            } else if(n_nodes == 4) {
                //P3
                make(elem.node_ref(0), poly_line[0]);
                make(elem.node_ref(2), poly_line[1]);
                make(elem.node_ref(3), poly_line[2]);
                make(elem.node_ref(1), poly_line[3]);
            } else {
                assert(false);
            }

        }

        template<typename T, int Dim>
        inline void make_non_affine(const libMesh::Elem &elem, moonolith::PolyLine<T, Dim> &poly_line)
        {
            make_non_affine(elem, poly_line.points);
        }

        template<typename T, int Dim>
        inline void make(const libMesh::Elem &elem, moonolith::Polygon<T, Dim> &poly)
        {
           assert(is_tri(elem.type()) || is_quad(elem.type()));
           auto n_nodes = is_tri(elem.type()) ? 3 : 4;

           poly.resize(n_nodes);

           for(auto i = 0; i < n_nodes; ++i) {
               const auto &p = elem.node_ref(i);
               auto &q = poly.points[i];
               make(p, q);
           }
        }


        template<typename T, int Dim>
        inline void make_non_affine(const libMesh::Elem &elem, moonolith::Polygon<T, Dim> &poly)
        {
            assert(is_tri(elem.type()) || is_quad(elem.type()));
            auto n_corner_nodes = is_tri(elem.type()) ? 3 : 4;
            auto n_nodes = elem.n_nodes();

            if(n_nodes == n_corner_nodes) {
                //fallback to affine case
                return make(elem, poly);
            }

            //P2 only
            assert(n_nodes = n_corner_nodes * 2);

            poly.resize(n_nodes);

            for(auto i = 0; i < n_corner_nodes; ++i) {
                const auto i2 = i*2;
                const auto &p = elem.node_ref(i);
                const auto &q = elem.node_ref(i + n_corner_nodes);
                make(p, poly.points[i2]);
                make(q, poly.points[i2 + 1]);
            }
        }

        template<typename T>
        inline void make_tetrahedron(const libMesh::Elem &elem, moonolith::Polyhedron<T> &poly)
        {
            static const std::size_t n_nodes = 4;

            assert(elem.dim() == 3);
            assert(is_tet(elem.type()));

            poly.el_ptr.resize(n_nodes + 1);
            poly.el_index.resize(12);
            poly.points.resize(n_nodes);

            poly.el_ptr[0] = 0;
            poly.el_ptr[1] = 3;
            poly.el_ptr[2] = 6;
            poly.el_ptr[3] = 9;
            poly.el_ptr[4] = 12;

            for(std::size_t i = 0; i < n_nodes; ++i) {
                const  auto &p = elem.node_ref(i);
                poly.points[i].x = p(0);
                poly.points[i].y = p(1);
                poly.points[i].z = p(2);
            }

            poly.el_index[0] = 0;
            poly.el_index[1] = 2;
            poly.el_index[2] = 1;
       
            poly.el_index[3] = 0;
            poly.el_index[4] = 3;
            poly.el_index[5] = 2;
       
            poly.el_index[6] = 0;
            poly.el_index[7] = 1;
            poly.el_index[8] = 3;
       
            poly.el_index[9] = 1;
            poly.el_index[10] = 2;
            poly.el_index[11] = 3;   

            poly.type = moonolith::Polyhedron<T>::TET;     
        }

        template<typename T>
        void make_pyramid(const libMesh::Elem &elem, moonolith::Polyhedron<T> &poly)
        {
            static const std::size_t n_nodes = 5;
            assert(elem.dim() == 3);
            assert(is_pyramid(elem.type()));

            poly.el_ptr.resize(5+1);
            poly.el_index.resize(16);
            poly.points.resize(n_nodes);

            for(std::size_t i = 0; i < n_nodes; ++i) {
                const  auto &p = elem.node_ref(i);
                poly.points[i].x = p(0);
                poly.points[i].y = p(1);
                poly.points[i].z = p(2);
            }

            poly.el_ptr[0] = 0;
            poly.el_ptr[1] = 3;
            poly.el_ptr[2] = 6;
            poly.el_ptr[3] = 9;
            poly.el_ptr[4] = 12;
            poly.el_ptr[5] = 16;

            //face 0
            poly.el_index[0] = 0;
            poly.el_index[1] = 1;
            poly.el_index[2] = 4;

            //face 1
            poly.el_index[3] = 1;
            poly.el_index[4] = 2;
            poly.el_index[5] = 4;

            //face 2
            poly.el_index[6] = 3;
            poly.el_index[7] = 0;
            poly.el_index[8] = 4;

            //face 3
            poly.el_index[9]  = 2;
            poly.el_index[10] = 3;
            poly.el_index[11] = 4;

            //face 4
            poly.el_index[12] = 0;
            poly.el_index[13] = 3;
            poly.el_index[14] = 2;
            poly.el_index[15] = 1;

            poly.type = moonolith::Polyhedron<T>::UNSTRUCTURED;
        }

        template<typename T>
        void make_prism(const libMesh::Elem &elem, moonolith::Polyhedron<T> &poly)
        {
            assert(elem.dim() == 3);
            assert(is_prism(elem.type()));

            static const std::size_t n_nodes = 6;

            poly.el_ptr.resize(5+1);
            poly.el_index.resize(18);
            poly.points.resize(n_nodes);

            for(std::size_t i = 0; i < n_nodes; ++i) {
                const  auto &p = elem.node_ref(i);
                poly.points[i].x = p(0);
                poly.points[i].y = p(1);
                poly.points[i].z = p(2);
            }

            poly.el_ptr[0] = 0;
            poly.el_ptr[1] = 3;
            poly.el_ptr[2] = 6;
            poly.el_ptr[3] = 10;
            poly.el_ptr[4] = 14;
            poly.el_ptr[5] = 18;

            //face 0
            poly.el_index[0] = 0;
            poly.el_index[1] = 1;
            poly.el_index[2] = 2;

            //face 1
            poly.el_index[3] = 4;
            poly.el_index[4] = 3;
            poly.el_index[5] = 5;

            //face 2
            poly.el_index[6] = 0;
            poly.el_index[7] = 2;
            poly.el_index[8] = 5;
            poly.el_index[9] = 3;

            //face 3
            poly.el_index[10] = 1;
            poly.el_index[11] = 4;
            poly.el_index[12] = 5;
            poly.el_index[13] = 2;

            //face 4
            poly.el_index[14] = 0;
            poly.el_index[15] = 1;
            poly.el_index[16] = 4;
            poly.el_index[17] = 3;

            poly.type = moonolith::Polyhedron<T>::UNSTRUCTURED;
        }

        template<typename T>
        void make_hex(const libMesh::Elem &elem, moonolith::Polyhedron<T> &poly)
        {
            assert(elem.dim() == 3);
            assert(is_hex(elem.type()));

            static const std::size_t n_nodes = 8;
            poly.el_ptr.resize(6+1);
            poly.el_index.resize(6*4);
            poly.points.resize(n_nodes);

            for(std::size_t i = 0; i < n_nodes; ++i) {
                const  auto &p = elem.node_ref(i);
                poly.points[i].x = p(0);
                poly.points[i].y = p(1);
                poly.points[i].z = p(2);
            }

            poly.el_ptr[0] = 0;
            poly.el_ptr[1] = 4;
            poly.el_ptr[2] = 8;
            poly.el_ptr[3] = 12;
            poly.el_ptr[4] = 16;
            poly.el_ptr[5] = 20;
            poly.el_ptr[6] = 24;

            //face 0
            poly.el_index[0] = 0;
            poly.el_index[1] = 1;
            poly.el_index[2] = 5;
            poly.el_index[3] = 4;

            //face 1
            poly.el_index[4] = 1;
            poly.el_index[5] = 2;
            poly.el_index[6] = 6;
            poly.el_index[7] = 5;

            //face 2
            poly.el_index[8]  = 3;
            poly.el_index[9]  = 7;
            poly.el_index[10] = 6;
            poly.el_index[11] = 2;

            //face 3
            poly.el_index[12] = 0;
            poly.el_index[13] = 4;
            poly.el_index[14] = 7;
            poly.el_index[15] = 3;

            //face 4
            poly.el_index[16] = 2;
            poly.el_index[17] = 1;
            poly.el_index[18] = 0;
            poly.el_index[19] = 3;

            //face 5
            poly.el_index[20] = 6;
            poly.el_index[21] = 7;
            poly.el_index[22] = 4;
            poly.el_index[23] = 5;

            poly.type = moonolith::Polyhedron<T>::HEX;
        }

        template<typename T>
        inline void make(const libMesh::Elem &elem, moonolith::Polyhedron<T> &poly)
        {
            if(is_tet(elem.type())) {
                make_tetrahedron(elem, poly);
                return;
            }

            if(is_hex(elem.type())) {
                make_hex(elem, poly);
                return;
            }

            if(is_pyramid(elem.type())) {
                make_pyramid(elem, poly);
                return;
            }

            if(is_prism(elem.type())) {
                make_prism(elem, poly);
                return;
            }

            assert(false);
        }


    template<typename T, int Dim>
    inline void normal(const libMesh::Elem &elem, moonolith::Vector<T, Dim> &nn)
    {
        using namespace libMesh;
        Point o, u, v, n;

        if(Dim == 2) {
            assert(elem.n_nodes() >= 2);
            o = elem.point(0);
            u = elem.point(1);
            u -= o;
            n(0) =  u(1);
            n(1) = -u(0);

        } else {
            assert(Dim == 3);
            o = elem.point(0);
            u = elem.point(1);
            v = elem.point(2);
            u -= o;
            v -= o;
            n = u.cross(v);
        }

        n *= 1./n.norm();

        make(n, nn);
    }

    template<typename T, int Dim>
    inline void convert(
        const libMesh::QGauss &q_in,
        const moonolith::Vector<T, Dim> &point_shift,
        const T &point_rescale,
        const T &weight_rescale,
        moonolith::Quadrature<T, Dim> &q_out)
    {
        const auto &p_in = q_in.get_points();
        const auto &w_in = q_in.get_weights();

        const std::size_t n_qp = q_in.n_points();
        
        q_out.resize(n_qp);

        for(std::size_t k = 0; k < n_qp; ++k) {
            q_out.weights[k] = w_in[k] * weight_rescale;

            const auto &pk_in = p_in[k];
            auto &pk_out = q_out.points[k];

            for(int i = 0; i < Dim; ++i) {
                pk_out[i] = (pk_in(i) + point_shift[i]) * point_rescale;
            }
        }
    }

    template<typename T, int Dim>
    inline void convert(
        const moonolith::Quadrature<T, Dim> &q_in,
        const moonolith::Vector<T, Dim> &point_shift,
        const T &point_rescale,
        const T &weight_rescale,
        QMortar &q_out)
    {
        const auto &p_in = q_in.points;
        const auto &w_in = q_in.weights;

        auto &p_out = q_out.get_points();
        auto &w_out = q_out.get_weights();

        const std::size_t n_qp = q_in.n_points();
        
        q_out.resize(n_qp);

        for(std::size_t k = 0; k < n_qp; ++k) {
            w_out[k] = w_in[k] * weight_rescale;

            const auto &pk_in = p_in[k];
            auto &pk_out = p_out[k];

            for(int i = 0; i < Dim; ++i) {
                pk_out(i) = pk_in[i] * point_rescale + point_shift[i];
            }
        }
    }

    template<typename T, int Dim>
    inline void convert(
        const moonolith::Quadrature<T, Dim> &q_in,
        const T &weight_rescale,
        QMortar &q_out)
    {
        static const moonolith::Vector<T, Dim> zero;
        return convert(q_in, zero, 1., weight_rescale, q_out);
    }

    // inline void make_line_transform(
    //     const libMesh::Elem &elem,
    //     moonolith::AffineTransform<double, 1, 2> &trafo)
    // {
    //     const auto &q0 = elem.node_ref(0);
    //     const auto &q1 = elem.node_ref(1);

    //     moonolith::Vector<double, 2> p0, p1;

    //     p0.x = q1(0);
    //     p0.y = q0(1);

    //     p1.x = q1(0);
    //     p1.y = q1(1);

    //     make(p0, p1, trafo);
    // }

    // inline void make_triangle_transform(
    //     const libMesh::Elem &elem,
    //     moonolith::AffineTransform<double, 2, 3> &trafo)
    // {
    //     const auto &q0 = elem.node_ref(0);
    //     const auto &q1 = elem.node_ref(1);
    //     const auto &q2 = elem.node_ref(2);

    //     moonolith::Vector<double, 3> p0, p1, p2;

    //     p0.x = q0(0);
    //     p0.y = q0(1);
    //     p0.z = q0(2);

    //     p1.x = q1(0);
    //     p1.y = q1(1);
    //     p1.z = q1(2);

    //     p2.x = q2(0);
    //     p2.y = q2(1);
    //     p2.z = q2(2);

    //     make(p0, p1, p2, trafo);
    // }

    inline void make_transform(
        const libMesh::Elem &elem,
        std::shared_ptr<moonolith::Transform<double, 2, 3>> &trafo)
    {
        trafo = std::make_shared<Transform2>(elem);
    }   

    inline void make_transform(
        const libMesh::Elem &elem,
        std::shared_ptr<moonolith::Transform<double, 1, 2>> &trafo)
    {
        trafo = std::make_shared<Transform1>(elem);
    }   

    inline void make_transform(const libMesh::Elem &elem, moonolith::AffineTransform<double, 2, 2> &trafo)
    {
        libMesh::Point p0(0.0, 0.0, 0.0);
        libMesh::Point p1(1.0, 0.0, 0.0);
        libMesh::Point p2(0.0, 1.0, 0.0);

        // figure out the ordering of the map
        p0 = libMesh::FE<2, libMesh::LAGRANGE>::map(&elem, p0);
        p1 = libMesh::FE<2, libMesh::LAGRANGE>::map(&elem, p1);
        p2 = libMesh::FE<2, libMesh::LAGRANGE>::map(&elem, p2);

        moonolith::Vector<double, 2> q0, q1, q2;

        q0.x = p0(0);
        q0.y = p0(1);

        q1.x = p1(0);
        q1.y = p1(1);

        q2.x = p2(0);
        q2.y = p2(1);

        moonolith::make(q0, q1, q2, trafo);
    }

    inline void make_transform(const libMesh::Elem &elem, moonolith::AffineTransform<double, 3, 3> &trafo)
    {
        libMesh::Point p0(0.0, 0.0, 0.0);
        libMesh::Point p1(1.0, 0.0, 0.0);
        libMesh::Point p2(0.0, 1.0, 0.0);
        libMesh::Point p3(0.0, 0.0, 1.0);

        // figure out the ordering of the map
        p0 = libMesh::FE<3, libMesh::LAGRANGE>::map(&elem, p0);
        p1 = libMesh::FE<3, libMesh::LAGRANGE>::map(&elem, p1);
        p2 = libMesh::FE<3, libMesh::LAGRANGE>::map(&elem, p2);
        p3 = libMesh::FE<3, libMesh::LAGRANGE>::map(&elem, p3);

        moonolith::Vector<double, 3> q0, q1, q2, q3;

        q0.x = p0(0);
        q0.y = p0(1);
        q0.z = p0(2);

        q1.x = p1(0);
        q1.y = p1(1);
        q1.z = p1(2);

        q2.x = p2(0);
        q2.y = p2(1);
        q2.z = p2(2);

        q3.x = p3(0);
        q3.y = p3(1);
        q3.z = p3(2);

        moonolith::make(q0, q1, q2, q3, trafo);
    }

    inline void make_transform(const libMesh::Elem &elem, moonolith::AffineTransform<double, 2, 3> &trafo)
    {
        libMesh::Point p0(0.0, 0.0);
        libMesh::Point p1(1.0, 0.0);
        libMesh::Point p2(0.0, 1.0);

        // figure out the ordering of the map
        p0 = libMesh::FE<2, libMesh::LAGRANGE>::map(&elem, p0);
        p1 = libMesh::FE<2, libMesh::LAGRANGE>::map(&elem, p1);
        p2 = libMesh::FE<2, libMesh::LAGRANGE>::map(&elem, p2);

        moonolith::Vector<double, 3> q0, q1, q2;

        q0.x = p0(0);
        q0.y = p0(1);
        q0.z = p0(2);

        q1.x = p1(0);
        q1.y = p1(1);
        q1.z = p1(2);

        q2.x = p2(0);
        q2.y = p2(1);
        q2.z = p2(2);

        moonolith::make(q0, q1, q2, trafo);
    }

    inline void make_transform(const libMesh::Elem &elem, moonolith::AffineTransform<double, 1, 2> &trafo)
    {
        libMesh::Point p0(-1.0, 0.0);
        libMesh::Point p1(1.0, 0.0);

        // figure out the ordering of the map
        p0 = libMesh::FE<1, libMesh::LAGRANGE>::map(&elem, p0);
        p1 = libMesh::FE<1, libMesh::LAGRANGE>::map(&elem, p1);

        moonolith::Vector<double, 2> q0, q1;

        q0.x = p0(0);
        q0.y = p0(1);

        q1.x = p1(0);
        q1.y = p1(1);

        moonolith::make(q0, q1, trafo);
    }


    template<class E>
    void make_triangle_1(const libMesh::Elem &in, E &out)
    {
        //reverse engineer the ordering from ref points to physical points
        libMesh::Point p(0.0, 0.0, 0.0);
        libMesh::Point ref_p(0.0, 0.0, 0.0);

        ///////////////////////////
        p = libMesh::FE<2, libMesh::LAGRANGE>::map(&in, ref_p);
        make(p, out.node(0));
        ///////////////////////////

        ref_p(0) = 1.0;
        p = libMesh::FE<2, libMesh::LAGRANGE>::map(&in, ref_p);
        make(p, out.node(1));

        ///////////////////////////
        ref_p(0) = 0.0;
        ref_p(1) = 1.0;
        p = libMesh::FE<2, libMesh::LAGRANGE>::map(&in, ref_p);
        make(p, out.node(2));
    }

    template<class E>
    void make_triangle_2(const libMesh::Elem &in, E &out)
    {
        make_triangle_1(in, out);

        //reverse engineer the ordering from ref points to physical points
        libMesh::Point p(0.0, 0.0, 0.0);
        libMesh::Point ref_p(0.5, 0.0, 0.0);

        ///////////////////////////
        p = libMesh::FE<2, libMesh::LAGRANGE>::map(&in, ref_p);
        make(p, out.node(3));

        ///////////////////////////
        // ref_p(0) = 0.0;
        ref_p(1) = 0.5;
        p = libMesh::FE<2, libMesh::LAGRANGE>::map(&in, ref_p);
        make(p, out.node(4));

        ///////////////////////////
        ref_p(0) = 0.0;
        ref_p(1) = 0.5;
        p = libMesh::FE<2, libMesh::LAGRANGE>::map(&in, ref_p);
        make(p, out.node(5));
    }

    template<int Dim>
    void make_element(const libMesh::Elem &in, moonolith::Triangle<double, 1, Dim> &out)
    {
        make_triangle_1(in , out);
    }

    template<int Dim>
    void make_element(const libMesh::Elem &in, moonolith::Triangle<double, 2, Dim> &out)
    {
        make_triangle_2(in, out);
        out.set_affine(in.has_affine_map());
    }   


    template<class E>
    void make_quad_1(const libMesh::Elem &in, E &out)
    {
        //reverse engineer the ordering from ref points to physical points
        libMesh::Point p(0.0, 0.0, 0.0);
        libMesh::Point ref_p(-1.0, -1.0, 0.0);

        ///////////////////////////
        p = libMesh::FE<2, libMesh::LAGRANGE>::map(&in, ref_p);
        make(p, out.node(0));
        ///////////////////////////

        ref_p(0) = 1.0;
        p = libMesh::FE<2, libMesh::LAGRANGE>::map(&in, ref_p);
        make(p, out.node(1));

        ///////////////////////////
        ref_p(1) = 1.0;
        p = libMesh::FE<2, libMesh::LAGRANGE>::map(&in, ref_p);
        make(p, out.node(2));

        ///////////////////////////
        ref_p(0) = -1.0;
        p = libMesh::FE<2, libMesh::LAGRANGE>::map(&in, ref_p);
        make(p, out.node(3));
    }

    template<class E>
    void make_quad_2(const libMesh::Elem &in, E &out)
    {
        make_quad_1(in, out);
        //reverse engineer the ordering from ref points to physical points
        libMesh::Point p(0.0, 0.0, 0.0);
        libMesh::Point ref_p(0.0, -1.0, 0.0);

        ///////////////////////////
        p = libMesh::FE<2, libMesh::LAGRANGE>::map(&in, ref_p);
        make(p, out.node(3));
        ///////////////////////////

        ref_p(0) = 1.0;
        ref_p(1) = 0.0;
        p = libMesh::FE<2, libMesh::LAGRANGE>::map(&in, ref_p);
        make(p, out.node(4));

        ///////////////////////////
        ref_p(0) = 0.0;
        ref_p(1) = 1.0;
        p = libMesh::FE<2, libMesh::LAGRANGE>::map(&in, ref_p);
        make(p, out.node(5));

        ///////////////////////////
        ref_p(0) = -1.0;
        ref_p(1) = 0.0;
        p = libMesh::FE<2, libMesh::LAGRANGE>::map(&in, ref_p);
        make(p, out.node(6));
    }

    template<int Dim>
    void make_element(const libMesh::Elem &in, moonolith::Quad<double, 1, Dim> &out)
    {
        make_quad_1(in, out);
    }

    template<int Dim>
    void make_element(const libMesh::Elem &in, moonolith::Quad<double, 2, Dim> &out)
    {
        make_quad_2(in, out);
        out.set_affine(in.has_affine_map());
    }

    template<class E>
    void make_seg_1(const libMesh::Elem &in, E &out)
    {
        //reverse engineer the ordering from ref points to physical points
        libMesh::Point p(0.0, 0.0, 0.0);
        libMesh::Point ref_p(-1.0, 0.0, 0.0);

        ///////////////////////////
        p = libMesh::FE<1, libMesh::LAGRANGE>::map(&in, ref_p);
        make(p, out.node(0));
        ///////////////////////////

        ref_p(0) = 1.0;
        p = libMesh::FE<1, libMesh::LAGRANGE>::map(&in, ref_p);
        make(p, out.node(1));
        ///////////////////////////
    }

    template<class E>
    void make_seg_2(const libMesh::Elem &in, E &out)
    {
        make_seg_1(in, out);
        //reverse engineer the ordering from ref points to physical points
        libMesh::Point p(0.0, 0.0, 0.0);
        libMesh::Point ref_p(0.0, 0.0, 0.0);

        ///////////////////////////
        p = libMesh::FE<1, libMesh::LAGRANGE>::map(&in, ref_p);
        make(p, out.node(2));
        ///////////////////////////
    }

    template<int Dim>
    void make_element(const libMesh::Elem &in, moonolith::Segment<double, 1, Dim> &out)
    {
        make_seg_1(in, out);
    }

    template<int Dim>
    void make_element(const libMesh::Elem &in, moonolith::Segment<double, 2, Dim> &out)
    {
        make_seg_2(in, out);
        out.set_affine(in.has_affine_map());
    }

    template<class E> using LMS = moonolith::OwnedElemShape<E>;

    template<int Dim>
    std::unique_ptr< moonolith::Shape<double, Dim-1, Dim> > make_shape(const libMesh::Elem &elem, const libMesh::FEType &type);

    template<>
    inline std::unique_ptr< moonolith::Shape<double, 1, 2> > make_shape<2>(const libMesh::Elem &elem, const libMesh::FEType &type)
    {
        using moonolith::Segment;

        if(is_edge(elem.type())) {

             if(type.order == 2) {
                 auto s = moonolith::make_unique< LMS< Segment<double, 2, 2> > >();
                 make_element(elem, s->elem());
                 s->init();
                 return s;
             }

             if(type.order == 1) {
                 auto s = moonolith::make_unique< LMS< Segment<double, 1, 2> > >();
                 make_element(elem, s->elem());
                 s->init();
                 return s;
             }
        }

        assert(false);
        std::cerr << " unsupported type" << std::endl;
        return nullptr;
    }

    template<>
    inline std::unique_ptr< moonolith::Shape<double, 2, 3> > make_shape<3>(const libMesh::Elem &elem, const libMesh::FEType &type)
    {
       using moonolith::Triangle;
       using moonolith::Quad;
      
       if(is_tri(elem.type())) {

            if(type.order == 2) {
                auto s = moonolith::make_unique< LMS< Triangle<double, 2, 3> > >();
                make_element(elem, s->elem());
                s->init();
                return s;
            }

            if(type.order == 1) {
                auto s = moonolith::make_unique< LMS< Triangle<double, 1, 3> > >();
                make_element(elem, s->elem());
                s->init();
                return s;
            }

       } else if(is_quad(elem.type())) {

            if(type.order == 2) {
                auto s = moonolith::make_unique< LMS< Quad<double, 2, 3> > >();
                make_element(elem, s->elem());
                s->init();
                return s;
            }

            if(type.order == 1) {
                auto s = moonolith::make_unique< LMS< Quad<double, 1, 3> > >();
                make_element(elem, s->elem());
                s->init();
                return s;
            }

       }

       assert(false);
       std::cerr << " unsupported type" << std::endl;
       return nullptr;
    }

}

#endif //UTOPIA_LIBMESH_TOMOONOLITH_CONVERTIONS_HPP
