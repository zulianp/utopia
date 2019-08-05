#ifndef UTOPIA_CONTACT_SIM_PARAMS_HPP
#define UTOPIA_CONTACT_SIM_PARAMS_HPP

#include <string>

#include "libmesh/point.h"
#include "libmesh/dense_vector.h"

namespace utopia {

    typedef struct {
        std::string mesh_path;
        int boundary_tag_1;
        int boundary_tag_2;
        int boundary_tag_3;

        double search_radius;

        double dirichlet_value_1;
        double dirichlet_value_2;
        double dirichlet_value_3;

        bool strict_search_radius;
    } ContactSimParams;

    static const double LARGE_VAL = 1e8;

    static const ContactSimParams contact18	  			  = { "../data/contact18.e", 1, 2, -1, 0.1, -0.2, 0.2, LARGE_VAL, false };
    static const ContactSimParams contact8	  			  = { "../data/contact8.e", 1, 2, -1, 0.2, 0, 0, LARGE_VAL,false };
    static const ContactSimParams contact8_tris	  		  = { "../data/contact8_tris.e", 1, 2, -1, 0.1, -0.2, 0.2, LARGE_VAL,false };
    static const ContactSimParams contact_3D_tets	      = { "../data/contact_3D_tets.e", 1, 2, -1, 0.1, -0.2, 0.2, LARGE_VAL,false };
    static const ContactSimParams contact_cylinder	  	  = { "../data/contact_cylinder.e", 1, 2, -1, 0.5, -0.1, 0.1, LARGE_VAL,false };
    static const ContactSimParams contact_sphere	  	  = { "../data/contact_sphere.e", 1, 2, -1, 0.3, -0.2, 0.2, LARGE_VAL,false };
    static const ContactSimParams contact_cuboids	  	  = { "../data/contact_cuboids.e", 1, 2, -1, 0.1, -0.2, 0.2, LARGE_VAL,false };
    static const ContactSimParams triple_contact_circle	  = { "../data/triple_contact_circle.e", 1, 2, -1, 0.2, -0.2, 0.2, LARGE_VAL,false };
    static const ContactSimParams contact_circles	  	  = { "../data/contact_circles_no_tag.e", 1, 2, -1, 0.8, -0.2, 0.2, LARGE_VAL,false };
    static const ContactSimParams contact2	  			  = { "../data/contact3D.e", 1, 2, -1, 0.1, -0.2, 0.2, LARGE_VAL,false };
    static const ContactSimParams contact_quads	  		  = { "../data/contact_quads.e", 1, 2, -1, 0.1, -0.2, 0.2, LARGE_VAL,false };
    static const ContactSimParams multi_contact	  		  = { "../data/multi_contact.e", 1, 2, -1, 0.05, 0.0, 0.0, LARGE_VAL,false };
    static const ContactSimParams multi_contact_3D_3	  = { " /Users/patrick/Desktop/Code/moonolith_meshes/exodus/contact_ti_coarse.e", 1, 0, -1, 0.1, -0.2, 0.2, LARGE_VAL,false };
    static const ContactSimParams leaflets_contact	  	  = { "../data/leaflets.e", 101, 102, 103, 0.1, 0.05, 0.05, 0.05, true };
    static const ContactSimParams multi_contact_3D	  	  = { "/Users/patrick/Desktop/Code/moonolith_meshes/exodus/contact_mb_4.e", 1, 2, -1, 0.28, -0.59, 1.1, LARGE_VAL, false };
    static const ContactSimParams multi_contact_3D_2	  = { "/Users/patrick/Desktop/Code/moonolith/cutlibpp/data/multiple_contact/multiple_contact_fine.e", 1, 2, -1, 0.3, -2., 0., LARGE_VAL, false };
    static const ContactSimParams contact_least_squares   = { "../data/contact_least_squares.e", 1, 2, -1, 0.1, -0.2, 0.2, LARGE_VAL, false } ;
    static const ContactSimParams contact_least_squares_2 = { "../data/contact_least_squares_2.e", 1, 2, -1, 0.1, -0.2, 0.2, LARGE_VAL, false } ;
    static const ContactSimParams multi_contact_quads 	  = { "../data/multi_contact_quads.e", 1, 2, -1, 0.01, -0.03, 0.03, LARGE_VAL, false };
    static const ContactSimParams hip_femure_contact	  = { "/Users/patrick/Downloads/ASCII_bone/all_sidesets.e", 1, 2, -1, 10.0, 10.0, 0.0, LARGE_VAL, false };
    static const ContactSimParams implant_contact	      = { "../data/implant.e", 1, 11, -1, 4.0, -4.0, 0.0, LARGE_VAL, false };
    static const ContactSimParams contact_cubes	      	  = { "../data/multibody.e", 1, 2, -1, 0.1, -0.2, 0.0, LARGE_VAL, false };
    static const ContactSimParams hertz_contact			  = { "../data/hertz_contact.e", 1, 2, -1, 0.4, -0.07, 0.07, LARGE_VAL, false };
    static const ContactSimParams hertz_contact_coarse	  = { "../data/hertz_contact_coarse.e", 1, 2, -1, 0.4, -0.1, 0.1, LARGE_VAL, false };

    inline static void upper_boundary_cond(const libMesh::Point & p, libMesh::DenseVector<libMesh::Real> & output)
    {
        if(output.size() < 3) {
            output.resize(3);
        }

        output.zero();
        output(1) = -0.2;
    }

    inline static void lower_boundary_cond(const libMesh::Point & p, libMesh::DenseVector<libMesh::Real> & output)
    {
        if(output.size() < 3) {
            output.resize(3);
        }

        output.zero();
        output(1) = 0.2;
    }

}

#endif //UTOPIA_CONTACT_SIM_PARAMS_HPP
