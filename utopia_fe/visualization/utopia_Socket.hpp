#ifndef UTOPIA_SOCKET_HPP
#define UTOPIA_SOCKET_HPP

#include "utopia.hpp"
#include "utopia_fe_base.hpp"

#include <string>
#include <vector>

//forward declarations
namespace libMesh {
	class MeshBase;
	class Point;
}

namespace utopia {
	class Box;

	void __attribute__ ((used)) plot_polygon(const int dims, const int n_points, const double *points, const std::string &name = "polygon/p");
	void __attribute__ ((used)) plot_lines(const int dims, const int n_points, const double *points, const std::string &name = "lines/l");
	void __attribute__ ((used)) plot_points(const int dims, const int n_points, const double *points, const std::string &name = "points/p");
	void __attribute__ ((used)) quiver(const int dims, const int n_points, const double *points, const double *vectors, const std::string &name = "quiver/q");
	void __attribute__ ((used)) plot_quad_points(const int dims, const std::vector<libMesh::Point> &points, const std::string &name = "quad_points/qp");
	void __attribute__ ((used)) plot_mesh(const libMesh::MeshBase &mesh, const std::string &name = "mesh/m");
	void __attribute__ ((used)) plot_mesh_f(const libMesh::MeshBase &mesh, const double * function, const std::string &name = "mesh_f/m");
	void __attribute__ ((used)) plot_polygon_mesh(const libMesh::MeshBase &mesh, const std::string &name = "p_mesh/m");
	void __attribute__ ((used)) plot_box(const Box &box, const std::string &name = "box/b");
	// void plot_mesh_surf_elements(const libMesh::MeshBase &mesh, const std::string &name);

	void plot_scaled_normal_field(
		const libMesh::MeshBase &mesh,
		const UVector &normals,
		const UVector &scale,
		const std::string &name = "normal_field/nf");
}

#endif //UTOPIA_SOCKET_HPP
