#ifndef UTOPIA_SOCKET_HPP
#define UTOPIA_SOCKET_HPP

#include <string>
#include <vector>
#include <libmesh/point.h>
#include <libmesh/mesh.h>

namespace utopia {
	class Socket {
	private:
		const std::string _host, _service;

		
	public:
		typedef char byte;
		
		Socket(const std::string &_host, const std::string &_service)
			: _host(_host), _service(_service)
		{}

		bool write(const byte * buff, const int size);
	};

	class Box;

	void plot_polygon(const int dims, const int n_points, const double *points, const std::string &name = "p");
	void plot_lines(const int dims, const int n_points, const double *points, const std::string &name = "l");
	void plot_points(const int dims, const int n_points, const double *points, const std::string &name = "p");
	void quiver(const int dims, const int n_points, const double *points, const double *vectors, const std::string &name = "q");
	void plot_quad_points(const int dims, const std::vector<libMesh::Point> &points, const std::string &name = "qp");
	void plot_mesh(const libMesh::MeshBase &mesh, const std::string &name);
	void plot_mesh_f(const libMesh::MeshBase &mesh, const double * function, const std::string &name);
	void plot_polygon_mesh(const libMesh::MeshBase &mesh, const std::string &name);
	void plot_box(const Box &box, const std::string &name = "box");
	// void plot_mesh_surf_elements(const libMesh::MeshBase &mesh, const std::string &name);
}

#endif //UTOPIA_SOCKET_HPP
