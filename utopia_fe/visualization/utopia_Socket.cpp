
#include "utopia_Socket.hpp"

#include "libmesh/elem.h"
#include <libmesh/point.h>
#include <libmesh/mesh.h>

#include <sstream>
#include <string>
#include <ostream>
#include <iostream>

#include "MortarAssemble.hpp"
#include "Box.hpp"

#include "utopia_fe_config.hpp"

#ifdef WITH_BOOST

#include "utopia_NormalTangentialCoordinateSystem.hpp"


#include <boost/asio.hpp>
#include <boost/array.hpp>

static const char * host = "localhost";
static const char * port = "24200";

using boost::asio::io_service;
using boost::asio::ip::tcp;
using boost::system::error_code;

using std::cout;
using std::endl;

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

	bool Socket::write(const Socket::byte * buff, const int size)
	{
		try {

			boost::asio::io_service io_service;
			tcp::resolver resolver(io_service);
			tcp::resolver::query query(_host, _service);
			tcp::resolver::iterator endpoint_iterator = resolver.resolve(query);

			tcp::socket socket(io_service);
			boost::asio::connect(socket, endpoint_iterator); 
			boost::system::error_code error;

			boost::asio::write(socket, boost::asio::buffer(buff, size), boost::asio::transfer_all(), error);

		} catch (std::exception& e) {

			std::cerr << e.what() << std::endl;
			return false;
		}
		return true;
	}

	class Header {
	public:

		Header()
		: name("object"), 
		type("undefined"), 
		color_rgb("255, 0, 0"), 
		shader("Color Shader"),
		mode("facesWithFrames"),
		// mode("faces"),
		is_polygon(false)
		{}

		void write(std::ostream &os) {
			std::stringstream ss;
			
			ss << "{\n \"type\": \"" << type <<  "\", \"name\": \"" << name << "\"," <<
			" \"color\" : [" << color_rgb << "], " <<
			"\"shader\": \"" << shader << "\", " <<
			"\"mode\": \"" << mode << "\", " <<
			"\"polygon\":" << (is_polygon? "true" : "false") << "\n} ";	

			std::string header = ss.str();		
			const size_t len = header.size(); 
			os.write((const char *)&len, sizeof(len));	
			os.write(header.c_str(), header.size());
		}

		std::string name;
		std::string type;
		std::string color_rgb;
		std::string shader;
		std::string mode;
		bool is_polygon;
	};

	static void aux_plot(Header &h, const int dims, const int n_points, const double *points)
	{
		std::stringstream os;
		h.write(os);

		os.write((char *)&dims, sizeof(int));
		os.write((char *)&n_points, sizeof(int));

		for(int i = 0; i < dims; ++i) {
			for(int j = 0; j < n_points; ++j) {
				os.write((const char *)&points[j * dims + i], sizeof(double));
			}
		}

		if(dims < 3) {
			for(int i = dims; i < 3; ++i) {
				for(int j = 0; j < n_points; ++j) {
					double zero = 0;
					os.write((const char *)&zero, sizeof(double));
				}
			}
		}

		Socket socket(host, port);
		std::string buffer = os.str();

		// std::cout << "\n" << buffer << std::endl;
		socket.write(buffer.c_str(), buffer.size());
	}


	void plot_polygon(const int dims, const int n_points, const double *points, const std::string &name)
	{
		Header h;
		h.is_polygon = true;
		h.type 		 = "lines";
		h.name       =  name;
		h.color_rgb  = "0, 0, 0";

		aux_plot(h, dims, n_points, points);
	}

	void plot_lines(const int dims, const int n_points, const double *points, const std::string &name)
	{
		Header h;
		h.is_polygon = false;
		h.type 		 = "lines";
		h.name       = name;
		h.color_rgb  = "0, 255, 0";

		aux_plot(h, dims, n_points, points);
	}

	void plot_points(const int dims, const int n_points, const double *points, const std::string &name)
	{
		Header h;
		h.is_polygon = false;
		h.type 		 = "points";
		h.name       = name;
		h.color_rgb  = "164, 174, 0";

		aux_plot(h, dims, n_points, points);
	}

	void plot_quad_points(const int dims, const std::vector<libMesh::Point> &points, const std::string &name)
	{
		std::vector<double> pts(points.size() * dims);

		for(size_t i = 0; i < points.size(); ++i) {
			for(int d = 0; d < dims; ++d) {
				pts[i * dims + d] = points[i](d);
			}
		}

		plot_points(dims, points.size(), &pts[0], name);

		// for(auto p : pts) {

		// 	std::cout << p << "\n";
		// }

		// std::cout << std::endl;

	}

	void quiver(const int dims, const int n_points, const double *points, const double *vectors, const std::string &name)
	{
		Header h;
		h.is_polygon = false;
		h.type 		 = "lines";
		h.name       =  name;
		h.color_rgb  = "0, 0, 255";

		std::stringstream os;
		h.write(os);

		const int n_points_2 = n_points * 2;

		os.write((char *)&dims, sizeof(int));
		os.write((char *)&n_points_2, sizeof(int));

		for(int i = 0; i < dims; ++i) {
			for(int j = 0; j < n_points; ++j) {
				os.write((const char *)&points[j * dims + i], sizeof(double));
				const double displaced = points[j * dims + i] + vectors[j * dims + i];
				os.write((const char *)&displaced, sizeof(double));
			}
		}

		if(dims < 3) {
			for(int i = dims; i < 3; ++i) {
				for(int j = 0; j < n_points; ++j) {
					double zero = 0;
					os.write((const char *)&zero, sizeof(double));
					os.write((const char *)&zero, sizeof(double));
				}
			}
		}

		Socket socket(host, port);
		std::string buffer = os.str();

		// std::cout << "\n" << buffer << std::endl;
		socket.write(buffer.c_str(), buffer.size());
	}

	template<class It>
	void stringfy(const It &begin, const It &end, std::ostream &os)
	{
		for(It it = begin; it != end; ++it) {
			os << std::to_string(*it); 

			if(it + 1 != end) {
				os << ", ";
			}
		}
	}

	void plot_polygon_mesh(const libMesh::MeshBase &mesh, const std::string &name)
	{
		libMesh::DenseMatrix<libMesh::Real> polygon;
		for(auto e_it = mesh.active_local_elements_begin(); e_it != mesh.active_local_elements_end(); ++e_it) {
			make_polygon(**e_it, polygon);
			plot_polygon(polygon.n(), polygon.m(), &polygon.get_values()[0], name + "/" + std::to_string((*e_it)->id()));
		}
	}

	void plot_mesh_aux(const libMesh::MeshBase &mesh,const std::string &shader, const std::string &name, std::ostream &os)
	{
		using namespace libMesh;

		std::vector<int> el_ptr, el_index;
		std::vector<double> points;

		el_ptr.resize(mesh.n_active_local_elem() + 1);
		el_index.reserve(mesh.n_active_local_elem() * 8);
		points.reserve(mesh.n_local_nodes() * mesh.mesh_dimension());

		el_ptr[0] = 0;

		int index = 0;
		for(auto e_it = mesh.active_local_elements_begin(); e_it != mesh.active_local_elements_end(); ++e_it, ++index) {
			auto &el = **e_it;
			for(uint i = 0; i < el.n_nodes(); ++i) {
				el_index.push_back(el.node_id(i)); 
			}

			el_ptr[index + 1] = el_ptr[index] + el.n_nodes();
		}


		for(auto n_it = mesh.local_nodes_begin(); n_it != mesh.local_nodes_end(); ++n_it) {
			auto &p = (**n_it);

			for(uint i = 0; i < mesh.mesh_dimension(); ++i) {
				points.push_back(p(i));
			}
		}


		Header h;
		h.is_polygon = false;
		h.type 		 = "mesh_json";
		h.name       = name;
		h.color_rgb  = "0, 0, 0";
		h.shader     =  shader;


		
		h.write(os);

		std::stringstream ss;
		ss << "{\n" << "\"dims\" : ";
		ss << std::to_string(mesh.mesh_dimension());
		ss << (",\n\"points\" : [");

		for(int i = 0; i < points.size(); ++i) {
			ss << std::to_string(points[i]);
			if(i != points.size()-1) {
				ss << (",");
			}
		}

		ss << "],\n";

		////////////////////////////////
		ss << "\"elptr\" : [";
		stringfy(el_ptr.begin(), el_ptr.end(), ss);
		ss << ("],\n");
		////////////////////////////////
		ss << "\"elindex\" : [";
		stringfy(el_index.begin(), el_index.end(), ss);
		ss << ("]\n}");
		////////////////////////////////

		std::string json_mesh = ss.str();
		long length = json_mesh.length();
		os.write((char *)&length, sizeof(long));
		os << json_mesh;
	}

	void plot_mesh(const libMesh::MeshBase &mesh, const std::string &name)
	{
		Socket socket(host, port);
		std::stringstream os;
		plot_mesh_aux(mesh, "Gooch Shader", name, os);
		std::string buffer = os.str();
		socket.write(buffer.c_str(), buffer.size());
	}

	void plot_mesh_f(const libMesh::MeshBase &mesh, const double * function, const std::string &name)
	{
		Socket socket(host, port);
		std::stringstream os;
		plot_mesh_aux(mesh, "Function shader gouraud", name, os);
		
		long r = mesh.n_local_nodes();
		long c = 1;
		os.write((Socket::byte *)&r, sizeof(long));
		os.write((Socket::byte *)&c, sizeof(long));
		os.write((Socket::byte *)function, r * c * sizeof(double));

		std::string buffer = os.str();
		socket.write(buffer.c_str(), buffer.size());
	}


	void plot_box(const Box &box, const std::string &name)
	{

		const auto &b_min = box.get_min();
		const auto &b_max = box.get_max();

		int n_coords = (pow(2, box.get_dims()) * box.get_dims());
		std::vector<double> points;

		switch(box.get_dims()) {
			case 2:
			{
				points.resize(n_coords);

				//0
				points[0] = b_min(0);
				points[1] = b_min(1);
				
				//1
				points[2] = b_max(0);
				points[3] = b_min(1);
				
				//2
				points[4] = b_max(0);
				points[5] = b_max(1);
				
				//3
				points[6] = b_min(0);
				points[7] = b_max(1);

				plot_polygon(2, 4, &points[0], name);
				break;
			}
			case 3:
			{
				std::vector<double> p2(3), p3(3), p4(3), p5(3), p6(3), p7(3);
				
				//2
				p2[0] = b_max(0);
				p2[1] = b_min(1);
				p2[2] = b_min(2);

				//3
				p3[0] = b_min(0);
				p3[1] = b_min(1);
				p3[2] = b_max(2);

				//4
				p4[0] = b_max(0);
				p4[1] = b_min(1);
				p4[2] = b_max(2);

				//5
				p5[0] = b_min(0);
				p5[1] = b_max(1);
				p5[2] = b_min(2);

				//6
				p6[0] = b_max(0);
				p6[1] = b_max(1);
				p6[2] = b_min(2);

				//7
				p7[0] = b_min(0);
				p7[1] = b_max(1);
				p7[2] = b_max(2);

				points.reserve(n_coords * 2);
				
				//lower square
				points.insert(points.end(), b_min.get_values().begin(), b_min.get_values().end());
				points.insert(points.end(), p2.begin(), p2.end());

				points.insert(points.end(), b_min.get_values().begin(), b_min.get_values().end());
				points.insert(points.end(), p3.begin(), p3.end());

				points.insert(points.end(), p3.begin(), p3.end());
				points.insert(points.end(), p4.begin(), p4.end());

				points.insert(points.end(), p2.begin(), p2.end());
				points.insert(points.end(), p4.begin(), p4.end());

				//vertical axis
				points.insert(points.end(), b_min.get_values().begin(), b_min.get_values().end());
				points.insert(points.end(), p5.begin(), p5.end());

				points.insert(points.end(), p2.begin(), p2.end());
				points.insert(points.end(), p6.begin(), p6.end());

				points.insert(points.end(), p3.begin(), p3.end());
				points.insert(points.end(), p7.begin(), p7.end());

				points.insert(points.end(), p4.begin(), p4.end());
				points.insert(points.end(), b_max.get_values().begin(), b_max.get_values().end());

				//upper square
				points.insert(points.end(), b_max.get_values().begin(), b_max.get_values().end());
				points.insert(points.end(), p6.begin(), p6.end());

				points.insert(points.end(), b_max.get_values().begin(), b_max.get_values().end());
				points.insert(points.end(), p7.begin(), p7.end());

				points.insert(points.end(), p7.begin(), p7.end());
				points.insert(points.end(), p5.begin(), p5.end());

				points.insert(points.end(), p5.begin(), p5.end());
				points.insert(points.end(), p6.begin(), p6.end());

				plot_lines(3, 24, &points[0], name);
				break;
			}
			default:
			{
				break;
			}
		}
	}


	void plot_scaled_normal_field(const libMesh::MeshBase &mesh,
		const DVectord &normals,
		const DVectord &scale,
		const std::string &name)
	{

		using namespace libMesh;
		int mesh_dim = mesh.mesh_dimension();


		DenseVector<double> local_normal;
		DenseVector<Real> local_scale;

		std::vector<double> all_points, all_normals;

		std::vector<double> point(mesh_dim, 0.);
		for(auto n_it = mesh.active_nodes_begin(); n_it != mesh.active_nodes_end(); ++n_it) {
			Node &n = **n_it;

			std::vector<dof_id_type> node_dof_ids;

			for(int d = 0; d < mesh_dim; ++d) {
				auto dof_id = n.dof_number(0, d, 0);
				node_dof_ids.push_back(dof_id);

				point[d] = n(d);
			}

			get_vector(normals, node_dof_ids, local_normal);
			get_vector(scale,   node_dof_ids, local_scale);

			for(int d = 0; d < mesh_dim; ++d) {
				local_normal(d) *= local_scale(0);
			}

			if(local_normal.l2_norm() < 1e-16) continue;

			all_points.insert(all_points.end(),
				point.begin(),
				point.end());

			all_normals.insert(all_normals.end(),
				local_normal.get_values().begin(),
				local_normal.get_values().end());
		}

		if(all_points.empty()) {
			return;
		}

		quiver(mesh_dim,
			all_points.size()/mesh_dim,
			&all_points[0],
			&all_normals[0],
			name);
	}
}

#else 
//empty functions
namespace utopia {

	void plot_polygon(const int, const int, const double *, const std::string &)
	{
		std::cerr << "[Warning] plot function not implemented, make sure to have a proper boost installation" << std::endl;
	}

	void plot_lines(const int, const int, const double *, const std::string &)
	{
		std::cerr << "[Warning] plot function not implemented, make sure to have a proper boost installation" << std::endl;
	}

	void plot_points(const int, const int, const double *, const std::string &)
	{
		std::cerr << "[Warning] plot function not implemented, make sure to have a proper boost installation" << std::endl;
	}

	void quiver(const int, const int, const double *, const double *vectors, const std::string &)
	{
		std::cerr << "[Warning] plot function not implemented, make sure to have a proper boost installation" << std::endl;
	}

	void plot_quad_points(const int, const std::vector<libMesh::Point> &, const std::string &)
	{
		std::cerr << "[Warning] plot function not implemented, make sure to have a proper boost installation" << std::endl;
	}

	void plot_mesh(const libMesh::MeshBase &, const std::string &)
	{
		std::cerr << "[Warning] plot function not implemented, make sure to have a proper boost installation" << std::endl;
	}

	void plot_mesh_f(const libMesh::MeshBase &, const double * function, const std::string &)
	{
		std::cerr << "[Warning] plot function not implemented, make sure to have a proper boost installation" << std::endl;
	}

	void plot_polygon_mesh(const libMesh::MeshBase &, const std::string &)
	{
		std::cerr << "[Warning] plot function not implemented, make sure to have a proper boost installation" << std::endl;
	}

	void plot_box(const Box &, const std::string &)
	{
		std::cerr << "[Warning] plot function not implemented, make sure to have a proper boost installation" << std::endl;
	}


	void plot_scaled_normal_field(
		const libMesh::MeshBase &,
		const DVectord &,
		const DVectord &,
		const std::string &)
	{
		std::cerr << "[Warning] plot function not implemented, make sure to have a proper boost installation" << std::endl;
	}
}

#endif //WITH_BOOST
