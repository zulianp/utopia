#ifndef MFEM_L2P_HASH_GRID_HPP
#define MFEM_L2P_HASH_GRID_HPP 

#include <vector>
#include "Box.hpp"
#include "libmesh/mesh.h"

namespace utopia {
	inline bool is_simplex(const int type)
	{
		return type == static_cast<int>(libMesh::TRI3) || type == static_cast<int>(libMesh::TRI6) ||
			   type == static_cast<int>(libMesh::TET4) || type == static_cast<int>(libMesh::TET10);
	}

	class HashGrid {
	public:
		long hash(const libMesh::DenseVector<libMesh::Real> &point) const;
		long hash(const std::vector<long> &coord) const;
		void hash_range(const libMesh::DenseVector<libMesh::Real> &min, const libMesh::DenseVector<libMesh::Real> &max, std::vector<long> &hashes);
		HashGrid(const Box &box, const std::vector<int> &dims);

		inline long n_cells() const
		{
			return n_cells_;
		}

		void print(std::ostream &os) const;

	private:
		Box box_;
		libMesh::DenseVector<libMesh::Real> range_;
		std::vector<int> dims_;
		long n_cells_;
	};

	bool fix_normal_orientation(const libMesh::Elem &elem, int side, libMesh::Point &n);
	void enlarge_box_from_side(const int dim, const libMesh::Elem &side, Box &box, const libMesh::Real blow_up);

	bool hash_grid_detect_intersection(const Box &src_box,
									   const std::vector<int> &dims,
									   const std::vector<Box> &src_boxes, 
									   const std::vector<Box> &dest_boxes, 
									   std::vector<int> &pairs,
									   const libMesh::Real tol = 0.0
									   );

	void build_boxes(const libMesh::MeshBase &mesh, std::vector<Box> &element_boxes, Box &mesh_box);
	void boundary_build_boxes(const libMesh::MeshBase &mesh, std::vector<Box> &element_boxes, Box &mesh_box, std::vector<long> &map, const libMesh::Real blow_up);
	bool hash_grid_detect_intersections(const libMesh::MeshBase &src, const libMesh::MeshBase &dest, std::vector<int> &pairs);
	bool boundary_hash_grid_detect_intersections(const libMesh::MeshBase &mesh, std::vector<int> &pairs, const libMesh::Real blow_up);
}

#endif //MFEM_L2P_HASH_GRID_HPP
