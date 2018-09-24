#ifndef UTOPIA_UI_MESH_HPP
#define UTOPIA_UI_MESH_HPP

#include "utopia_ui.hpp"
#include "utopia_libmesh.hpp"

#include "libmesh/mesh_refinement.h"

#include <memory>

namespace utopia {
	template<class Mesh>
	class UIMesh {};// final : public Serializable { };

	template<>
	class UIMesh<libMesh::DistributedMesh> : public Serializable {
	public:
		template<class... Args>
		UIMesh(Args &&...args)
		: mesh_(std::make_shared<libMesh::DistributedMesh>(std::forward<Args...>(args...))) 
		{}

		void read(InputStream &is) override {
			std::string mesh_type = "square";
			std::string path = "";
	
			int refinements = 0;
			
			double span[3] = { 0., 0., 0. };

			double min_coords[3] = {0., 0., 0.};
			double max_coords[3] = {1., 1., 1.};

			int n[3] = {5, 5, 5};

			double scale = 1.;

			std::string elem_type = "quad";

			is.read("type", mesh_type);
			is.read("elem-type", elem_type);
			is.read("order", order);
			is.read("path", path);

			is.read("refinements", refinements);
			
			is.read("span-x", span[0]);
			is.read("span-y", span[1]);
			is.read("span-z", span[2]);
			
			is.read("min-x", min_coords[0]);
			is.read("min-y", min_coords[1]);
			is.read("min-z", min_coords[2]);

			is.read("max-x", max_coords[0]);
			is.read("max-y", max_coords[1]);
			is.read("max-z", max_coords[2]);

			is.read("n-x", n[0]);
			is.read("n-y", n[1]);
			is.read("n-z", n[2]);

			is.read("scale", scale);


			if(mesh_type == "file") {
				mesh_->read(path);
			} else if(mesh_type == "square") {
				libMesh::MeshTools::Generation::build_square(*mesh_,
					n[0], n[1],
					min_coords[0], max_coords[0],
					min_coords[1], max_coords[1],
					get_type(elem_type, order, 2)
					);
			} else if(mesh_type == "cube") {
				libMesh::MeshTools::Generation::build_cube(*mesh_,
					n[0], n[1], n[2],
					min_coords[0], max_coords[0],
					min_coords[1], max_coords[1],
					min_coords[2], max_coords[2],
					get_type(elem_type, order, 3)
					);
			} else if(mesh_type == "aabb") {
				libMesh::DistributedMesh temp_mesh(mesh_->comm());
				temp_mesh.read(path);

				auto bb = bounding_box(temp_mesh);

				if(temp_mesh.spatial_dimension() == 3) {
					libMesh::MeshTools::Generation::build_cube(
						*mesh_,
						n[0], n[1], n[2],
						bb.min()(0) - span[0], bb.max()(0) + span[0],
						bb.min()(1) - span[1], bb.max()(1) + span[1],
						bb.min()(2) - span[2], bb.max()(2) + span[2],
						get_type(elem_type, order, 3)
						);

				} else {
					libMesh::MeshTools::Generation::build_square(
						*mesh_,
						n[0], n[1],
						bb.min()(0) - span[0], bb.max()(0) + span[0],
						bb.min()(1) - span[1], bb.max()(1) + span[1],
						get_type(elem_type, order, 2)
						);
				}
			}

			scale_mesh(scale, *mesh_);

			refine(refinements, *mesh_);

			if(mesh_type == "file" && order == 2) {
				mesh_->all_second_order();
			}
		}

		inline libMesh::DistributedMesh &mesh()
		{
			assert(mesh_);
			return *mesh_;
		}

		inline std::shared_ptr<libMesh::DistributedMesh> mesh_ptr()
		{
			return mesh_;
		}

		inline bool empty() const {
			return bool(mesh_);
		}

	private:
		int order = 1;
		std::shared_ptr<libMesh::DistributedMesh> mesh_;

		////////////////////////////////////////////////

		libMesh::ElemType get_type(
			const std::string &elem_type,
			const int order,
			const int dim) const
		{
		    if(dim == 3) {
		        libMesh::ElemType type = libMesh::HEX8;

		        if(elem_type == "tet") {
		            type = libMesh::TET4;
		        }

		        if(order == 2) {
		            type = libMesh::HEX20;

		            if(elem_type == "tet") {
		                type = libMesh::TET10;
		            }
		        }

		        return type;

		    } else if(dim == 2) {
		        libMesh::ElemType type = libMesh::QUAD4;

		        if(elem_type == "tri") {
		            type = libMesh::TRI3;
		        }
		        
		        if(order == 2) {
		            type = libMesh::QUAD8;

		            if(elem_type == "tri") {
		                type = libMesh::TRI6;
		            }
		        }

		        return type;
		    }

		    return libMesh::TRI3;
		}

		static void refine(const int n_refs, libMesh::MeshBase &mesh)
		{
		    if(n_refs <= 0) return;
		    
		    libMesh::MeshRefinement mesh_refinement(mesh);
		    mesh_refinement.make_flags_parallel_consistent();
		    mesh_refinement.uniformly_refine(n_refs);
		}

		static void scale_mesh(const double &scale_factor, libMesh::MeshBase &mesh)
		{
			if(scale_factor == 1.) return;
			assert(scale_factor > 0.);

			// for(auto it = mesh.local_nodes_begin(); it != mesh.local_nodes_end(); ++it) {
			for(auto it = mesh.nodes_begin(); it != mesh.nodes_end(); ++it) {

				for(int i = 0; i < LIBMESH_DIM; ++i) {
					(**it)(i) *= scale_factor;
				}
			}
		}
	};
}


#endif //UTOPIA_UI_MESH_HPP
