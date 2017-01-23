
#include <assert.h>
#include "HashGrid.hpp"
#include <libmesh/elem.h>
#include "utopia.hpp"

#include "MortarAssemble.hpp"

namespace utopia {
	void HashGrid::print(std::ostream &os) const
	{
		os << "box:\n";
		box_.print(os);

		os << "range:\n";
		range_.print(os);

		os << "dims: ";
		for(auto d : dims_) {
			os << d << " ";
		}

		os << "\n";
		os << "n_cells: " << n_cells_ << "\n";
	}

	long HashGrid::hash(const libMesh::DenseVector<libMesh::Real> &point) const
	{
		double x = (point(0) - box_.get_min(0))/range_(0);
		long result = floor(x * dims_[0]);

		long totalDim = dims_[0];

		for(int i = 1; i < range_.size(); ++i) {
			result *= dims_[i];

			x = (point(i) - box_.get_min(i))/range_(i);
			result += floor(x * dims_[i]);
			totalDim *= dims_[i];
		}

		if(result >= totalDim || result < 0) {
			printf("error -> %d\n", result);
		}

		return result;
	}

	long HashGrid::hash(const std::vector<long> &coord) const
	{
		long result   = coord[0];
		long totalDim = dims_[0];

		for(int i = 1; i < range_.size(); ++i) {
			result *= dims_[i];
			result += coord[i];
			totalDim *= dims_[i];
		}

		return result;
	}

	void HashGrid::hash_range(const libMesh::DenseVector<libMesh::Real> &min, const libMesh::DenseVector<libMesh::Real> &max, std::vector<long> &hashes)
	{
		const int dim = min.size();
		std::vector<long> imin(dim), imax(dim);

			//generate tensor indices
		for(int i = 0; i < dim; ++i) {
			double x = (min(i) - box_.get_min(i))/range_(i);
			imin[i] = floor(x * dims_[i]);
		}

		for(int i = 0; i < dim; ++i) {
			double x = (max(i) - box_.get_min(i))/range_(i);
			imax[i] = floor(x * dims_[i]);
		}

		std::vector<long> offsets(dim);
		for(int i = 0; i < dim; ++i) {
			offsets[i] = imax[i] - imin[i];
		}

			//FIXME make more general for greater dim
		if(dim == 1) {
			std::vector<long> coord(1);
			for(int i = imin[0]; i <= imax[0]; ++i) {
				coord[0] = i;
				hashes.push_back(hash(coord)); 
			}

		} else if(dim == 2) {
			std::vector<long> coord(2);
			for(int i = imin[0]; i <= imax[0]; ++i) {
				for(int j = imin[1]; j <= imax[1]; ++j) {
					coord[0] = i;
					coord[1] = j;
					hashes.push_back(hash(coord)); 
				}
			}
		} else if(dim == 3) {
			std::vector<long> coord(3);
			for(int i = imin[0]; i <= imax[0]; ++i) {
				for(int j = imin[1]; j <= imax[1]; ++j) {
					for(int k = imin[2]; k <= imax[2]; ++k) {
						coord[0] = i;
						coord[1] = j;
						coord[2] = k;
						hashes.push_back(hash(coord)); 
					}
				}
			}
		} else {
			assert(false && "dim > 3 not supported yet!");
		}

		assert(!hashes.empty());
	}

	HashGrid::HashGrid(const Box &box, const std::vector<int> &dims)
	: box_(box), dims_(dims), n_cells_(1)
	{
		box_.enlarge(1e-8);

		range_  = box_.get_max();
		range_ -= box_.get_min();

		for(int i = 0; i < dims_.size(); ++i) {
			n_cells_ *= dims_[i];
		}
	}

	void build_boxes(const libMesh::MeshBase &mesh, std::vector<Box> &element_boxes, Box &mesh_box)
	{
		const int dim = mesh.mesh_dimension();
		element_boxes.resize(mesh.n_active_local_elem());
		mesh_box.reset(dim);

		auto e_begin = mesh.active_local_elements_begin();
		auto e_end   = mesh.active_local_elements_end();

		libMesh::DenseMatrix<libMesh::Real> polygon;
		libMesh::Point p;

		long index = 0;
		for(auto e_it = e_begin; e_it != e_end; ++e_it, ++index) {
			const auto &e = **e_it;

			element_boxes[index].reset(dim);

			if(e.has_affine_map()) {
				for(int i = 0; i < e.n_nodes(); ++i) {				
					element_boxes[index] += e.point(i);
					mesh_box += element_boxes[index];
				}
			} else {
				
				// if(e.type() == libMesh::TRI6) {
					utopia::make_polygon(e, polygon);

					for(int i = 0; i < polygon.m(); ++i) {	
						p(0) = polygon(i, 0);
						p(1) = polygon(i, 1);

						element_boxes[index] += p;
						mesh_box += element_boxes[index];
					}

				// } else {
				// 	std::cerr << "[Warning] non-affine element not necessarily bounded by box\n" << std::endl;
				// 	for(int i = 0; i < e.n_nodes(); ++i) {				
				// 		element_boxes[index] += e.point(i);
				// 		mesh_box += element_boxes[index];
				// 	}
				// }
			}
		}
	}

	bool fix_normal_orientation(const libMesh::Elem &elem, int side, libMesh::Point &n)
	{
		using namespace libMesh;
		Point c = elem.centroid();

		auto side_ptr = elem.build_side_ptr(side);
		Point dir = side_ptr->point(0) - c;

		if(dir.contract(n) < 0) {
			n *= -1;
			return true;
		}

		return false;
	}

	void compute_side_normal(const int dim, const libMesh::Elem &side, libMesh::Point &n)
	{
		using namespace libMesh;
		Point o, u, v;

		if(dim == 2) {
			assert(side.n_nodes() == 2);
			o = side.point(0);
			u = side.point(1);
			u -= o;
			n(0) =  u(1);
			n(1) = -u(0);
			
		} else {
			assert(dim == 3);
			o = side.point(0);
			u = side.point(1);
			v = side.point(2);	
			u -= o;
			v -= o;
			n = u.cross(v);
		}

		n *= 1./n.norm();
	}

	void enlarge_box_from_side(const int dim, const libMesh::Elem &side, Box &box, const libMesh::Real blow_up)
	{

		libMesh::Point c, n, p;

		compute_side_normal(dim, side, n);

		// c = side.centroid();

		for(uint i = 0; i < side.n_nodes(); ++i) {
			c = side.point(i);

			n *= blow_up;
			p = c;
			p += n;

			box += p;
			p = c;
			n *= 0.01;
			p -= n;
			box += p;
		}
	}

	void boundary_build_boxes(const libMesh::MeshBase &mesh, std::vector<Box> &element_boxes, Box &mesh_box, std::vector<long> &map, const libMesh::Real blow_up)
	{
		const int dim = mesh.mesh_dimension();
		element_boxes.reserve(mesh.n_active_local_elem());
		map.reserve(element_boxes.capacity());
		mesh_box.reset(dim);

		auto e_begin = mesh.active_local_elements_begin();
		auto e_end   = mesh.active_local_elements_end();

		libMesh::Point o, u, v, n, c, p;
		long map_index = 0;
		for(auto e_it = e_begin; e_it != e_end; ++e_it, ++map_index) {
			const auto &e = **e_it;
			if(!e.on_boundary()) { 
				continue;
			}

			const long index = map.size();

			element_boxes.emplace_back();
			element_boxes[index].reset(dim);

			for(uint side = 0; side < e.n_sides(); ++side) {
				if(e.neighbor_ptr(side) != nullptr) continue;

				auto side_ptr = e.build_side_ptr(side);	
				compute_side_normal(dim, *side_ptr, n);
				
				if(fix_normal_orientation(e, side, n)) {
					std::cerr << "[Warning] face with wrong orientation detected\n" << std::endl;
				}

				assert( n.contract(side_ptr->centroid()-e.centroid()) > 0 );

				for(uint i = 0; i < side_ptr->n_nodes(); ++i) {
					c = side_ptr->point(i);
				// n.print(); 
				// std::cout << "\n";

					n *= blow_up;
					p = c;
					p += n;

					// c.print();
					// std::cout << "\n";
					
					element_boxes[index] += p;
					p = c;
					n *= 0.1;
					p -= n;
					element_boxes[index] += p;
				}
			}

			mesh_box += element_boxes[index];
			map.push_back(map_index);

			// element_boxes[index].print(std::cout);


		}
	}

	bool are_neighbors(const libMesh::Elem &e_first, const libMesh::Elem &e_second)
	{
		if(e_first.has_neighbor(&e_second)) {
			return true;
		}

		if(is_simplex(e_first.type())) {
			for(uint side = 0; side < e_first.n_sides(); ++side) {
				auto neigh_ptr = e_first.neighbor_ptr(side);
				if(neigh_ptr == nullptr) continue;
				if(e_second.has_neighbor(neigh_ptr)) {
					return true;
				}
			}
		}

		if(is_simplex(e_second.type())) {
			for(uint side = 0; side < e_second.n_sides(); ++side) {
				if(e_second.neighbor_ptr(side) == nullptr) continue;
				if(e_first.has_neighbor(e_second.neighbor_ptr(side))) {
					return true;
				}
			}
		}

		for(uint node_first = 0; node_first < e_first.n_nodes(); ++node_first) {
			for(uint node_second = 0; node_second < e_first.n_nodes(); ++node_second) {
				if(e_first.node_id(node_first) == e_second.node_id(node_second)) {
					return true;
				}
			}
		}

		return false;
	}

	bool boundary_hash_grid_detect_intersections(const libMesh::MeshBase &mesh, std::vector<int> &pairs, const libMesh::Real blow_up)
	{
		utopia::Chrono c;

		bool verbose = true;

		const int dim = mesh.mesh_dimension();
		Box box(dim);

		std::vector<Box> boxes;
		std::vector<long> map;

		c.start();
		boundary_build_boxes(mesh, boxes, box, map, blow_up);
		c.stop();

		if(verbose) {
			std::cout << "build boxes: ";
			c.describe(std::cout);
			box.print(std::cout);
		}

		const int n_x_dim = std::max(1, int(pow(mesh.n_active_local_elem(), 1./dim)));	
		
		std::vector<int> dims(dim);
		for(int i = 0; i < dim; ++i) {
			dims[i] = n_x_dim;
		}

		c.start();
		if(!hash_grid_detect_intersection(box, dims, boxes, boxes, pairs, 1e-16)) {
			return false;
		}

		if(verbose) {
			std::cout << "finding candidates: ";
			c.describe(std::cout);
		}

		//remap indices
		for(auto it = pairs.begin(); it != pairs.end(); /*inside*/) {
		assert(*it < map.size());
		*it++ = map[*it];

		assert(*it < map.size());
		*it++ = map[*it];
	}

		//clean-up indices
	std::vector<int> cleaned_up_pairs;
	cleaned_up_pairs.reserve(pairs.size());

		for(auto it = pairs.begin(); it != pairs.end(); /*inside*/) {
	int first  = *it++;
	int second = *it++;

			//remove self
	if(first == second) continue;

	const auto &e_first  = *mesh.elem(first);
	const auto &e_second = *mesh.elem(second);

	if(are_neighbors(e_first, e_second)) {
				// std::cout << "neigs\n";
		continue;
	}

			// std::cout << first << ", " << second << std::endl;

	cleaned_up_pairs.push_back(first);
	cleaned_up_pairs.push_back(second);
}

if(verbose) {
	std::cout << "n_boxes		 = " << boxes.size() << "\n";
	std::cout << "before_n_pairs = " << (pairs.size()/2) << "\n";
	std::cout << "n_pairs		 = " << (cleaned_up_pairs.size()/2) << "\n";
}

pairs = std::move(cleaned_up_pairs);
return !pairs.empty();
}

bool hash_grid_detect_intersections(const libMesh::MeshBase &src, const libMesh::MeshBase &dest, std::vector<int> &pairs)
{
	const int dim = dest.mesh_dimension();
	Box src_box(dim);
	Box dest_box(dim);

	std::vector<Box> src_boxes;
	std::vector<Box> dest_boxes;

	build_boxes(src,  src_boxes,  src_box);
	build_boxes(dest, dest_boxes, dest_box);

	const int n_x_dim = std::max(1, int(pow(src.n_active_local_elem(), 1./dim)));	
	std::vector<int> dims(dim);
	for(int i = 0; i < dim; ++i) {
		dims[i] = n_x_dim;
	}

	return hash_grid_detect_intersection(src_box, dims, src_boxes, dest_boxes, pairs);
}

bool hash_grid_detect_intersection(const Box &src_box,
	const std::vector<int> &dims,
	const std::vector<Box> &src_boxes, 
	const std::vector<Box> &dest_boxes, 
	std::vector<int> &pairs,
	const libMesh::Real tol
	) 
{
	HashGrid hgrid(src_box, dims);
	std::vector< std::vector<int> > src_table(hgrid.n_cells());

		// hgrid.print(std::cout);

	std::vector<long> hashes;
	for(int i = 0; i < src_boxes.size(); ++i) {
		const Box &b = src_boxes[i];
		hashes.clear();
		hgrid.hash_range(b.get_min(), b.get_max(), hashes);

		for(auto h : hashes) {
			if(h < 0) continue;
			
			src_table[h].push_back(i);
		}
	}

	std::vector<int> candidates;
	for(int i = 0; i < dest_boxes.size(); ++i) {
		const Box &b = dest_boxes[i];

		hashes.clear();
		hgrid.hash_range(b.get_min(), b.get_max(), hashes);
		assert(!hashes.empty());

		candidates.clear();

		for(auto h : hashes) {
			for(int j : src_table[h]) {
				const Box &bj = src_boxes[j];

				if(b.intersects(bj, tol)) {
					candidates.push_back(j);
				}
			}
		}

		std::sort(candidates.begin(), candidates.end());
		auto last = std::unique(candidates.begin(), candidates.end());
		candidates.erase(last, candidates.end());

		for(auto c : candidates) {
			pairs.push_back(c);
			pairs.push_back(i);
		}
	}

	return !pairs.empty();
}
}