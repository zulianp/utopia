#ifndef UTOPIA_GRID_HPP
#define UTOPIA_GRID_HPP


#include "utopia.hpp"

#include "utopia_intersector.hpp"
#include "utopia_make_unique.hpp"
#include "moonolith_static_math.hpp"
// #include "utopia_Socket.hpp"

#include <functional>
#include <array>
#include <vector>

namespace utopia {

	namespace private_ {
		template<int Dim, typename IntT>
		struct GridDofs {};
	}

	//grid defined in unit-cube
	template<int Dim>
	class Grid {
	public:
		using Scalar  = double;
		using Integer = long;
		using SizeT   = std::size_t;

		using Vector  = moonolith::Vector<Scalar, Dim>;
		using Array   = std::array<Integer, Dim>;
		using Mapping = std::function<Vector (const Vector &)>;
		using DOFMapping = std::function<Integer (const Integer &)>;
		using Index   = std::vector<Integer>;

		//////////////////// Points, Nodes, and DOFs //////////////////

		inline Vector point(const Integer &hash) const
		{
			Array index;
			node_index_from_hash(hash, index);
			return point(index);
		}

		inline Vector point(const Array &index) const
		{
			Vector ret;
			for(int i = 0; i < Dim; ++i) {
				assert(index[i] <= dims[i]);
				ret[i] = static_cast<Scalar>(index[i])/(dims[i]);
			}

			return map(ret);
		}

		inline Integer node_hash_from_index(const Array &coord) const
		{
			Integer result = coord[0];

			for(int i = 1; i < Dim; ++i) {
				result *= dims[i] + 1;
				result += coord[i];
			}

			return result;
		}

		inline Integer node_dof_from_hash(const Integer &hash) const
		{
			return dof_map(hash);
		}

		inline Integer node_dof_from_index(const Array &coord) const
		{
			return node_dof_from_hash(node_hash_from_index(coord));
		}

		inline void node_index_from_hash(const Integer hash, Array &coord) const
		{
			Integer current = hash;
			const Integer last = Dim - 1;

			for(Integer i = last; i >= 0; --i) {
				const Integer next = current / (dims[i] + 1);
				coord[i] = current - next * (dims[i] + 1);
				current = next;	
			}

			assert(hash == this->node_hash_from_index(coord));
		}

		template<typename IntT>
		void dofs(const Integer cell_hash, std::vector<IntT> &dof_indices) const
		{
			utopia::private_::GridDofs<Dim, IntT>::get(*this, cell_hash, dof_indices);
		}

		//////////////////// element access ////////////////////////

		inline Array element_index(const Integer hash) const
		{
			Array ret;
			element_index_from_hash(hash, ret);
			return ret;
		}

		inline void element_aabb(const Integer element_hash, Vector &emin, Vector &emax) const
		{
			auto imin = element_index(element_hash);
			auto imax = imin;

			for(auto &i : imax) { ++i; }

			emin = point(imin);
			emax = point(imax);
		}

		inline Integer element_hash(const Vector &y) const
		{
			Vector x = inverse_map(y);
			Integer result = floor(x[0] * dims[0]);

			Integer total_dim = dims[0];

			for(int i = 1; i < Dim; ++i) {
				result *= dims[i];
				result += floor(x[i] * dims[i]);
				total_dim *= dims[i];
			}

			if(result >= total_dim || result < 0) {
				printf("error -> %ld\n", result);
			}

			assert(result < n_elements());
			return result;
		}

		inline Integer element_hash(const Array &coord) const
		{
			return element_hash_from_index(coord);
		}

		//z and y major
		inline Integer element_hash_from_index(const Array &coord) const
		{
			assert(element_is_valid(coord));

			Integer result = coord[0];

			for(int i = 1; i < Dim; i++) {
				result *= dims[i];
				result += coord[i];
			}

			assert(result < n_elements());
			return result;
		}

		inline void element_index_from_hash(const Integer hash, Array &coord) const
		{
			Integer current = hash;
			const Integer last = Dim - 1;

			for(Integer i = last; i >= 0; --i) {
				const Integer next = current / dims[i];
				coord[i] = current - next * dims[i];
				current = next;	
			}

			assert(hash == this->element_hash_from_index(coord));
			assert(element_is_valid(coord));
		}

		void elements_in_range(
			const Vector &min,
			const Vector &max,
			Index &hashes) const
		{
			hashes.clear();

			Array imin, imax;

			auto im_min = inverse_map(min);
			auto im_max = inverse_map(max);

			//generate tensor indices
			for(int i = 0; i < Dim; ++i) {
				const Integer id = floor(im_min[i] * dims[i]);
				imin[i] = std::min( std::max(id, Integer(0)), Integer(dims[i] - 1) );
			}

			for(int i = 0; i < Dim; ++i) {
				const Integer id = floor(im_max[i] * dims[i]);
				imax[i] = std::min( std::max(id, Integer(0)), Integer(dims[i] - 1) );
			}

			assert(element_is_valid(imin));
			assert(element_is_valid(imax));

			//FIXME make more general for greater dim and extract to class
			Array coord;
			if(Dim == 1) {
				for(int i = imin[0]; i <= imax[0]; ++i) {
					coord[0] = i;

					auto h = element_hash(coord); assert(h < n_elements());
					hashes.push_back(h);
				}

			} else if(Dim == 2) {
				for(int i = imin[0]; i <= imax[0]; ++i) {
					for(int j = imin[1]; j <= imax[1]; ++j) {
						coord[0] = i;
						coord[1] = j;

						auto h = element_hash(coord); assert(h < n_elements());
						hashes.push_back(h);
					}
				}
			} else if(Dim == 3) {
				for(int i = imin[0]; i <= imax[0]; ++i) {
					for(int j = imin[1]; j <= imax[1]; ++j) {
						for(int k = imin[2]; k <= imax[2]; ++k) {
							coord[0] = i;
							coord[1] = j;
							coord[2] = k;

							auto h = element_hash(coord); assert(h < n_elements());
							hashes.push_back(h);
						}
					}
				}
			} else {
				assert(false && "dim > 3 not supported yet!");
			}

			assert(!hashes.empty());
		}


		Grid()
		{
			map = [](const Vector &x) -> Vector { return x; };
			inverse_map = map;
			dof_map = [](const Integer &hash) -> Integer { return hash; };
			std::fill(std::begin(dims), std::end(dims), 0);
		}

		inline Integer n_elements(const Integer dim) const
		{
			return dims[dim];
		}

		inline Integer n_elements() const
		{
			Integer ret = 1;

			for(auto d : dims) {
				ret *= d;
			}

			return ret;
		}

		inline bool element_is_valid(const Array &index) const
		{
			for(int i = 0; i < Dim; ++i) {
				if(index[i] < 0 || index[i] >= dims[i]) {
					return false;
				}
			}

			return true;
		}

		inline Integer n_nodes(const Integer dim) const
		{
			return dims[dim] + 1;
		}

		inline Integer n_nodes() const
		{
			Integer ret = 1;

			for(auto d : dims) {
				ret *= (d + 1);
			}

			return ret;
		}

		void describe(std::ostream &os) const
		{
			for(auto d : dims) {
				os << d << " ";
			}

			os << "\n";
		}

		//fields
		Array dims;
		Mapping map;
		Mapping inverse_map;
		DOFMapping dof_map;
	};

	namespace private_ {

		template<typename IntT>
		struct GridDofs<1, IntT> {

			static void get(const Grid<1> &grid, const Grid<1>::Integer cell_hash, std::vector<IntT> &dof_indices)
			{
				dof_indices.resize(2);

				auto imin = grid.element_index(cell_hash);
				auto imax = imin;

				for(auto &i : imax) {++i; }

				dof_indices[0] = grid.node_dof_from_index(imin);
				dof_indices[1] = grid.node_dof_from_index(imax);
			}
		};

		template<typename IntT>
		struct GridDofs<2, IntT> {

			static void get(const Grid<2> &grid, const Grid<2>::Integer cell_hash, std::vector<IntT> &dof_indices)
			{
				dof_indices.resize(4);

				auto imin = grid.element_index(cell_hash);
				auto imax = imin;

				for(auto &i : imax) {++i; }

				dof_indices[0] = grid.node_dof_from_index(imin);
				dof_indices[1] = grid.node_dof_from_index({ imax[0], imin[1] });
				dof_indices[2] = grid.node_dof_from_index(imax);
				dof_indices[3] = grid.node_dof_from_index({ imin[0], imax[1] });
			}
		};

		template<typename IntT>
		struct GridDofs<3, IntT> {
			static void get(const Grid<3> &grid, const Grid<3>::Integer cell_hash, std::vector<IntT> &dof_indices)
			{	
				dof_indices.resize(8);

				auto imin = grid.element_index(cell_hash);
				auto imax = imin;

				for(auto &i : imax) {++i; }

					/*
					* HEX8:   7        6
					*         o--------o
					*        /:       /|
					*       / :      / |
					*    4 /  :   5 /  |
					*     o--------o   |
					*     |   o....|...o 2
					*     |  .3    |  /
					*     | .      | /
					*     |.       |/
					*     o--------o
					*     0        1
					*/

				dof_indices[0] = grid.node_dof_from_index(imin);
				dof_indices[1] = grid.node_dof_from_index({ imax[0], imin[1], imin[2] });
				dof_indices[2] = grid.node_dof_from_index({ imax[0], imax[1], imin[2] });
				dof_indices[3] = grid.node_dof_from_index({ imin[0], imax[1], imin[2] });
				dof_indices[4] = grid.node_dof_from_index({ imin[0], imin[1], imax[2] });
				dof_indices[5] = grid.node_dof_from_index({ imax[0], imin[1], imax[2] });
				dof_indices[6] = grid.node_dof_from_index(imax);
				dof_indices[7] = grid.node_dof_from_index({ imin[0], imax[1], imax[2] });
				
			}
		};
	}

}

#endif //UTOPIA_GRID_HPP
