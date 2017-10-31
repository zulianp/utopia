#ifndef UTOPIA_HOMEMADE_MESH_HPP
#define UTOPIA_HOMEMADE_MESH_HPP

#include "utopia_Base.hpp"
#include "utopia_homemade_FEForwardDeclarations.hpp"

#include <vector>
#include <memory>

namespace utopia {

	class Mesh {
	public:
		class Impl;

		void make_triangle();
		int element_order(const int element_index) const;
		void node_indices(const int elem, std::vector<int> &index);
		int n_dims() const;
		Mesh();
		~Mesh();


		void *mesh_impl_ptr() const;


	private:
		//memory
		std::vector<int> el_ptr;
		std::vector<int> el_index;
		std::vector<int> el_type;
		std::vector<int> meta;
		std::vector<double> points;


		std::unique_ptr<Impl> impl_ptr;
	};

}

#endif //UTOPIA_HOMEMADE_MESH_HPP
