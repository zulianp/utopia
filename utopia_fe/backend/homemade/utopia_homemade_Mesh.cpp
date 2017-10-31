#include "utopia_homemade_Mesh.hpp"
#include "utopia_intersector.hpp"

namespace utopia {

	class Mesh::Impl {
	public:
		typedef Intersector::Mesh Detail;
		Detail detail;
	};

	int Mesh::n_dims() const
	{
		return impl_ptr->detail.n_dims;
	}

	void Mesh::make_triangle()
	{
		el_ptr.resize(2);
		el_ptr[0] = 0;
		el_ptr[1] = 3;

		el_index.resize(3);
		el_index[0] = 0;
		el_index[1] = 1;
		el_index[2] = 2;

		el_type.resize(1);
		el_type[0] = Intersector::ELEMENT_TYPE_TRIANGLE;

		points.resize(3 * 2);
		points[0] = 0.0;
		points[1] = 0.0;

		points[2] = 1.0;
		points[3] = 0.0;

		points[4] = 0.0;
		points[5] = 1.0;

		impl_ptr->detail.n_elements = 1;
		impl_ptr->detail.n_nodes  = 3;
		impl_ptr->detail.n_dims   = 2;
		impl_ptr->detail.el_ptr   = &el_ptr[0];
		impl_ptr->detail.el_index = &el_index[0];
		impl_ptr->detail.el_type  = &el_type[0];
		impl_ptr->detail.meta     = nullptr;// &meta[0];
		impl_ptr->detail.points   = &points[0];
	}

	int Mesh::element_order(const int element_index) const
	{
		switch(el_type[element_index]) {
			case Intersector::ELEMENT_TYPE_TRIANGLE:
			{
				return 1;
			}
			default: 
			{
				return 0;
			}
		}
	}

	void Mesh::node_indices(const int elem, std::vector<int> &index)
	{
		std::size_t begin = impl_ptr->detail.el_ptr[elem];
		std::size_t end   = impl_ptr->detail.el_ptr[elem+1];
		std::size_t n_nodes = end - begin;
		index.resize(n_nodes);
		std::copy(&impl_ptr->detail.el_index[begin], &impl_ptr->detail.el_index[end], index.begin());
	}

	void *Mesh::mesh_impl_ptr() const
	{
		return static_cast<void *>(&impl_ptr->detail);
	}

	Mesh::Mesh()
	: impl_ptr(new Impl()) { }

	Mesh::~Mesh() {}
}
