#include "utopia_homemade_FE.hpp"
#include "utopia_homemade_Mesh.hpp"

namespace utopia {

	class FE::Impl {
	public:
		Intersector fe_backend;
		FEObject fe;
	};

	FE::FE()
	: impl_ptr(new Impl())
	{ }


	void FE::init(const int current_element, Mesh &mesh, int quadrature_order)
	{
		Intersector::Mesh * m_impl = static_cast<Intersector::Mesh *>(mesh.mesh_impl_ptr());



		// impl_ptr->fe_backend.make_fe_object_from_quad_2(
		// 	*mesh.impl, current_element, 3, quad_points, quad_weights, &fe_);
		impl_ptr->fe_backend.make_fe_object(*m_impl, current_element, quadrature_order, &impl_ptr->fe);

		const auto n_funs = impl_ptr->fe.n_shape_functions;
		assert(n_funs >= 1);

		grad.resize(impl_ptr->fe.n_quad_points);
		fun.resize(grad.size());
		dx.resize(grad.size());

		for(int q = 0; q < impl_ptr->fe.n_quad_points; ++q) {
			grad[q].resize(n_funs);
			fun[q].resize(n_funs);

			dx[q] = impl_ptr->fe.dx[q];

			for(int i = 0; i < n_funs; ++i) {
				grad[q][i] = zeros(impl_ptr->fe.n_dims);

				Write<ElementVector> w(grad[q][i]);
				fun[q][i] = impl_ptr->fe.fun[i][q];
				
				for(int k = 0; k < impl_ptr->fe.n_dims; ++k) {
					grad[q][i].set(k, impl_ptr->fe.grad[i][q * impl_ptr->fe.n_dims + k]);
				}
			}
		}
	}

	int FE::n_shape_functions() const
	{
		return impl_ptr->fe.n_shape_functions;
	}

	FE::~FE() {}

}