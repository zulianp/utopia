#ifndef UTOPIA_LIBMESH_SHAPE_HPP
#define UTOPIA_LIBMESH_SHAPE_HPP

#include "utopia_Project.hpp"
#include "utopia_MortarAssemble.hpp"

#include <libmesh/fe_type.h>
#include <libmesh/fe.h>

namespace utopia {
	template<typename Scalar, int Dim>
	class LibMeshShape : public Shape<Scalar, Dim> {
	public:
		using Vector = utopia::Vector<Scalar, Dim>;
		virtual ~LibMeshShape() {}

		LibMeshShape(const libMesh::Elem &elem, const libMesh::FEType type)
		: elem_(elem), type_(type), q_(1), max_iter_(10), tol_(1e-14)
		{
			init();
		}

		bool intersect(
			const Ray<Scalar, Dim> &ray,
			Scalar &t) override
		{
			{
				Write<Vectord> w_n(p_), w_p(n_);

				for(int i = 0; i < Dim; ++i) {
					p_.set(i, ray.origin[i]);
					n_.set(i, ray.normal[i]);
				}
			}

			//valid initial guess
			x_ref_.set(0.);


			Scalar g_2 = 0.;
			Scalar nn = dot(n_, n_);

			for(int i = 0; i < max_iter_; ++i) {
				reinit(x_ref_);

				//get main fe quantities
				get_point(x_);
				get_jacobian(J_);
				
				JtJ_ = transpose(J_) * J_;
				Jtn_ = transpose(J_) * n_;

				r_ = p_ + t * n_;

				g_1_ = -(transpose(J_) * (x_ - r_));
				g_2  = dot(n_, r_ - x_);

				//build gradient
				{
					Write<Vectord> w_g(g_);

					for(int d = 0; d < Dim-1; ++d) {
						g_.set(d, g_1_.get(d));
					}

					g_.set(Dim-1, g_2);

				}

				v_ = x_ - r_;
				get_hessian(v_, H_fe_);

				H_fe_ += JtJ_;
				
				//build hessian
				{
					Write<Matrixd> w_H(H_);
					Read<Vectord> r_Jtn(Jtn_);

					//H_11
					each_read(H_fe_, [&](const SizeType i, const SizeType j, const Scalar value) {
						H_.set(i, j, value);
					});

					//Off diag blocks
					for(int d = 0; d < Dim-1; ++d) {
						const Scalar val = - Jtn_.get(d);
						
						H_.set(0, d, val);
						H_.set(d, 0, val);
					}

					//H_dd
					H_.set(Dim-1, Dim-1, nn);
				}

				solver_.solve(H_, g_, h_);

				//update solution
				{
					Read<Vectord> r_h(h_);
					Write<Vectord> w_x_ref(x_ref_);

					const Scalar alpha = -1.;

					for(int d = 0; d < Dim-1; ++d) {
						x_ref.add(d, alpha * h_.get(d));
					}

					t += alpha * h_.get(Dim-1);
				}

				Scalar norm_h = norm2(h_);

				if(norm2 <= tol_) {
					return true;
				}
			}

			return false;
		}

	private:
		const libMesh::Elem &elem_;
		libMesh::FEType type_;
		std::unique_ptr<libMesh::FEBase> fe_;
		QMortar q_;

		//mats and vecs
		Vectord g_, g_1_, x_, x_ref_, p_, n_, r_, h_, Jtn_, v_;
		Matrixd J_, H_, H_fe_, JtJ_;

		LUDecomposition<Matrixd, Vectord> solver_;

		int max_iter_;
		Scalar tol_;

		inline void get_point(Vectord &x) const
		{
			// x.set(0.);

			Write<Vectord> w_(x);

			const auto &xyz = fe_->get_xyz();
			const auto &phi = fe_->get_phi();
			const auto n_shape_fun = phi.size();

			// for(std::size_t k = 0; k < n_shape_fun; ++k) {
				for(int i = 0; i < Dim; ++i) {
					x.set(i, xyz[0](i));// * phi[k][0]));
				}
			// }
		}

		inline void get_hessian(const Vectord &v, Matrixd &H) const
		{
			H.set(0.);

			Read<Vectord> r_(v);
			Write<Matrixd> w_(H);
			//FIXME ref-element hessian?
			const auto &d2phi = fe_->get_d2phi();
			const auto n_shape_fun = d2phi.size();

			for(std::size_t k = 0; k < n_shape_fun; ++k) {
				const auto dot_prod = /* element coords */ xyz[0](0) * v.get(0);

				for(int i = 1; i < Dim; ++ i) {
					dot_prod += /* element coords */ xyz[0](i) * v.get(i);
				}

				for(int i = 0; i < Dim; ++i) {
					for(int j = 0; j < Dim - 1; ++j) {
						auto val = dot_prod * d2phi[k][0](i, j);
						J.add(i, j, val);
					}
				}
			}

		}

		inline void get_jacobian(Matrixd &J) const
		{
			J.set(0.);

			Write<Matrixd> w_(J);
			//FIXME ref-element gradient?
			const auto &dphi = fe_->get_dphi();
			const auto n_shape_fun = dphi.size();

			for(std::size_t k = 0; k < n_shape_fun; ++k) {
				for(int i = 0; i < Dim; ++i) {
					for(int j = 0; j < Dim - 1; ++j) {
						auto val = /* element coords */ xyz[0](i) * dphi[k][0][j];
						J.add(i, j, val);
					}
				}
			}
		}


		void init()
		{
			fe_ = libMesh::FEBase::build(elem_.dim(), type_);
			fe_->get_xyz();
			fe_->get_dphi();
			fe_->get_d2phi();
			fe_->get_phi();

			q_.get_weights()[0] = 1.;

			g_ = zeros(Dim);
			x_ = zeros(Dim);
			x_ref_ = zeros(Dim);
			p_ = zeros(Dim);
			n_ = zeros(Dim);
			h_ = zeros(Dim);

			g_1_ = zeros(Dim-1);

			J_ = zeros(Dim, Dim-1);
			H_ = zeros(Dim, Dim);

			H_fe_ = zeros(Dim-1, Dim-1);
		}

	protected:
		virtual void reinit(const Vectord &x_ref)
		{
			Read<Vectord> r_(x_ref);

			for(int i = 0; i < Dim; ++i) {
				q_.get_points()(i) = x_ref.get(i);
			}

			fe_->attach_quadrature_rule(&q_);
			fe_->reinit(*elem_);
		}

	};

	template<typename Scalar, int Dim>
	class LibMeshSideShape final : public Shape<Scalar, Dim> {
	public:
		LibMeshSideShape(const libMesh::Elem &elem, const libMesh::FEType type, const int side)
		: elem_(elem), type_(type), side_(side)
		{}

		bool intersect(
			const Ray<Scalar, Dim> &ray,
			Scalar &t) override
		{
			return false;
		}

	private:
		const libMesh::Elem &elem_;
		libMesh::FEType type_;
		int side_;
	};



}


#endif //UTOPIA_LIBMESH_SHAPE_HPP
