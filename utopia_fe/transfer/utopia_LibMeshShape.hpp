#ifndef UTOPIA_LIBMESH_SHAPE_HPP
#define UTOPIA_LIBMESH_SHAPE_HPP

#include "utopia_Project.hpp"
#include "MortarAssemble.hpp"

#include <libmesh/fe_type.h>
#include <libmesh/fe.h>

namespace utopia {
	template<typename Scalar, int Dim>
	class LibMeshShape : public Shape<Scalar, Dim> {
	public:
		using Vector = utopia::Vector<Scalar, Dim>;
		virtual ~LibMeshShape() {}

		LibMeshShape(const libMesh::Elem &elem, const libMesh::FEType type)
		: elem_(elem), type_(type), q_(Dim), max_iter_(100), tol_(1e-14)
		{
			init();
		}

		inline bool intersect(
			const Ray<Scalar, Dim> &ray,
			Scalar &t) override
		{
			return intersect_gradient_descent(ray, t);
		}



	private:
		const libMesh::Elem &elem_;
		libMesh::FEType type_;
		std::unique_ptr<libMesh::FEBase> fe_;
		QMortar q_;

		//mats and vecs
		Vectord g_, g_1_, x_, x_ref_, p_, n_, r_, h_, Jtn_, v_, u_;
		Matrixd J_, H_, H_fe_, JtJ_;

		LUDecomposition<Matrixd, Vectord> solver_;

		int max_iter_;
		Scalar tol_;

		inline void get_point(Vectord &x) const
		{
			// x.set(0.);

			Write<Vectord> w_(x);

			const auto &xyz = fe_->get_xyz();
			// const auto &phi = fe_->get_phi();
			// const auto n_shape_fun = phi.size();

			// for(std::size_t k = 0; k < n_shape_fun; ++k) {
				for(int i = 0; i < Dim; ++i) {
					x.set(i, xyz[0](i));// * phi[k][0]));
				}
			// }
		}

		// inline void get_hessian(const Vectord &v, Matrixd &H) const
		// {
		// 	H.set(0.);

		// 	Read<Vectord> r_(v);
		// 	Write<Matrixd> w_(H);
		// 	//FIXME ref-element hessian?
		// 	const auto &d2phi = fe_->get_d2phi();
		// 	const auto n_shape_fun = d2phi.size();

		// 	for(std::size_t k = 0; k < n_shape_fun; ++k) {
		// 		const auto dot_prod = /* element coords */ xyz[0](0) * v.get(0);

		// 		for(int i = 1; i < Dim; ++ i) {
		// 			dot_prod += /* element coords */ xyz[0](i) * v.get(i);
		// 		}

		// 		for(int i = 0; i < Dim; ++i) {
		// 			for(int j = 0; j < Dim - 1; ++j) {
		// 				auto val = dot_prod * d2phi[k][0](i, j);
		// 				J.add(i, j, val);
		// 			}
		// 		}
		// 	}

		// }

		inline void get_jacobian(Matrixd &J) const
		{
			// J.set(0.);

			Write<Matrixd> w_(J);
			//FIXME ref-element gradient?


			// const auto n_shape_fun = jacs.size();

			// for(std::size_t k = 0; k < n_shape_fun; ++k) {

			// for(int j = 0; j < Dim - 1; ++j) {
				const auto &ddxi = fe_->get_dxyzdxi();
				for(int i = 0; i < Dim; ++i) {
					auto val = ddxi[0](i);
					J.set(i, 0, val);
				}

				if(Dim > 2) {
					const auto &ddeta = fe_->get_dxyzdeta();
					for(int i = 0; i < Dim; ++i) {
						auto val = ddeta[0](i);
						J.set(i, 1, val);
					}
				}
			// }
		}


		void init()
		{
			fe_ = libMesh::FEBase::build(elem_.dim(), type_);
			fe_->get_xyz();
			fe_->get_dxyzdxi();

			if(Dim > 2) {
				fe_->get_dxyzdeta();
			}

			q_.resize(1);
			q_.get_weights()[0] = 1.;

			g_ = zeros(Dim);
			x_ = zeros(Dim);
			x_ref_ = zeros(Dim - 1);
			p_ = zeros(Dim);
			n_ = zeros(Dim);
			h_ = zeros(Dim);
			u_ = zeros(Dim);

			g_1_ = zeros(Dim-1);

			J_ = zeros(Dim, Dim-1);
			H_ = zeros(Dim, Dim);

			H_fe_ = zeros(Dim-1, Dim-1);
		}

		bool intersect_gradient_descent(
			const Ray<Scalar, Dim> &ray,
			Scalar &t)
		{
			{
				Write<Vectord> w_n(p_), w_p(n_);

				for(int i = 0; i < Dim; ++i) {
					p_.set(i, ray.o[i]);
					n_.set(i, ray.dir[i]);
				}
			}

			//valid initial guess
			x_ref_.set(0.);
			u_.set(0.);

			Scalar g_2 = 0.;
			Scalar nn = dot(n_, n_);

			auto norm_g_prev = 10000.;

			for(int i = 0; i < max_iter_; ++i) {
				reinit(x_ref_);

				//get main fe quantities
				get_point(x_);
				get_jacobian(J_);

				// std::cout << "current fe: \n";
				// std::cout << "point: " << std::endl;
				// disp(x_);

				// std::cout << "Jacobian: " << std::endl;
				// disp(J_);
				
				JtJ_ = transpose(J_) * J_;
				Jtn_ = transpose(J_) * n_;

				r_ = p_ + t * n_;

				// std::cout << "r: " << std::endl;
				// disp(r_);

				g_1_ = (transpose(J_) * (x_ - r_));
				g_2  = dot(n_, r_ - x_);

				//build overall gradient
				{
					Write<Vectord> w_g(g_);

					for(int d = 0; d < Dim-1; ++d) {
						g_.set(d, g_1_.get(d));
					}

					g_.set(Dim-1, g_2);
				}

				// g_ *= 2.;
				u_ -= g_;

				{
					
					//update ref coordinates
					Read<Vectord> r_u(u_);
					Write<Vectord> w_x_ref(x_ref_);

					for(int d = 0; d < Dim-1; ++d) {
						x_ref_.set(d, u_.get(d));
					}

					//update ray intersection
					t = u_.get(Dim - 1);
				}

				const Scalar norm_g = norm2(g_);

				if(norm_g_prev < norm_g) {
					std::cout << "diverged" << std::endl;
					return false;
				}

				// std::cout << "iter: " << i << " : " << norm_g << std::endl;
				// std::cout << "current \n";
				// disp(u_);

				if(norm_g <= tol_) {
					return true;
				}

				norm_g_prev = norm_g;
			}

			return false;
		}



		// bool intersect_newton(
		// 	const Ray<Scalar, Dim> &ray,
		// 	Scalar &t)
		// {
		// 	{
		// 		Write<Vectord> w_n(p_), w_p(n_);

		// 		for(int i = 0; i < Dim; ++i) {
		// 			p_.set(i, ray.origin[i]);
		// 			n_.set(i, ray.normal[i]);
		// 		}
		// 	}

		// 	//valid initial guess
		// 	x_ref_.set(0.);


		// 	Scalar g_2 = 0.;
		// 	Scalar nn = dot(n_, n_);

		// 	for(int i = 0; i < max_iter_; ++i) {
		// 		reinit(x_ref_);

		// 		//get main fe quantities
		// 		get_point(x_);
		// 		get_jacobian(J_);
				
		// 		JtJ_ = transpose(J_) * J_;
		// 		Jtn_ = transpose(J_) * n_;

		// 		r_ = p_ + t * n_;

		// 		g_1_ = -(transpose(J_) * (x_ - r_));
		// 		g_2  = dot(n_, r_ - x_);

		// 		//build gradient
		// 		{
		// 			Write<Vectord> w_g(g_);

		// 			for(int d = 0; d < Dim-1; ++d) {
		// 				g_.set(d, g_1_.get(d));
		// 			}

		// 			g_.set(Dim-1, g_2);

		// 		}

		// 		v_ = x_ - r_;
		// 		get_hessian(v_, H_fe_);

		// 		H_fe_ += JtJ_;
				
		// 		//build hessian
		// 		{
		// 			Write<Matrixd> w_H(H_);
		// 			Read<Vectord> r_Jtn(Jtn_);

		// 			//H_11
		// 			each_read(H_fe_, [&](const SizeType i, const SizeType j, const Scalar value) {
		// 				H_.set(i, j, value);
		// 			});

		// 			//Off diag blocks
		// 			for(int d = 0; d < Dim-1; ++d) {
		// 				const Scalar val = - Jtn_.get(d);
						
		// 				H_.set(0, d, val);
		// 				H_.set(d, 0, val);
		// 			}

		// 			//H_dd
		// 			H_.set(Dim-1, Dim-1, nn);
		// 		}

		// 		solver_.solve(H_, g_, h_);

		// 		//update solution
		// 		{
		// 			Read<Vectord> r_h(h_);
		// 			Write<Vectord> w_x_ref(x_ref_);

		// 			const Scalar alpha = -1.;

		// 			for(int d = 0; d < Dim-1; ++d) {
		// 				x_ref.add(d, alpha * h_.get(d));
		// 			}

		// 			t += alpha * h_.get(Dim-1);
		// 		}

		// 		Scalar norm_h = norm2(h_);

		// 		if(norm2 <= tol_) {
		// 			return true;
		// 		}
		// 	}

		// 	return false;
		// }

	protected:
		virtual void reinit(const Vectord &x_ref)
		{
			Read<Vectord> r_(x_ref);

			for(int i = 0; i < Dim - 1; ++i) {
				q_.get_points()[i] = x_ref.get(i);
			}

			fe_->attach_quadrature_rule(&q_);
			fe_->reinit(&elem_);
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
