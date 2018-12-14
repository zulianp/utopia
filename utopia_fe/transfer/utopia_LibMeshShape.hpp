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
        
        LibMeshShape(const libMesh::Elem &elem, const libMesh::FEType type, const bool use_newton = true)
        : elem_(elem), type_(type), q_(Dim), max_iter_(100), tol_(1e-14), use_newton_(use_newton), verbose_(false)
        {
            init();
        }
        
        inline void verbose(const bool val)
        {
            verbose_ = val;
        }

        inline bool intersect(const Ray<Scalar, Dim> &ray,
                              Scalar &t) override
        {
        	if(use_newton_) {
        		return intersect_newton(ray, t);
        	} else {
            	return intersect_gradient_descent(ray, t);
            }
        }

        inline void ref(Vector &ref_point) const
        {
            Read<Vectord> r_(x_ref_);

            for(int d = 0; d < Dim-1; ++d) {
                ref_point[d] = x_ref_.get(d);
            }

            ref_point[Dim-1] = 0.;
        }


        inline bool make_quadrature(
            const Vector &ray_dir,
            const std::vector<Vector>  &composite_q_points,
            const std::vector<Scalar>  &composite_q_weights,
            QMortar &q
        )
        {
            std::vector<Scalar> gap;
            return make_quadrature(
                ray_dir,
                composite_q_points,
                composite_q_weights,
                q,
                gap);
        }

        inline bool make_quadrature(
            const Vector &ray_dir,
            const std::vector<Vector>  &composite_q_points,
            const std::vector<Scalar>  &composite_q_weights,
            QMortar &q,
            std::vector<Scalar> &gap
        )
        {
            Ray<Scalar, Dim> ray;
            ray.dir = ray_dir;

            const std::size_t n_qp = composite_q_weights.size();
            q.resize(n_qp);
            gap.resize(n_qp);

            Vector ref_point;
            for(std::size_t i = 0; i < n_qp; ++i) {
                q.get_weights()[i] = composite_q_weights[i];

                ray.o = composite_q_points[i];

                Scalar t = 0.;

                if(!intersect(ray, t)) {
                    std::cerr << "[Error] now what?!" << std::endl;
                    assert(false);
                    return false;
                }

                gap[i] = t;

                ref(ref_point);

                for(int d = 0; d < Dim; ++d) {
                    q.get_points()[i](d) = ref_point[d];
                } 
            }

            return true;
        }

    private:
        const libMesh::Elem &elem_;
        libMesh::FEType type_;
        std::unique_ptr<libMesh::FEBase> fe_;
        QMortar q_;
        
        //mats and vecs
        Vectord g_, g_1_, x_, x_ref_, p_, n_, r_, h_, Jtn_, v_, u_;
        Matrixd A_, J_, H_, H_fe_, JtJ_;
        
        LUDecomposition<Matrixd, Vectord> solver_;
        
        int max_iter_;
        Scalar tol_;
        bool use_newton_;
        bool verbose_;
        
        inline void get_point(Vectord &x) const
        {
            Write<Vectord> w_(x);
            const auto &xyz = fe_->get_xyz();
            for(int i = 0; i < Dim; ++i) {
                x.set(i, xyz[0](i));
            }
        }
        
        inline void get_hessian(const Vectord &v, Matrixd &H) const
        {        
        	Read<Vectord> r_(v);
        	Write<Matrixd> w_(H);

            const auto &d2d2x = fe_->get_d2xyzdxi2();

            auto val = 0.;

            for(int i = 0; i < Dim; ++i) {
            	val += v.get(i) * d2d2x[0](i);
            }

            H.set(0, 0, val);
         
        	if(Dim > 2) {
        		const auto &d2dxdy = fe_->get_d2xyzdeta2();
        		const auto &d2d2y  = fe_->get_d2xyzdxideta();

        		auto val01 = 0., val11 = 0.;

        		for(int i = 0; i < Dim; ++i) {
        			val01 += v.get(i) * d2dxdy[0](i);
        			val11 += v.get(i) * d2d2y[0](i);
        		}

        		H.set(0, 1, val01);
        		H.set(1, 0, val01);
        		H.set(1, 1, val11);
        	}
        }
        
        inline void get_jacobian(Matrixd &J) const
        {
            Write<Matrixd> w_(J);
            
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
        }
       
        void init()
        {
            fe_ = libMesh::FEBase::build(elem_.dim(), type_);
            fe_->get_xyz();
            fe_->get_dxyzdxi();
            
            if(Dim > 2) {
                fe_->get_dxyzdeta();
            }

            if(use_newton_) {
	            fe_->get_d2xyzdxi2();
	         
            	if(Dim > 2) {
            		fe_->get_d2xyzdeta2();
            		// fe_->get_d2xyzdzeta2();
            		fe_->get_d2xyzdxideta();
            	}

            	// fe_->get_d2xyzdxidzeta();
            	// fe_->get_d2xyzdetadzeta();
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
            A_ = zeros(Dim -1, Dim - 1);
            
            H_fe_ = zeros(Dim-1, Dim-1);
        }
        
        bool intersect_gradient_descent(const Ray<Scalar, Dim> &ray,
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
                
                JtJ_ = transpose(J_) * J_;
                Jtn_ = transpose(J_) * n_;
                
                r_ = p_ + t * n_;
                
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

                if(verbose_) { std::cout << "gradient_descent: iter: " << i << " " << norm_g << std::endl; }
                
                if(norm_g_prev < norm_g) {
                    std::cout << "diverged" << std::endl;
                    return false;
                }
                
                if(norm_g <= tol_) {
                    return true;
                }
                
                norm_g_prev = norm_g;
            }
            
            return false;
        }

        bool intersect_newton(const Ray<Scalar, Dim> &ray,
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

                if(verbose_) {
                    disp("x:");
                    disp(x_);
                    disp("----------");
                }
                
                JtJ_ = transpose(J_) * J_;
                Jtn_ = transpose(J_) * n_;
                
                r_ = p_ + t * n_;
                
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

                const Scalar norm_g = norm2(g_);

                if(verbose_) { std::cout << "newton: iter: " << i << " " << norm_g << std::endl; }
                
                if(norm_g_prev < norm_g) {
                    std::cout << "diverged" << std::endl;
                    return false;
                }
                
                if(norm_g <= tol_) {
                    return true;
                }

                get_hessian(r_, A_);

                {
                    Write<Matrixd> w_H(H_);
                    Read<Matrixd> r_A(A_);
                    Read<Vectord> r_Jtn(Jtn_);

                    for(int d1 = 0; d1 < Dim-1; ++d1) {
                        for(int d2 = 0; d2 < Dim-1; ++d2) {
                            H_.set(d1, d2, A_.get(d1, d2) + JtJ_.get(d1, d2));
                        }

                        H_.set(d1, Dim-1, -Jtn_.get(d1));
                        H_.set(Dim-1, d1, -Jtn_.get(d1));
                    }

                    H_.set(Dim-1, Dim-1, nn);
                }

                solver_.solve(H_, g_, h_); 
                u_ -= h_;

                if(verbose_) {

                    disp("----------");
                    
                    disp("J:");
                    disp(J_);
                    disp("----------");

                    disp("JtJ:");
                    disp(JtJ_);
                    disp("----------");

                    disp("A:");
                    disp(A_);
                    disp("----------");
                  
                    disp("H:");
                    disp(H_);
                    disp("----------");
                   
                    disp("g:");
                    disp(g_);
                    disp("----------");
                   
                    disp("h:");
                    disp(h_);
                    disp("----------");
                   
                    disp("u:");
                    disp(u_);
                    disp("----------");



                }
                
                {
                    //update ref coordinates
                    Read<Vectord> r_u(u_);
                    Write<Vectord> w_x_ref(x_ref_);
                    
                    for(int d = 0; d < Dim - 1; ++d) {
                        x_ref_.set(d, u_.get(d));
                    }
                    
                    //update ray intersection
                    t = u_.get(Dim - 1);
                }
                
                if(verbose_) {
                   disp("x_ref:");
                   disp(x_ref_);
                   disp("----------");
                }
                
                norm_g_prev = norm_g;
            }
            
            return false;
        }

        
    protected:
        virtual void reinit(const Vectord &x_ref)
        {
            Read<Vectord> r_(x_ref);
            
            for(int i = 0; i < Dim - 1; ++i) {
                q_.get_points()[0](i) = x_ref.get(i);
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
