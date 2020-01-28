#ifndef UTOPIA_PROJECTED_GAUSS_SEIDEL_QR_HPP
#define UTOPIA_PROJECTED_GAUSS_SEIDEL_QR_HPP

#include "utopia_ForwardDeclarations.hpp"
#include "utopia_BoxConstraints.hpp"
#include "utopia_IterativeSolver.hpp"
#include "utopia_Smoother.hpp"
#include "utopia_QPSolver.hpp"
#include "utopia_RowView.hpp"
#include "utopia_ProjectedGaussSeidel.hpp"
#include <cmath>

namespace utopia {
    //slow and innefficient implementation just for testing
    template<class Matrix, class Vector>
    class ProjectedGaussSeidelQR : public ProjectedGaussSeidel<Matrix, Vector>
    {
    public:

        DEF_UTOPIA_SCALAR(Matrix)

        ProjectedGaussSeidelQR()
        : use_line_search_(false), n_local_sweeps_(3)
        {
            ProjectedGaussSeidel<Matrix, Vector>::use_line_search(false);

        }

        ProjectedGaussSeidelQR(const ProjectedGaussSeidelQR &) = default;

        inline ProjectedGaussSeidelQR * clone() const override
        {
            auto ptr = new ProjectedGaussSeidelQR(*this);
            ptr->set_box_constraints(this->get_box_constraints());
            return ptr;
        }


        void read(Input &in) override
        {
            QPSolver<Matrix, Vector>::read(in);
            Smoother<Matrix, Vector>::read(in);

            in.get("use_line_search", use_line_search_);
            in.get("n_local_sweeps", n_local_sweeps_);
        }

        void set_R(const Matrix &R){
            R_ = R;
        }

        virtual bool smooth(const Vector &b, Vector &x) override
        {
            const Matrix &A = *this->get_operator();

            // init(A);
            SizeType it = 0;
            SizeType n_sweeps = this->sweeps();
            if(this->has_bound()) {
                std::cout<<"-- constrained step.... "<<std::endl;
               step(A, b, x);
            } else {
                std::cout<<"-- unconstrained step.... "<<std::endl;
                //while(unconstrained_step(A, b, x) && it++ < n_sweeps) {}
            }
            return it == SizeType(this->sweeps() - 1);
        } 

        bool apply(const Vector &b, Vector &x) override
        {
            if(this->verbose())
                this->init_solver("utopia ProjectedGaussSeidel", {" it. ", "|| u - u_old ||"});

            const Matrix &A = *this->get_operator();

            x_old = x;
            bool converged = false;
            const SizeType check_s_norm_each = 1;

            int iteration = 0;
            while(!converged) {
                if(this->has_bound()) {
                    step(A, b, x);
                } else {
                    this->unconstrained_step(A, b, x);
                }

                if(iteration % check_s_norm_each == 0) {
                    const Scalar diff = norm2(x_old - x);

                    // if(this->verbose()) {
                        PrintInfo::print_iter_status({static_cast<Scalar>(iteration), diff});
                    // }

                    converged = this->check_convergence(iteration, 1, 1, diff);
                }

                ++iteration;

                if(converged) break;

                x_old = x;
            }
            return converged;
        }        


        void print_usage(std::ostream &os) const override
        {
            QPSolver<Matrix, Vector>::print_usage(os);
            Smoother<Matrix, Vector>::print_usage(os);

            this->print_param_usage(os, "use_line_search", "bool", "Determines if line-search should be used.", "true");
            this->print_param_usage(os, "n_local_sweeps", "int", "Number of local sweeps.", "3");
        }

        bool step(const Matrix &A, const Vector &b, Vector &x) override
        {
            r = b - A * x;

            std::cout<<"doing step ... "<< std::endl;

            d = diag(A);
            d_inv = 1./d;            

            inactive_set_ = local_values(local_size(b).get(0), 0.0);
            //localize gap function for correction
            g = this->get_upper_bound() - R_*x;
            l = this->get_lower_bound() - R_*x;
            
            c *= 0.;
            
            Scalar g_i, l_i;
            
            Range rr = row_range(A);
            {
                Write<Vector> w_a(inactive_set_);
                ReadAndWrite<Vector> rw_c(c);
                Read<Vector> r_d_inv(d_inv), r_g(g), r_l(l);
                Read<Matrix> r_A(A);
                Read<Vector> r_x(x);
                Read<Vector> r_r(r);
                Read<Matrix> r_R(R_);           
                SizeType n_rows = local_size(R_).get(0); 
                //SizeType n_cols = local_size(R_).get(1); 

                //std::cout<<"n_cols: "<< this->n_local_sweeps() << "   n_rows: "<< n_rows << "  \n";

                // for(SizeType il = 0; il < 1; il++) 
                // {

                    for(auto i = rr.begin(); i != rr.end(); ++i) 
                    {
                        RowView<const Matrix> row_view(A, i);
                        decltype(i) n_values = row_view.n_values();

                        Scalar s = r.get(i);
                        
                        for(auto index = 0; index < n_values; ++index) 
                        {
                            const decltype(i) j = row_view.col(index);
                            const auto a_ij = row_view.get(index);

                            if(rr.inside(j) && j != i) 
                            {
                                s -= a_ij * c.get(j);
                            }
                        }

                        //update correction
//                        x.set(i, d_inv.get(i)*(b.get(i) - s) + x.get(i)) ;
                        c.set(i, d_inv.get(i)*s) ;

                        if  (i < n_rows)
                        {
                            RowView<const Matrix> row_viewR(R_, i);     
                            decltype(i) nnz_R = row_viewR.n_values();
                            //std::cout << "nnz_R: " << nnz_R << std::endl;

                            Scalar r_sum_c = 0.0;

                            Scalar r_ii = 0.0;
                            for(auto index = 0; index < nnz_R; ++index) 
                            {
                                r_ii = 0.0; 
                                const decltype(i) j = row_viewR.col(index);
                                const auto r_ij = row_viewR.get(index);

                                if(j < i)
                                {
                                    r_sum_c += r_ij * c.get(j);
                                }
                                else if(i == j)
                                {
                                    r_ii = r_ij;
                                }
                            }
    
                            // std::cout <<" r_sum_c" << r_sum_c << std::endl;
                            if (r_ii > 0)
                            {
                                g_i = (g.get(i) - r_sum_c)/r_ii; // g.get(i)*invR_ii
                                l_i = (l.get(i) - r_sum_c)/r_ii; // g.get(i)*invR_ii
                            }
                            else if (r_ii < 0)
                            {
                                l_i = (g.get(i) - r_sum_c)/r_ii; // g.get(i)*invR_ii
                                g_i = (l.get(i) - r_sum_c)/r_ii; // g.get(i)*invR_ii
                            }

                            //update correction
                            // std::cout << "l_i:" << l_i << "  g_i:" << g_i << "  c_i:" << c.get(i) << std::endl;
                            if (( g_i <= c.get(i) ) || (c.get(i) <= l_i))
                            {
                                c.set(i, std::max(std::min( c.get(i), g_i), l_i));
                                inactive_set_.set(i, 0.0);
                                //std::cout << "GS: activeset id:" << i << std::endl;
                            }
                            else
                            {
                                inactive_set_.set(i, 1.0);
                            }
                        }//std::cout << "row_sum:" << r_sum << std::endl;
                    }

                    // inactive_set_ *= 0.;
                    // for (auto i = rr.begin(); i != rr.end(); ++i)
                    // {
                    //     // std::cout<<"l.get(i): "<< l.get(i) <<std::endl;
                    //     // std::cout<<"g.get(i): "<< g.get(i) <<std::endl;
                    //     // std::cout<<"x.get(i): "<< x.get(i) <<std::endl;

                    //     if(( l.get(i) < x.get(i) ) && (x.get(i) < g.get(i)))
                    //     {
                    //         inactive_set_.set(i, 1.);
                    //     }
                    //     else
                    //     {
                    //         std::cout << "GS: activeset id:" << i << std::endl;
                    //     }
                    // }

                    std::cout<<"------------------------------------------------------- "<<std::endl;

                    //disp(inactive_set_);
                }

            // }



            x += c;
            return true;
        }

        void init(const Matrix &A)
        {
            d = diag(A);
            d_inv = 1./d;
            c = local_zeros(local_size(A).get(0));
            inactive_set_ = local_zeros(local_size(c));
        }


        virtual void update(const std::shared_ptr<const Matrix> &op) override
        {
            IterativeSolver<Matrix, Vector>::update(op);
            init(*op);
        }


        void use_line_search(const bool val)
        {
            use_line_search_ = val;
        }

        inline SizeType n_local_sweeps() const
        {
            return n_local_sweeps_;
        }

        inline void n_local_sweeps(const SizeType n_local_sweeps)
        {
            n_local_sweeps_ = n_local_sweeps;
        }

        inline void use_symmetric_sweep(const bool use_symmetric_sweep)
        {
            use_symmetric_sweep_ = use_symmetric_sweep;
        }

        const Vector& get_inactive_set()
        {
            return inactive_set_;
        }



        
    private:
        bool use_line_search_;
        bool use_symmetric_sweep_;
        SizeType n_local_sweeps_;

        Vector r, d, g,l ,c, d_inv, x_old, descent_dir;
        Vector inactive_set_;
        Vector is_c_;
        Matrix R_;
    };
}

#endif //UTOPIA_PROJECTED_GAUSS_SEIDEL_QR_HPP
