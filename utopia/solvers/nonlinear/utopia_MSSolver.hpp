#ifndef UTOPIA_MS_SOLVER_HPP
#define UTOPIA_MS_SOLVER_HPP

#include "utopia_Core.hpp"
#include "utopia_LinearSolver.hpp"
#include "utopia_Function.hpp"
#include "utopia_NonLinearSolver.hpp"
#include "utopia_LS_Strategy.hpp"
#include "utopia_LinearSolverInterfaces.hpp"

#include <iomanip>
#include <limits>
#include <memory>


namespace utopia
{
    template<class Matrix, class Vector, int Backend = Traits<Vector>::Backend>
    class MSSolver final : public NonLinearSolver<Matrix, Vector>
    {
        typedef UTOPIA_SCALAR(Vector)    Scalar;
        typedef UTOPIA_SIZE_TYPE(Vector) SizeType;
        typedef typename NonLinearSolver<Matrix, Vector>::Solver Solver;
        using LinearSolverT = LinearSolver<Matrix, Vector>;
        
    public:
        
        enum NormType {
            L2_NORM = 0,
            A_NORM = 1,
            A_SQUARED_NORM = 2,
            N_NORM_TYPES
        };
        
        class Normed {
        public:
            virtual ~Normed() {}
            virtual void gradient(Function<Matrix, Vector> &function, const Vector &x, Vector &gradient) = 0;
            virtual Scalar norm(const Vector &x) const = 0;
            virtual Scalar dot(const Vector &left, const Vector &right) const = 0;
            virtual void update(const std::shared_ptr<Matrix> &mat) {}
            virtual void transform_gradient(Vector &in, Vector &out) = 0;
            virtual bool needs_hessian() const { return true; }
        };
        
        class L2Normed final : public Normed {
        public:
            void gradient(Function<Matrix, Vector> &function, const Vector &x, Vector &gradient) override
            {
                function.gradient(x, gradient);
            }
            
            Scalar norm(const Vector &x) const override
            {
                return norm2(x);
            }
            
            Scalar dot(const Vector &left, const Vector &right) const override
            {
                return utopia::dot(left, right);
            }
            
            void transform_gradient(Vector &in, Vector &out) override
            {
                out = in;
            }
            
            bool needs_hessian() const override { return false; }
        };
        
        class ANormed final : public Normed {
        public:
            ANormed(const std::shared_ptr<LinearSolverT> &linear_solver)
            : M_inv(linear_solver)
            {}
            
            void update(const std::shared_ptr<Matrix> &mat) override
            {
                M = mat;
                M_inv->update(mat);
            }
            
            void gradient(Function<Matrix, Vector> &function, const Vector &x, Vector &gradient) override
            {
                function.gradient(x, temp);
                M_inv->apply(temp, gradient);
            }
            
            void transform_gradient(Vector &in, Vector &out) override
            {
                M_inv->apply(in, out);
            }
            
            Scalar norm(const Vector &x) const override
            {
                return std::sqrt(utopia::dot(x, (*M) * x));
            }
            
            Scalar dot(const Vector &left, const Vector &right) const override
            {
                return std::sqrt(utopia::dot(left, (*M) * right));
            }
            
            std::shared_ptr<Matrix> M;
            std::shared_ptr<LinearSolverT> M_inv;
            Vector temp;
        };
        
        class ASquaredNormed final : public Normed {
        public:
            ASquaredNormed(const std::shared_ptr<LinearSolverT> &linear_solver)
            : M_inv(linear_solver)
            {}
            
            void update(const std::shared_ptr<Matrix> &mat) override
            {
                M = mat;
                if(!M2) {
                    M2 = std::shared_ptr<Matrix>();
                }
                
                *M2 = transpose(*M) * (*M);
                M_inv->update(M2);
            }
            
            void gradient(Function<Matrix, Vector> &function, const Vector &x, Vector &gradient) override
            {
                function.gradient(x, temp);
                M_inv->apply(temp, gradient);
            }
            
            void transform_gradient(Vector &in, Vector &out) override
            {
                M_inv->apply(in, out);
            }
            
            Scalar norm(const Vector &x) const override
            {
                return utopia::norm2((*M) * x);
            }
            
            Scalar dot(const Vector &left, const Vector &right) const override
            {
                return utopia::dot((*M) * left, (*M) * right);
            }
            
            std::shared_ptr<Matrix> M;
            std::shared_ptr<Matrix> M2;
            std::shared_ptr<LinearSolverT> M_inv;
            Vector temp;
        };
        
        MSSolver(const std::shared_ptr <Solver> &linear_solver):
        NonLinearSolver<Matrix, Vector>(linear_solver),
        delta_(0.3),
        delta_prime_(0.35),
        epsilon0_(1.),
        T1_([this](const Scalar &epsilon) -> Scalar { return epsilon/epsilon0_; }),
        T2_([](const Scalar &epsilon) -> Scalar { return 0.35 * epsilon; }),
        G_([](const Scalar &epsilon, const Scalar &epsilon0) -> Scalar { return epsilon0; }),
        norm_type_(L2_NORM),
        convex_hull_n_gradients_(2)
        {
            normed_.resize(N_NORM_TYPES);
            normed_[L2_NORM] = std::make_shared<L2Normed>();
            normed_[A_NORM]  = std::make_shared<ANormed>(std::shared_ptr<LinearSolverT>(linear_solver->clone()));
            normed_[A_SQUARED_NORM] = std::make_shared<ASquaredNormed>(std::shared_ptr<LinearSolverT>(linear_solver->clone()));
        }
        
        // bool convex_hull_minmia
        
        bool B(Function<Matrix, Vector> &fun,
               Normed &normed,
               const Vector &left_in,
               const Vector &right_in,
               const Vector &dir,
               const Scalar &tol,
               const Scalar &step_size, //not needed
               const Scalar &val_left_in,
               const Scalar &val_right_in,
               Vector &h_g)
        {
            Scalar val_left = val_left_in;
            Scalar val_right = val_right_in;
            
            Vector leftpoint = left_in;
            Vector midpoint  = right_in;
            Vector rightpoint;
            
            std::cout << "step_size: " << step_size << std::endl;

            //we compute the directional derivative
            Vector g;
            fun.gradient(midpoint, g);
            Scalar d = dot(g, dir);
            Scalar prev_d = d;
            
            bool success = true;
            SizeType it = 0;
            while(d > tol) {
                Scalar f_val = 0.;
                fun.value(midpoint, f_val);
                
                if(empty(rightpoint)) {
                    rightpoint = right_in;
                    midpoint  = 0.5 * (left_in + right_in);
                } else {
                    auto left_diff = f_val - val_left;
                    auto right_diff = val_right - f_val;
                    
                    std::cout << left_diff << "<=" << right_diff << std::endl;

                    if(left_diff <= right_diff) {
                        leftpoint = midpoint;
                        midpoint = 0.5 * (midpoint + rightpoint);
                        val_left = f_val;
                    } else {
                        rightpoint = midpoint;
                        midpoint = 0.5 * (leftpoint + midpoint);
                        val_right = f_val;
                    }
                }
                
                fun.gradient(midpoint, g);
                d = dot(g, dir);
               
                ++it;
                
                if(it > 100) {
                    success = false;
                    break;
                }
            }
            
            normed.transform_gradient(g, h_g);
            return success;
        }
        
        void line_search(Function<Matrix, Vector> &fun, const Vector &x, const Vector &dir, Scalar &alpha)
        {
            
        }
        
        bool solve(Function<Matrix, Vector> &fun, Vector &x) override
        {
            using namespace utopia;
            
            gradients_.clear();
            
            Vector y, h_g, b_g;
            
            Scalar g_norm = 1.;
            SizeType it = 0;
            
            bool converged = false;
            
            this->init_solver("MSSolver", {" it. ", "|| g ||"});
            
            Scalar epsilonk = epsilon0_;
            auto normed = normed_[norm_type_];
            std::shared_ptr<Matrix> H = std::make_shared<Matrix>();
            
            while(!converged)
            {
                if(normed->needs_hessian()) {
                    fun.hessian(x, *H);
                    normed->update(H);
                }
                
                normed->gradient(fun, x, h_g);
                g_norm = normed->norm(h_g);
                epsilonk = G_(g_norm, epsilonk);
                
                // // print iteration status on every iteration
                if(this->verbose_) {
                    PrintInfo::print_iter_status(it, {g_norm});
                }
                
                // // check convergence and print interation info
                converged = this->check_convergence(it, g_norm, 1, 1);
                
                if(converged) {
                    break;
                }
                
                Scalar val_x = 0.;
                fun.value(x, val_x);
                
                bool stop = false;
                while(!stop) {
                    gradients_.clear();
                    gradients_.push_back(h_g);
                    
                    auto a_norm = g_norm;
                    auto a_norm2 = a_norm * a_norm;
                    
                    while(a_norm > T1_(epsilonk)) {
                        
                        const auto &dir = gradients_.back();
                        Scalar step_size = (epsilonk/g_norm);
                        assert(step_size > 0.);
                        
                        y = x - step_size * dir;
                        Scalar val_y = 0.;
                        fun.value(y, val_y);
                        const Scalar diff_val = val_y - val_x;
                        
                        if(diff_val < (-delta_ * g_norm * epsilonk)) {
                            Scalar alpha = epsilonk/g_norm;
                            line_search(fun, x, -dir, alpha);
                            
                            assert(alpha > 0.);
                            
                            x -= alpha * dir;
                            stop = true;
                            break;
                            
                        } else {
                            //use B to compute new gradient
                            
                            if(!B(fun,
                                  *normed,
                                  x, //left
                                  y, //right
                                  gradients_.back(), //dir
                                  delta_prime_ * a_norm2, //tol
                                  step_size,
                                  val_x, //val_left
                                  val_y, //val_right
                                  b_g
                                  )) {
                                assert(false);
                                return false;
                            }
                            
                            if(convex_hull_n_gradients_ == 2) {
                                //case with only 2 vectors
                                
                                //compute smallest element of the set of gradients
                                const auto &first_g = gradients_.back();
                                const auto n_bg = normed->norm(b_g);
                                const auto n_bg2 = n_bg * n_bg;
                                const auto dot_bg_fg = normed->dot(b_g, first_g);
                                
                                const auto lambda = (n_bg2 - dot_bg_fg)/(n_bg2 + a_norm2 - 2. * dot_bg_fg);
                                
                                //add small elements to list
                                gradients_.back() = (lambda * first_g + (1.-lambda) * b_g);
                            } else {
                                assert(false);
                                //case with more here:
                                //....
                            }
                            
                            a_norm = normed->norm(gradients_.back());
                            a_norm2 = a_norm * a_norm;
                            
                            if(epsilonk < this->atol() && a_norm < this->atol()) {
                                stop = true;
                                converged = true;
                                break;
                            }
                        }
                    }
                    
                    if(!stop) {
                        epsilonk = T2_(epsilonk);
                    }
                }
                
                
                it++;
            }
            
            this->print_statistics(it);
            return true;
        }
        
        void read(Input &in) override
        {
            NonLinearSolver<Matrix, Vector>::read(in);
            // in.get("dumping", delta_);
        }
        
    private:
        Scalar delta_, delta_prime_, epsilon0_;   /*!< Dumping parameter. */
        std::function<Scalar(const Scalar &)> T1_, T2_;
        std::function<Scalar(const Scalar &, const Scalar &)> G_;
        // std::shared_ptr<LSStrategy> ls_strategy_;     /*!< Strategy used in order to obtain step \f$ \delta_k \f$ */
        NormType norm_type_;
        
        std::vector<std::shared_ptr<Normed>> normed_;
        std::vector<Vector> gradients_;
        SizeType convex_hull_n_gradients_;
    };
    
}
#endif //UTOPIA_MS_SOLVER_HPP
