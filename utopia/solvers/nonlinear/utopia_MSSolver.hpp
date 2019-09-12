#ifndef UTOPIA_MS_SOLVER_HPP
#define UTOPIA_MS_SOLVER_HPP

#include "utopia_Core.hpp"
#include "utopia_LinearSolver.hpp"
#include "utopia_Function.hpp"
#include "utopia_NewtonBase.hpp"
#include "utopia_LS_Strategy.hpp"
#include "utopia_LinearSolverInterfaces.hpp"
#include "utopia_QPSolver.hpp"

#include <iomanip>
#include <limits>
#include <memory>


namespace utopia
{

    template<class Matrix, class Vector>
    class HilbertFunction {
    public:
        using Scalar = UTOPIA_SCALAR(Vector);

        virtual ~HilbertFunction() {}
        virtual void gradient(Function<Matrix, Vector> &function, const Vector &x, Vector &gradient) = 0;
        virtual Scalar norm(const Vector &x) const = 0;
        virtual Scalar dot(const Vector &left, const Vector &right) const = 0;
        virtual void update(const std::shared_ptr<Matrix> & /*mat*/) {}
        virtual void transform_gradient(Vector &in, Vector &out) = 0;
        virtual bool needs_hessian() const { return true; }
    };

    template<class GlobalMatrix, class GlobalVector>
    class IMSConvexHullSolver {
    public:
        using Scalar = UTOPIA_SCALAR(GlobalVector);
        virtual ~IMSConvexHullSolver() {}

        virtual void solve(HilbertFunction<GlobalMatrix, GlobalVector> &normed,
                   const Scalar a_norm2,
                   std::vector<GlobalVector> &gradients,
                   GlobalVector &in_out_b_g) = 0;
    };

    template<class GlobalMatrix, class GlobalVector, class Matrix, class Vector>
    class MSConvexHullSolver final : public IMSConvexHullSolver<GlobalMatrix, GlobalVector> {
    public:
        using Scalar = UTOPIA_SCALAR(Vector);

        void solve(HilbertFunction<GlobalMatrix, GlobalVector> &normed,
                   const Scalar a_norm2,
                   std::vector<GlobalVector> &gradients,
                   GlobalVector &in_out_b_g) override
        {
            const std::size_t n_gradients = gradients.size();
            const std::size_t n_gp1 = n_gradients + 1;
            const std::size_t n_gp2 = n_gradients + 2;
            const std::size_t n_2gp2 = 2 * n_gradients + 2;
            const auto m_g = n_gradients * 2 + 3;

            //Create Matrix (A e Id// e^T 0 0// 0 0 0), where the symmetric matrix A= (A_ij )is defined by A_ij = <gradients[i],gradients[j]> with gradients[n_gradients]:= b_g, Id is the Identity and e^T=(1,...,1)
            Matrix C = sparse(m_g , m_g , m_g); // C well be an element of the generalized gradient \partial F (\lambda). Note that the first n+2 rows are fixed.
            Scalar valb_gb_g = normed.dot(in_out_b_g, in_out_b_g);
            C.set(n_gradients, n_gradients, valb_gb_g/a_norm2);
            C.set(n_gp1,n_gradients, 1);
            C.set(n_gradients, n_gp1, 1);
            C.set(n_gradients, n_2gp2,-1);

            for(std::size_t i = 0; i < n_gradients; ++i) {
                Scalar valii    = normed.dot(gradients[i], gradients[i]);
                Scalar vali_b_g = normed.dot(gradients[i], in_out_b_g);
                C.set(i,i,valii/a_norm2);
                C.set(i, n_gradients, vali_b_g/a_norm2);
                C.set(n_gradients, i, vali_b_g/a_norm2);
                C.set(n_gp1,i, 1);
                C.set(i,n_gp1, 1);
                C.set(i,n_gp2 + i,-1);

                for(std::size_t j = i+1; j < n_gradients; ++j) {
                    Scalar valij = normed.dot(gradients[i], gradients[j]);
                    C.set(i, j, valij/a_norm2);
                    C.set(j, i, valij/a_norm2);
                }
            }

            auto lb = std::make_shared<Vector>(zeros(m_g));

            Vector lambda = zeros(m_g),
            r_k = zeros(m_g),
            p_k = zeros(m_g),
            Cp_k = zeros(m_g),
            val_newton_func = zeros(m_g),
            desc_dir = zeros(m_g);

            // compute a useful initial lambda. This lambda represents the smallest element of conv{b_g , gradients[0] }.
            lambda.set(n_gradients-1 , std::min(1.,std:: max(0.,  (C.get(n_gradients,n_gradients)-C.get(n_gradients-1,n_gradients))/( C.get(n_gradients-1,n_gradients-1)-2*C.get(n_gradients-1,n_gradients)+C.get(n_gradients,n_gradients)))));
            //if (lambda.get(0)<0)  {lambda.set(0, 0);}// xx is wrong in the alternative!!!!!
            //if (lambda.get(0)>1)  {lambda.set(0, 1);}// xx is wrong in the alternative!!!!!
            lambda.set(n_gradients, 1-lambda.get(n_gradients-1) ); // xx is wrong in the alternative!!!!!

            //compute the value of the function F(lambda)=val_newton_func and actualize the gradient C of this function

            val_newton_func.set(n_gp1,-1);
            for (std::size_t i=0; i< n_gp1; i++){
                val_newton_func.set(i, lambda.get(n_gp1)-lambda.get(n_gp2+i) );
                val_newton_func.set(n_gp1, lambda.get(i) + val_newton_func.get(n_gp1) ); //xx better?
                for(std::size_t j=0; j< n_gp1; j++){
                    val_newton_func.set(i, C.get(i,j)*lambda.get(j) + val_newton_func.get(i) ); //xx better?
                }
                val_newton_func.set(i+n_gp2, lambda.get(i)-std:: max(0.,lambda.get(i)-lambda.get(i+n_gp2)));
                if(lambda.get(i) > lambda.get(i+n_gp2)){
                    C.set(i+n_gp2, i, 0);
                    C.set(i+n_gp2, i+n_gp2, 1);
                }
                else{
                    C.set(i+n_gp2, i, 1);
                    C.set(i+n_gp2, i+n_gp2, 0);
                }
            }

            int counter1=0;
            while( dot(val_newton_func,val_newton_func) > .000000000001 && counter1 < 500 ){ //Begin Semismooth Newton Method to solve F(lambda)=0

                //solve C desc_dir=val_newton_func by solving C^TC desc_dir= C^T val_newton_func, since C is not symmetric and might not be regular.

                int counter2=0;
                Scalar be_k,al_k,res_k;
                desc_dir.set(0.);
                r_k= transpose(C) *val_newton_func;
                p_k=r_k;
                res_k=dot(r_k,r_k);
                while (res_k>.000000000001 && counter2<500)     //Begin loop of the CG-method xx
                {
                    counter2++;
                    Cp_k=transpose(C)*C*p_k;
                    al_k=res_k/dot(p_k,Cp_k);
                    desc_dir+=al_k*p_k;
                    r_k-=al_k*Cp_k;
                    be_k=1/res_k;
                    res_k=dot(r_k,r_k);
                    be_k*=res_k;
                    p_k=r_k+be_k*p_k;
                }                                                      //Ende CG-loop. Equation is solved.

                //make a Newton step

                lambda-= desc_dir;

                //compute the value of the function F(lambda)=val_newton_func and actualize the gradient C of this function

                val_newton_func.set(n_gp1,-1);
                for (std::size_t i=0; i< n_gp1; i++){
                    val_newton_func.set(i, lambda.get(n_gp1) - lambda.get(n_gp2+i) );
                    val_newton_func.set(n_gp1, lambda.get(i) + val_newton_func.get(n_gp1) ); //xx better?
                    for(std::size_t j=0; j< n_gradients+1; j++){
                        val_newton_func.set(i, C.get(i,j)*lambda.get(j) + val_newton_func.get(i) ); //xx better?
                    }
                    val_newton_func.set(i+n_gp2, lambda.get(i)-std:: max(0.,lambda.get(i)-lambda.get(i+n_gp2)));
                    if(lambda.get(i) > lambda.get(i+n_gp2)){
                        C.set(i+n_gp2, i, 0);
                        C.set(i+n_gp2, i+n_gp2, 1);
                    }
                    else{
                        C.set(i+n_gp2, i, 1);
                        C.set(i+n_gp2, i+n_gp2, 0);
                    }
                }

                counter1++;
            }// End Semismooth Newton Method
            // b_g is no longer needed and can be used as temporary variable
            in_out_b_g *= lambda.get(n_gradients);
            for (std::size_t i = 0; i < n_gradients; i++) {
                in_out_b_g += lambda.get(i) * gradients[i];
            }

            gradients.push_back(in_out_b_g);
            //  xx                               if (n_gradients > max_gradients)
            //  xx                                   del oldest gradient
        }

    public:
    };

    /**
     * @brief
     */
    template<class Matrix, class Vector, int Backend = Traits<Vector>::Backend>
    class MSSolver final : public NewtonBase<Matrix, Vector>
    {
        typedef UTOPIA_SCALAR(Vector)    Scalar;
        typedef UTOPIA_SIZE_TYPE(Vector) SizeType;
        using LinearSolverT = LinearSolver<Matrix, Vector>;

    public:

        enum NormType {
            L2_NORM = 0,
            A_NORM = 1,
            A_SQUARED_NORM = 2,
            N_NORM_TYPES
        };

        using HilbertFunction = utopia::HilbertFunction<Matrix, Vector>;

        class L2HilbertFunction final : public HilbertFunction {
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

        class AHilbertFunction final : public HilbertFunction {
        public:
            AHilbertFunction(const std::shared_ptr<LinearSolverT> &linear_solver)
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

        class ASquaredHilbertFunction final : public HilbertFunction {
        public:
            ASquaredHilbertFunction(const std::shared_ptr<LinearSolverT> &linear_solver)
            : M_inv(linear_solver)
            {}

            void update(const std::shared_ptr<Matrix> &mat) override
            {
                assert(mat);
                M = mat;

                if(!M2) {
                    M2 = std::make_shared<Matrix>();
                }

                assert(M2);

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

        MSSolver(const std::shared_ptr<LinearSolverT> &linear_solver):
        NewtonBase<Matrix, Vector>(linear_solver),
        delta_(0.3),
        delta_prime_(0.35),
        radius0_(3.),
        T1_([this](const Scalar &radius) -> Scalar { return radius/radius0_; }),
        T2_([](const Scalar &radius) -> Scalar { return 0.35 * radius; }),
        G_([](const Scalar &/*radius*/, const Scalar &radius0) -> Scalar { return radius0; }),
        norm_type_(L2_NORM),
        convex_hull_n_gradients_(2)
        {
            normed_.resize(N_NORM_TYPES);
            normed_[L2_NORM] = std::make_shared<L2HilbertFunction>();
            normed_[A_NORM]  = std::make_shared<AHilbertFunction>(std::shared_ptr<LinearSolverT>(linear_solver->clone()));
            normed_[A_SQUARED_NORM] = std::make_shared<ASquaredHilbertFunction>(std::shared_ptr<LinearSolverT>(linear_solver->clone()));
        }

        void set_convex_hull_n_gradients(const SizeType n)
        {
            convex_hull_n_gradients_ = n;
        }

        inline void set_norm_type(const NormType norm_type)
        {
            norm_type_ = norm_type;
        }

        // bool convex_hull_minmia

        bool B(Function<Matrix, Vector> &fun,
               HilbertFunction &normed,
               const Vector &left_in,
               const Vector &right_in,
               const Vector &dir,
               const Scalar &tol,
               const Scalar &/*step_size*/, //not needed
               const Scalar &val_left_in,
               const Scalar &val_right_in,
               Vector &h_g)
        {
            Scalar val_left = val_left_in;
            Scalar val_right = val_right_in;

            leftpoint = left_in;
            midpoint  = right_in;


            // std::cout << "step_size: " << step_size << std::endl;

            //we compute the directional derivative

            fun.gradient(midpoint, g_buff);
            Scalar d = dot(g_buff, dir);

            bool success = true;
            bool first = true;
            SizeType it = 0;
            while(d > tol) {
                Scalar f_val = 0.;
                fun.value(midpoint, f_val);

                if(first) {
                    rightpoint = right_in;
                    midpoint  = 0.5 * (left_in + right_in);
                    first = false;
                } else {
                    auto left_diff = f_val - val_left;
                    auto right_diff = val_right - f_val;

                    // std::cout << left_diff << "<=" << right_diff << std::endl;

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

                fun.gradient(midpoint, g_buff);
                d = dot(g_buff, dir);

                ++it;

                if(it > 100) {
                    success = false;
                    break;
                }
            }

            normed.transform_gradient(g_buff, h_g);
            return success;
        }

        bool line_search(Function<Matrix, Vector> &fun,
                         const Vector &x,
                         const Vector &dir,
                         Scalar &alpha)
        {
            x_new = x + alpha * dir;
            Scalar f_x = 0.;
            Scalar f_x_new = 0.;

            fun.value(x, f_x);
            fun.value(x_new, f_x_new);

#ifndef NDEBUG
            const Scalar f_x_old = f_x;
#endif //NDEBUG

            SizeType n_alpha = 0;
            while(f_x_new < f_x) {
                f_x = f_x_new;
                ++n_alpha;

                x_new += alpha * dir;
                fun.value(x_new, f_x_new);
            }

            alpha *= n_alpha;

            if(!n_alpha) {
                assert(false && "maybe use std line-search");
                return false;
            }

            assert(f_x <= f_x_old);
            return true;
        }

        bool solve(Function<Matrix, Vector> &fun, Vector &x) override
        {
            using namespace utopia;

            gradients_.clear();

            Vector y, h_g, b_g;

            Scalar g_norm = 1.;
            SizeType it = 0;

            bool converged = false;

            this->init_solver("MSSolver", {" it. ", "|| g ||", "f(x)"});

            Scalar radiusk = radius0_;
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
                radiusk = G_(g_norm, radiusk);

                // // print iteration status on every iteration

                // // check convergence and print interation info
                converged = this->check_convergence(it, g_norm, 1, 1);

                if(converged) {
                    break;
                }

                Scalar val_x = 0.;
                fun.value(x, val_x);


                if(this->verbose_) {
                    PrintInfo::print_iter_status(it, { g_norm, val_x });
                }


                bool stop = false;
                while(!stop) {
                    gradients_.clear();
                    gradients_.push_back(h_g);

                    auto a_norm = g_norm;
                    auto a_norm2 = a_norm * a_norm;

                    while(a_norm > T1_(radiusk)) {

                        const auto &dir = gradients_.back();
                        Scalar step_size = (radiusk/a_norm);
                        assert(step_size > 0.);

                        y = x - step_size * dir;
                        Scalar val_y = 0.;
                        fun.value(y, val_y);
                        const Scalar diff_val = val_y - val_x;

                        if(diff_val < (-delta_ * a_norm * radiusk)) {
                            Scalar alpha = radiusk/a_norm;
                            if(!line_search(fun, x, -dir, alpha)) {
                                //use std line-search
                            }

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

                                assert(false && "handle case?!");
                                return false;
                            }

                            if(!convex_hull_solver_ || convex_hull_n_gradients_ == 2) {
                                //case with only 2 vectors
                                //compute smallest element of the set of gradients
                                const auto &first_g = gradients_.back();
                                const auto n_bg = normed->norm(b_g);
                                const auto n_bg2 = n_bg * n_bg;
                                const auto dot_bg_fg = normed->dot(b_g, first_g);


                                //Seems harmful for convergence
                                // const auto lambda = std::min(1., std::max(0., (n_bg2 - dot_bg_fg)/(n_bg2 + a_norm2 - 2. * dot_bg_fg)));

                                //Seems harmful for B(x)
                                const auto lambda = (n_bg2 - dot_bg_fg)/(n_bg2 + a_norm2 - 2. * dot_bg_fg);

                                //add small elements to list
                                gradients_.back() = (lambda * first_g + (1.-lambda) * b_g);
                            } else {
                                convex_hull_solver_->solve(*normed, a_norm2, gradients_, b_g);

                                if(static_cast<SizeType>(gradients_.size()) >= convex_hull_n_gradients_) {
                                    gradients_.erase(gradients_.begin());
                                }
                            }

                            a_norm = normed->norm(gradients_.back());
                            a_norm2 = a_norm * a_norm;

                            if(radiusk < this->atol() && a_norm < this->atol()) {
                                stop = true;
                                converged = true;
                                break;
                            }
                        }
                    }

                    if(!stop) {
                        radiusk = T2_(radiusk);
                    }
                }


                it++;
            }

            this->print_statistics(it);
            return true;
        }

        void read(Input &in) override
        {
            NonLinearSolver<Vector>::read(in);
            // in.get("dumping", delta_);
        }

        void set_convex_hull_solver(const std::shared_ptr<IMSConvexHullSolver<Matrix, Vector>> &solver)
        {
            convex_hull_solver_ = solver;
        }

    private:
        Scalar delta_, delta_prime_, radius0_;   /*!< Dumping parameter. */
        std::function<Scalar(const Scalar &)> T1_, T2_;
        std::function<Scalar(const Scalar &, const Scalar &)> G_;
        // std::shared_ptr<LSStrategy> ls_strategy_;     /*!< Strategy used in order to obtain step \f$ \delta_k \f$ */
        NormType norm_type_;

        std::vector<std::shared_ptr<HilbertFunction>> normed_;
        std::vector<Vector> gradients_;
        SizeType convex_hull_n_gradients_;

        Vector x_new;
        Vector leftpoint;
        Vector midpoint;
        Vector rightpoint;
        Vector g_buff;

        std::shared_ptr<IMSConvexHullSolver<Matrix, Vector>> convex_hull_solver_;
    };

}
#endif //UTOPIA_MS_SOLVER_HPP
