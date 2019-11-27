#ifndef TR_GENERALIZED_CAUCHY_POINT_SUBPROBLEM
#define TR_GENERALIZED_CAUCHY_POINT_SUBPROBLEM

#include <string>
#include "utopia_BoxConstraints.hpp"
#include "utopia_QPSolver.hpp"

namespace  utopia
{

    template<class Matrix, class Vector>
    class GeneralizedCauchyPoint final: public OperatorBasedQPSolver<Matrix, Vector>
    {
        typedef UTOPIA_SCALAR(Vector) Scalar;


        public:
            GeneralizedCauchyPoint(): cp_memory_(5), initialized_(false),  loc_size_(0)
            {

            }

            ~GeneralizedCauchyPoint( ){}

            GeneralizedCauchyPoint * clone() const override
            {
                return new GeneralizedCauchyPoint(*this);
            }

            void read(Input &in) override
            {
                OperatorBasedQPSolver<Matrix, Vector>::read(in);
                in.get("memory_size", cp_memory_);
            }


            void print_usage(std::ostream &os) const override
            {
                OperatorBasedQPSolver<Matrix, Vector>::print_usage(os);
                this->print_param_usage(os, "memory_size", "int", "Memory (in terms of breakpoints) used during computation.", "5.0");
            }


            void memory_size(const SizeType & m)
            {
                cp_memory_ = m;
            }

            SizeType memory_size() const
            {
                return cp_memory_;
            }

            void update(const Operator<Vector> &A) override
            {
                SizeType loc_size_rhs = A.local_size().get(0);
                if(!initialized_ || !A.comm().conjunction(loc_size_ == loc_size_rhs)) {
                    init(loc_size_rhs);
                }
            }                


            bool solve(const Operator<Vector> &A, const Vector &rhs, Vector &sol) override
            {
                auto &box = this->get_box_constraints();
                update(A);
                rhs_minus_ = -1.0 *rhs; 
                return aux_solve(A, rhs_minus_, sol, box);
            }


            bool solve(const Matrix &A, const Vector &rhs, Vector &sol) override
            {
                auto &box = this->get_box_constraints();
                update(A);
                rhs_minus_ = -1.0 *rhs; 
                return aux_solve(A, rhs_minus_, sol, box);
            }


        private:
            bool aux_solve(const Operator<Vector> &H,  const Vector &g, Vector &s, const BoxConstraints<Vector> & constraints)
            {
                Scalar f_p, f_pp, t_current, t_next, dt, gd, delta_diff;
                // Vector break_points, sorted_break_points, active_set, e, Hd;

                const auto &ub = constraints.upper_bound();
                const auto &lb = constraints.lower_bound();


                bool converged = false;
                SizeType num_uniq_break_points, it=0;

                d_ = -1.0 * g;
                s = 0.0 * d_;

                this->get_breakpoints(d_, *lb, *ub, break_points_);
                vec_unique_sort_serial(break_points_, sorted_break_points_, this->memory_size());
                num_uniq_break_points = this->get_number_of_sorted_break_points(sorted_break_points_);
                t_current = 0.0;
                this->get_breakpoint_active_set(break_points_, t_current, active_set_);
                e_ = e_mul(active_set_, d_);
                d_ = d_ - e_;
                gd = dot(g, d_);
                H.apply(d_, Hd_);

                while(it < num_uniq_break_points && !converged)
                {

                    f_p = gd + dot(s, Hd_);
                    f_pp = dot(d_, Hd_);

                    t_next = (it==num_uniq_break_points)? 9e9 : this->get_next_break_point(sorted_break_points_, it);

                    if(f_pp ==0 || !std::isfinite(f_pp))
                        return true;


                    dt = - f_p/f_pp;
                    delta_diff = t_next - t_current;

                    if(f_p >=0)
                        converged = true;
                    else if(f_pp >0 && dt < delta_diff)
                    {
                        s += dt * d_;
                        converged = true;
                    }

                    if(converged ==true)
                        return true;

                    t_current = t_next;
                    this->get_breakpoint_active_set(break_points_, t_current, active_set_);
                    e_ = e_mul(active_set_, d_);

                    s = s + delta_diff * d_;
                    d_ = d_ - e_;

                    gd = gd - dot(g, e_);

                    Vector help;
                    H.apply(e_, help);

                    Hd_ = Hd_ - help;
                    it++;
                }

                return true;
            }


            SizeType get_number_of_sorted_break_points(const Vector & sorted_break_points)
            {
                // Vector help = local_values(1, 0.0);
                t_help_.set(0.0); 
                SizeType val = size(sorted_break_points).get(0);

                {
                    Write<Vector> w(t_help_);

                    if(mpi_world_rank()==0){
                        t_help_.set(0, val);
                    }
                }

                return sum(t_help_);
            }



            void get_breakpoints(const Vector & d, const Vector & lb, const Vector & ub, Vector &break_points) const
            {

                if(empty(break_points) || local_size(break_points)!=local_size(d))
                    break_points = local_values(local_size(d).get(0), 0);

                {
                  Read<Vector> r_ub(ub), r_lb(lb), r_d(d);
                  Write<Vector> wt(break_points);

                  each_write(break_points, [&ub, &lb, &d](const SizeType i) -> double {
                        Scalar li =  lb.get(i); Scalar ui =  ub.get(i); Scalar di =  d.get(i);
                        if(di > 0.0)
                            return ui/di;
                        else if(di < 0)
                            return li/di;
                        else
                            return 0.0;
                        }  );
                }
            }


            void get_breakpoint_active_set(const Vector & break_points, const Scalar & t_break, Vector & active_set) const
            {

                if(empty(active_set) || local_size(active_set) != local_size(break_points))
                    active_set = local_values(local_size(break_points).get(0), 0.0);

                {
                    Read<Vector> ab(break_points);

                    each_write(active_set, [&break_points, t_break](const SizeType i) -> double
                    {
                        Scalar t =  break_points.get(i);
                        return (approxeq(t, t_break)) ? 1.0 : 0.0;

                    }  );

                }
            }

            Scalar get_next_break_point(const Vector & sorted_break_points, const SizeType & index)
            {
                // Vector t_help = local_values(1, 0.0);
                t_help_.set(0.0); 
                Scalar value = 0.0;

                {
                    Read<Vector> r(sorted_break_points);

                    auto rr = range(sorted_break_points);
                    for (SizeType i = rr.begin(); i != rr.end(); ++i)
                    {
                        if(i==index)
                            value = sorted_break_points.get(i);
                    }
                }

                {
                    Write<Vector> w(t_help_);

                    auto rr = range(t_help_);
                    for (SizeType i = rr.begin(); i != rr.end(); ++i)
                            t_help_.set(i, value);
                }

                return sum(t_help_);
            }


        private:


            void init(const SizeType & ls)
            {
                t_help_ = local_values(1, 0.0);

                auto zero_expr          = local_zeros(ls);
                break_points_           = zero_expr; 
                sorted_break_points_    = zero_expr; 
                active_set_             = zero_expr; 
                e_                      = zero_expr; 
                Hd_                     = zero_expr; 
                d_                      = zero_expr; 


                initialized_ = true;    
                loc_size_ = ls;      
            }


            SizeType cp_memory_;    // memory size
            Vector t_help_, break_points_, sorted_break_points_, active_set_, e_, Hd_, d_, rhs_minus_; 

            bool initialized_; 
            SizeType loc_size_;                 

    };
}

#endif //TR_GENERALIZED_CAUCHY_POINT_SUBPROBLEM
