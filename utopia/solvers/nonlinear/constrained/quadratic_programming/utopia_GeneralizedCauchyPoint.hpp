#ifndef TR_GENERALIZED_CAUCHY_POINT_SUBPROBLEM
#define TR_GENERALIZED_CAUCHY_POINT_SUBPROBLEM

#include <string>
#include "utopia_BoxConstraints.hpp"
#include "utopia_QPSolver.hpp"

namespace  utopia 
{

    template<class Matrix, class Vector>
    class GeneralizedCauchyPoint final: public MatrixFreeQPSolver<Vector>, public QPSolver<Matrix, Vector>
    {
        typedef UTOPIA_SCALAR(Vector) Scalar;


        public:
            GeneralizedCauchyPoint(): cp_memory_(5)
            {
                
            }
            
            ~GeneralizedCauchyPoint( ){}

            GeneralizedCauchyPoint * clone() const override
            {
                return new GeneralizedCauchyPoint(*this);
            }


            void set_memory_size(const SizeType & m)
            {
                cp_memory_ = m;
            }

            SizeType get_memory_size() const 
            {
                return cp_memory_;
            }


            bool solve(const Operator<Vector> &A, const Vector &rhs, Vector &sol) override
            {
                auto &box = this->get_box_constraints(); 
                return aux_solve(A, -1.0 *rhs, sol, box); 
            }


            bool solve(const Matrix &A, const Vector &rhs, Vector &sol) override
            {
                auto A_op_ptr = utopia::op_ref(A);
                auto &box = this->get_box_constraints(); 
                return aux_solve(*A_op_ptr, -1.0 *rhs, sol, box); 
            }


        private:


            bool aux_solve(const Operator<Vector> &H,  const Vector &g, Vector &s, const BoxConstraints<Vector> & constraints)
            {
                Scalar f_p, f_pp, t_current, t_next, dt, gd, delta_diff;
                Vector break_points, sorted_break_points, active_set, e, Hd; 

                const auto &ub = constraints.upper_bound();
                const auto &lb = constraints.lower_bound();

                bool converged = false; 
                SizeType num_uniq_break_points, it=0; 

                Vector d = -1.0 * g; 
                s = 0 * d; 

                this->get_breakpoints(d, *lb, *ub, break_points); 
                vec_unique_sort_serial(break_points, sorted_break_points, this->get_memory_size()); 
                num_uniq_break_points = this->get_number_of_sorted_break_points(sorted_break_points); 
                t_current = 0.0; 
                this->get_breakpoint_active_set(break_points, t_current, active_set); 
                e = e_mul(active_set, d); 
                d = d - e; 
                gd = dot(g, d); 
                H.apply(d, Hd);

                while(it < num_uniq_break_points && !converged)
                {

                    f_p = gd + dot(s, Hd); 
                    f_pp = dot(d, Hd); 

                    t_next = (it==num_uniq_break_points)? 9e9 : this->get_next_break_point(sorted_break_points, it); 

                    if(f_pp ==0 || !std::isfinite(f_pp))
                        return true; 


                    dt = - f_p/f_pp; 
                    delta_diff = t_next - t_current; 

                    if(f_p >=0)
                        converged = true; 
                    else if(f_pp >0 && dt < delta_diff)
                    {
                        s += dt * d;
                        converged = true; 
                    }

                    if(converged ==true)
                        return true; 

                    t_current = t_next; 
                    this->get_breakpoint_active_set(break_points, t_current, active_set); 
                    e = e_mul(active_set, d); 
               
                    s = s + delta_diff * d; 
                    d = d - e; 

                    gd = gd - dot(g, e); 
                    
                    Vector help; 
                    H.apply(e, help); 

                    Hd = Hd - help; 
                    it++; 
                }

                return true; 
            }


            SizeType get_number_of_sorted_break_points(const Vector & sorted_break_points) const
            {
                Vector help = local_values(1, 0.0); 
                SizeType val = size(sorted_break_points).get(0); 

                {
                    Write<Vector> w(help); 

                    if(mpi_world_rank()==0)
                        help.set(0, val); 
                }

                return sum(help); 
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

            Scalar get_next_break_point(const Vector & sorted_break_points, const SizeType & index) const
            {
                Vector t_help = local_values(1, 0.0); 
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
                    Write<Vector> w(t_help); 

                    auto rr = range(t_help);
                    for (SizeType i = rr.begin(); i != rr.end(); ++i)
                            t_help.set(i, value); 
                }        

                return sum(t_help); 
            }
             

        private:  
            SizeType cp_memory_;    // memory size
        
    };
}

#endif //TR_GENERALIZED_CAUCHY_POINT_SUBPROBLEM
