#ifndef UTOPIA_GENERALIZED_CAUCHY_POINT_HPP
#define UTOPIA_GENERALIZED_CAUCHY_POINT_HPP

#include "utopia_TRSubproblem.hpp"
#include "utopia_CauchyPoint.hpp"
#include "utopia_Parameters.hpp"
#include "utopia_LinearSolverInterfaces.hpp"

namespace utopia
{

	template<class Matrix, class Vector>
    class GeneralizedCauchyPoint 
    {
        typedef UTOPIA_SCALAR(Vector) Scalar;
        typedef UTOPIA_SIZE_TYPE(Vector) SizeType;
        
        public:

            GeneralizedCauchyPoint(const Parameters & params = Parameters()): cp_memory_(5)
            { }

            // inline GeneralizedCauchyPoint * clone() const override
            // {
            //     return new GeneralizedCauchyPoint(*this);
            // }

            void set_memory_size(const SizeType & m)
            {
                cp_memory_ = m;
            }

            SizeType get_memory_size() const 
            {
                return cp_memory_;
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



    void get_breakpoints(const Vector & d, const Vector & x, const Vector & lb, const Vector & ub, Vector &break_points, const Scalar & delta = 9e9) const
    {
        if(empty(break_points) || local_size(break_points)!=local_size(x))
            break_points = local_values(local_size(x).get(0), 0);

        {
          Read<Vector> r_ub(ub), r_lb(lb), r_x(x), r_d(d);
          Write<Vector> wt(break_points); 

          each_write(break_points, [ub, lb, x, d, delta](const SizeType i) -> double { 
                Scalar li =  lb.get(i); Scalar ui =  ub.get(i); Scalar xi =  x.get(i);  Scalar di =  d.get(i);  
                if(di > 0)
                {
                    Scalar val = ui - xi; 
                    return (std::min(val, delta))/di; 
                }
                else if(di < 0)
                {
                    Scalar val = li - xi; 
                    return (std::max(val, -1.0 * delta))/di; 
                }
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

            each_write(active_set, [break_points, t_break](const SizeType i) -> double 
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
            SizeType cp_memory_; // memory size


    };
}

#endif //UTOPIA_GENERALIZED_CAUCHY_POINT_HPP