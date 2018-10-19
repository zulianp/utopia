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


        void compute_breakpoints(const Vector & g, const Vector & x, const Vector & lb, const Vector & ub, Vector &t) const
        {
            auto inf = std::numeric_limits<Scalar>::infinity(); 

            if(empty(t) || local_size(t)!=local_size(x))
                t = local_values(local_size(x).get(0), inf);

            // TODO:: add checks if there are both bounds available 
            {
              Read<Vector> r_ub(ub), r_lb(lb), r_x(x), r_d(g);
              Write<Vector> wt(t); 

              each_write(t, [ub, lb, x, g, inf](const SizeType i) -> double { 
                          Scalar li =  lb.get(i); Scalar ui =  ub.get(i); Scalar xi =  x.get(i);  Scalar gi =  g.get(i);  
                          if(gi < 0)
                            return (xi - ui)/gi; 
                          else if(gi > 0)
                            return (xi - li)/gi; 
                        else 
                            return inf; 
            }  );
          }
        }

        void get_d_corresponding_to_ti(const Vector & t, const Vector & g, Vector &d, const Scalar & t_current) const
        {
            d = -1.0 * g; 

            {   // begin lock
                Write<Vector>  wd(d);
                Read<Vector>   rt(t);

                Range rr = range(d);

                for (SizeType i = rr.begin(); i != rr.end(); ++i)
                {
                    Scalar ti = t.get(i); 
                    if(std::abs(ti - t_current) < 1e-12)
                        d.set(i, 0.0); 
                }

            } // end of lock
        }

        void get_initial_feasible_set(const Vector & break_points, Vector & feasible_set) const
        {
            feasible_set = local_values(local_size(break_points).get(0), 0.0); 

            {   // begin lock
                Read<Vector>  wd(break_points);
                Write<Vector> r(feasible_set); 

                Range rr = range(break_points);

                for (SizeType i = rr.begin(); i != rr.end(); ++i)
                {
                    if(break_points.get(i) > 0)
                        feasible_set.set(i, 1); 
                }
       
            } // end of lock
        }



        Scalar get_next_t(const Vector & sorted_break_points, const SizeType & index) const
        {
            Vector t_help = local_values(1, 0.0); 
            Scalar value=0.0; 

            // this is horrible solution, but lets fix it later 
            {
                Read<Vector> r(sorted_break_points); 

                auto rr = range(sorted_break_points);
                for (SizeType i = rr.begin(); i != rr.end(); ++i)
                {
                    if(i==index)
                        value = sorted_break_points.get(i); 
                }
            }

            // this is horrible solution, but lets fix it later 
            {
                Write<Vector> w(t_help); 

                auto rr = range(t_help);
                for (SizeType i = rr.begin(); i != rr.end(); ++i)
                        t_help.set(i, value); 
            }        

            return sum(t_help); 
        }


        Scalar project_direction_on_boundary(const Vector & x, const Vector & d, const Vector & ub, const Vector & lb, Vector & x_cp, const SizeType & active_index) const
        {
            Vector x_local = local_values(1, 0.0); 
            Scalar val=0; 

            {   // begin lock
                Write<Vector>  wd(x_cp);

                Read<Vector> r1(ub); 
                Read<Vector> r2(lb); 
                Read<Vector> r3(d); 
                Read<Vector> r4(x); 

                Range rr = range(x_cp);

                for (SizeType i = rr.begin(); i != rr.end(); ++i)
                {
                    if(i==active_index)
                    {
                        val = (d.get(i)>0)? ub.get(i) : lb.get(i); 
                        x_cp.set(i, val); 
                    }
                }
       
            } // end of lock


            ///////////////////////////////////////////////////////////
            {   // begin lock
                Write<Vector>  wx(x_local);
                Range rr_x_local = range(x_local);

                x_local.set(rr_x_local.begin(), val); 
       
            } // end of lock

            return sum(x_local); 
        }


        bool get_global_active_index(const Vector & break_points, Vector & feasible_set, const Scalar & t_current, SizeType & index) const
        {
            Scalar counter=0.0; 

            // TODO:: put inf
            Vector indices = local_values(1, 9e9); 
            Vector counter_vec = local_values(1, 0); 

            {   // begin lock
                Read<Vector>  wd(break_points);
                Read<Vector>  rv(feasible_set); 
                Write<Vector>  rv2(indices); 

                Range rr = range(break_points);
                Range ind_range = range(indices);

                for (SizeType i = rr.begin(); i != rr.end(); ++i)
                {
                    Scalar value = break_points.get(i); 
                    if(std::abs(value - t_current) < 1e-12 && feasible_set.get(i)==1)
                    {
                        indices.set(ind_range.begin(), i); 
                        break; 
                    }
                }
            }

            index = min(indices); 

            if(index<0 || index > size(break_points).get(0))
                utopia_error("L-BFGS-B::get_global_active_index: index not valid. "); 


            {
                Write<Vector>  rv(feasible_set); 
                Range rr = range(feasible_set);

                for (SizeType i = rr.begin(); i != rr.end(); ++i)
                {
                    if(i==index)
                        feasible_set.set(i, 0); 
                }            

            }

            {
                Read<Vector>  rv(feasible_set); 
                Range rr = range(feasible_set);

                // horrible solution - looops should be merged 
                for (SizeType i = rr.begin(); i != rr.end(); ++i)
                {
                    Scalar value = break_points.get(i); 
                    if(std::abs(value - t_current) < 1e-12 && feasible_set.get(i)==1)
                    {
                        counter++; 
                    }
                }
            }    

            {
                Write<Vector> wvv(counter_vec); 
                Range r_counter = range(counter_vec);

                for (SizeType i = r_counter.begin(); i != r_counter.end(); ++i)
                    counter_vec.set(i, counter); 

            } // end of lock

            SizeType counter_global_sum = sum(counter_vec); 
            return (counter_global_sum>1) ? true: false; 
        }


    SizeType get_number_of_sorted_break_points(const Vector & sorted_break_points) const
    {
        Vector help = local_values(1, 0.0); 
        SizeType val = size(sorted_break_points).get(0); 

        // this is horrible solution, but lets fix it later 
        {
            Write<Vector> w(help); 

            if(mpi_world_rank()==0)
                help.set(0, val); 
        }

        return sum(help); 
    }



    // use approxeq for all stuff where u compare 
    void add_d_to_x(const Vector & x, Vector & x_cp, const Vector & feasible_set,  const Vector & d, const Scalar & tau) const
    {

        {   // begin lock
            Write<Vector>  wd(x_cp);
            Read<Vector>  r1(x);
            Read<Vector>  r2(feasible_set);
            Read<Vector>  r3(d);

            Range rr = range(x_cp);

            for (SizeType i = rr.begin(); i != rr.end(); ++i)
            {
                if(feasible_set.get(i) == 1.0)
                {
                    x_cp.set(i, x.get(i) + tau * d.get(i)); 
                }
                    
            }
   
        } // end of lock

    }


    void zero_dir_component(Vector & d, const SizeType & index) const
    {
        {   // begin lock
            Write<Vector>  wd(d);
            Range rr = range(d);

            for (SizeType i = rr.begin(); i != rr.end(); ++i)
            {
                if(i==index)
                    d.set(i, 0); 
            }
   
        } // end of lock

    }

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////// CHECKed //////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


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