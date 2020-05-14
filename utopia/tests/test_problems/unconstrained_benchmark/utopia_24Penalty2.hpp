#ifndef UTOPIA_SOLVER_PENALTY_II
#define UTOPIA_SOLVER_PENALTY_II

#include <cassert>
#include <vector>
#include "utopia_Function.hpp"

namespace utopia
{
    template<class Matrix, class Vector>
    class PenaltyII24 final: public UnconstrainedTestFunction<Matrix, Vector>
    {
    public:
        using Traits   = utopia::Traits<Vector>;
        using Scalar   = typename Traits::Scalar;
        using SizeType = typename Traits::SizeType;
        using Comm     = typename Traits::Communicator;

        PenaltyII24()
        {
            auto v_layout = serial_layout(dim());
            x_exact_.zeros(v_layout); // not known
            x_init_.values(v_layout, 0.5);
        }

        std::string name() const override
        {
            return "Penalty II";
        }

        SizeType dim() const override
        {
            return 10;
        }

        bool exact_sol_known() const override
        {
            return false;
        }


        bool value(const Vector &x, Scalar &result) const override
        {
           if( x.comm().size() > 1){
                utopia_error("Function is not supported in parallel... \n");
                return false;
            }

            const SizeType n = this->dim();
            assert(size(x).get(0) == n);

            Scalar alpha = 0.00001;

            Scalar t1 = -1.0;
            Scalar t2 = 0.0;
            Scalar t3 = 0.0;
            Scalar d2 = 1.0;
            Scalar s2 = 0.0;
            Scalar b = 0.0;

            {
                Read<Vector>read(x);

                for(auto j=1; j <=n; j++)
                {
                    t1 = t1 + (( n - j + 1 ) * x.get(j-1)*x.get(j-1));
                    Scalar s1 = std::exp(x.get(j-1)/10.0);

                    if (1<j)
                    {
                        Scalar s3 = s1 + s2 - (d2 * (std::exp(0.1)+1.0));
                        t2 += (s3 * s3);
                        Scalar a = ( s1 - (1.0 /std::exp(0.1)));
                        t3 +=  a*a;
                    }

                    s2 = s1;
                    d2 *= std::exp(0.1);
                }

                b = x.get(0) - 0.2;
            }


            result = alpha * ( t2 + t3 ) + (t1 * t1) + (b*b);

            return true;
        }

        bool gradient(const Vector &x, Vector &g) const override
        {
           if( x.comm().size() > 1){
                utopia_error("Function is not supported in parallel... \n");
                return false;
            }

            const SizeType n = this->dim();
            assert(size(x).get(0) == n);

            if(empty(g))
            {
                g.zeros(layout(x));
            }


            Scalar alpha = 0.00001;
            Scalar x0 = 0;

            Scalar t1 = -1.0;
            {
                Read<Vector>read(x);
                for(auto j=1; j <=n; j++)
                {
                    t1 = t1 + ( n - j + 1 ) * (x.get(j-1)*x.get(j-1));
                }
                x0 = x.get(0);
            }

            std::vector<Scalar> grad(n);

            {
                Read<Vector>read(x);

                Scalar d2 = 1.0;
                Scalar th = 4.0 * t1;
                Scalar s2 = 0.0;

                for(auto j=1; j <=n; j++)
                {
                    grad[j-1] = ( n - j + 1 ) * x.get(j-1) * th;
                    Scalar s1 = std::exp(x.get(j-1)/10.0);

                    if(1<j)
                    {
                      Scalar s3 = s1 + s2 - (d2 *(std::exp(0.1)+1.0));
                      grad[j-1] += alpha*s1*(s3 + s1 - (1.0 /std::exp(0.1)))/5.0;
                      grad[j-2] += alpha*s2*s3/5.0;
                    }

                    s2 = s1;
                    d2 = d2 * std::exp(0.1);
                }
            }


            {
                Write<Vector>w(g);
                auto r = range(g);

                for(auto i = r.begin(); i != r.end(); ++i)
                {
                    Scalar val;
                    if(i==0)
                    {
                        val  = grad[i] + (2.0 * (x0 - 0.2));
                    }
                    else
                    {
                        val = grad[i];
                    }

                    g.set(i, val);
                }
            }



            return true;
        }

        bool hessian(const Vector &x, Matrix &H) const override
        {
           if( x.comm().size() > 1){
                utopia_error("Function is not supported in parallel... \n");
                return false;
            }

            const SizeType n = this->dim();
            assert(local_size(x).get(0) == n);

            Scalar alpha = 0.00001;
            Scalar t1 = - 1.0;

            H.dense(square_matrix_layout(layout(x)), 0.0);

            {
                Read<Vector>read1(x);

                for(auto j=1; j <=n; j++)
                {
                    t1 += (n-j+1)* x.get(j-1)* x.get(j-1);
                }
            }

            Scalar d1 = std::exp(0.1);
            Scalar d2 = 1.0;
            Scalar s2 = 0.0;
            Scalar th = 4.0 * t1;

            std::vector<std::vector<Scalar> > hess(n, std::vector<Scalar>(n));

            {
                Read<Vector>read1(x);

                for(auto j=1; j <=n; j++)
                {
                    Scalar factor = ( n - j + 1 );
                    Scalar a = factor * x.get(j-1);
                    hess[j-1][j-1] = (8.0 * a*a) + (factor * th);
                    Scalar s1 = std::exp(x.get(j-1)/10.0);

                    if ( 1 < j )
                    {

                        Scalar s3 = s1 + s2 - (d2 *(d1 + 1.0 ));
                        hess[j-1][j-1] +=  alpha*s1*( s3 + s1 - 1.0/d1 + (2.0 * s1))/50.0;
                        hess[j-2][j-2] += alpha*s2*(s2+s3)/50.0;

                        for(auto k=1; k <=j-1; k++)
                        {
                            t1 = std::exp(k/10.0);
                            hess[j-1][k-1] = 8.0 * factor*( n - k + 1 ) * x.get(j-1) * x.get(k-1);
                        }

                        hess[j-1][j-2] += alpha * s1 * s2 / 50.0;
                    }

                s2 = s1;
                d2 = d1 * d2;

                }
            }

            hess[0][0] += 2.0;


            {
                Write<Matrix> wr1(H);

                for(auto i=0; i < n; i++)
                {
                    for(auto j=0; j < n; j++)
                    {
                        H.set(i, j,  hess[i][j]);
                    }
                }

            }

            // this could be done way much nicer...
            Matrix D = diag(diag(H));
            H = H + transpose(H) - D;


            return true;
        }

        Vector initial_guess() const override
        {
            return x_init_;
        }

        const Vector & exact_sol() const override
        {
            return x_exact_;
        }

        Scalar min_function_value() const override
        {
            return 2.93660e-4;
        }



    private:
        Vector x_init_;
        Vector x_exact_;

    };

}

#endif //UTOPIA_SOLVER_PENALTY_II
