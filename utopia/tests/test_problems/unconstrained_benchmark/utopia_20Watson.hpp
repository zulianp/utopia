#ifndef UTOPIA_SOLVER_WATSON_20
#define UTOPIA_SOLVER_WATSON_20

#include <vector>
#include <assert.h>
#include "utopia_Function.hpp"


namespace utopia
{
    template<class Matrix, class Vector>
    class Watson20 final: public UnconstrainedTestFunction<Matrix, Vector>
    {
    public:
        DEF_UTOPIA_SCALAR(Matrix)
        typedef UTOPIA_SIZE_TYPE(Vector) SizeType;

        Watson20()
        {
            assert(!utopia::is_parallel<Matrix>::value || mpi_world_size() == 1 && "does not work for parallel matrices");

            x_init_ = zeros(9);
            x_exact_ = zeros(9);

            {
                const Write<Vector> write2(x_exact_);
                const Write<Vector> write1(x_init_);

                x_exact_.set(0, -0.000015);
                x_exact_.set(1, 0.999790);
                x_exact_.set(2, 0.014764);
                x_exact_.set(3, 0.146342);
                x_exact_.set(4, 1.000821);
                x_exact_.set(5, -2.617731);
                x_exact_.set(6, 4.104403);
                x_exact_.set(7, -3.143612);
                x_exact_.set(8, 1.052627);
            }
        }

        std::string name() const override
        {
            return "Watson";
        }

        SizeType dim() const override
        {
            return 9.0;
        }

        bool exact_sol_known() const override
        {
            return false;  // just because we can not fit into precision
        }


        bool value(const Vector &x, typename Vector::Scalar &result) const override
        {
            if( mpi_world_size() > 1){
                utopia_error("Function is not supported in parallel... \n");
                return false;
            }

            assert(x.size() == this->dim());
            {
                const Read<Vector> read(x);

                result = 0.0;
                for(SizeType i = 1; i <=29; i++)
                {
                    Scalar s1 = 0.0;
                    Scalar d = 1.0;
                    for(SizeType j = 2; j <= this->dim(); j++)
                    {
                      s1 = s1 + ((j - 1.) * d * x.get(j-1));
                      d = d * i / 29.0;
                    }

                    Scalar s2 = 0.0;
                    d = 1.0;
                    for(SizeType j = 1; j <= this->dim(); j++)
                    {
                        s2 = s2 + (d * x.get(j-1));
                        d = d*i/29.0;
                    }

                    Scalar help = ( s1 - (s2 * s2) - 1.0 );
                    result += help*help;
                }

                Scalar help2 = x.get(1) - (x.get(0)*x.get(0)) - 1.0;
                result += x.get(0)*x.get(0) + help2*help2;
            }

            return true;
        }

        bool gradient(const Vector &x, Vector &g) const override
        {
            if( mpi_world_size() > 1){
                utopia_error("Function is not supported in parallel... \n");
                return false;
            }

            assert(x.size() == this->dim());
            g = zeros(this->dim());

            {
                const Read<Vector> read(x);
                const Write<Vector> w(g);

                std::vector<Scalar> g_help(this->dim());

                for(SizeType i = 1; i <=29; i++)
                {
                    Scalar s1 = 0.0;
                    Scalar d1 = i / 29.0;
                    Scalar d2 = 1.0;
                    for(SizeType j = 2; j <= this->dim(); j++)
                    {
                        s1 += (( j - 1. ) * d2 * x.get(j-1));
                        d2 = d2 * i / 29.0;
                    }

                    Scalar s2 = 0.0;
                    d2 = 1.0;
                    for(SizeType j = 1; j <= this->dim(); j++)
                    {
                        s2 += (d2 * x.get(j-1));
                        d2 = d2 * i / 29.0;
                    }

                    Scalar t = s1 - (s2 * s2) - 1.0;
                    Scalar s3 = 2.0 * s2 * i / 29.0; ;
                    d2 = 2.0 / d1;

                    for(SizeType j = 1; j <= this->dim(); j++)
                    {
                        g_help[j-1] += d2*((j-1.)-s3)*t;
                        d2 = d2 *i / 29.0;
                    }
                }

                Scalar t1 = x.get(1) - (x.get(0) *x.get(0)) - 1.0;

                g_help[0] += (2.0 * x.get(0))  - (4.0 * x.get(0) * t1);
                g_help[1] += 2.0 * t1;

                for(auto i=0; i < this->dim(); i++)
                {
                    g.set(i, g_help[i]);
                }
            }


            return true;
        }

        bool hessian(const Vector &x, Matrix &H) const override
        {
            if( mpi_world_size() > 1){
                utopia_error("Function is not supported in parallel... \n");
                return false;
            }

            assert(x.size() == this->dim());
            H = zeros(this->dim(), this->dim());

            std::vector<std::vector<Scalar> > hess(this->dim(), std::vector<Scalar>(this->dim()));

            {
                const Read<Vector> read(x);
                const Write<Matrix> write(H);


                for(SizeType i = 1; i <=29; i++)
                {
                    Scalar d1 = i / 29.0;
                    Scalar d2 = 1.0;
                    Scalar s1 = 0.0;
                    Scalar s2 = x.get(0);

                    for(SizeType j = 2; j <= this->dim(); j++)
                    {
                        s1 += (j-1.)*d2*x.get(j-1);
                        d2 *= d1;
                        s2 += d2 * x.get(j-1);
                    }

                    Scalar t = 2.0 * ( s1 - (s2 * s2) - 1.0 ) * d1 * d1;
                    Scalar s3 = 2.0 * d1 * s2;
                    d2 = 1.0 / d1;

                    for(SizeType j = 1; j <= this->dim(); j++)
                    {
                        Scalar t1 = j - 1. - s3;
                        hess[j-1][j-1] += 2.0 * ((t1*t1) - t)*d2*d2;
                        Scalar d3 = 1.0 / d1;

                        for(SizeType k = 1; k < j; k++)
                        {
                            hess[j-1][k-1] += 2.0*(t1*((k-1.)-s3)-t)*d2*d3;
                            d3 *= d1;
                        }
                        d2 *= d1;
                    }
                }

                Scalar t3 = x.get(1) - (x.get(0)*x.get(0)) - 1.0;
                hess[0][0] += 2.0 - (4.0 * ( t3 - (2.0*x.get(0)*x.get(0))));
                hess[1][1] += 2.0;
                hess[1][0] -= 4.0*x.get(0);


                for(auto i=0; i < this->dim(); i++)
                {
                    for(auto j=0; j < this->dim(); j++)
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
            return 1.39976e-06;
        }

    private:
        Vector x_init_;
        Vector x_exact_;

    };

}

#endif //UTOPIA_SOLVER_WATSON_20
