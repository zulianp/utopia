#ifndef UTOPIA_SOLVER_EXTENDED_POWELL_SINGULAR_22
#define UTOPIA_SOLVER_EXTENDED_POWELL_SINGULAR_22

#include <vector>
#include <assert.h>
#include "utopia_Function.hpp"


namespace utopia 
{
    template<class Matrix, class Vector>
    class ExtendedPowell22 final: public UnconstrainedTestFunction<Matrix, Vector> 
    {
    public:
        DEF_UTOPIA_SCALAR(Matrix)
        typedef UTOPIA_SIZE_TYPE(Vector) SizeType;

        ExtendedPowell22(const SizeType & n): n_(n)
        {

            // n needs to be multiply of 4 
            if(n_%4 != 0)
            {
                utopia_error("ExtendedPowell22:: problem size must be multiply of 4 \n"); 
            }

            
            assert(!utopia::is_parallel<Matrix>::value || mpi_world_size() == 1 && "does not work for parallel matrices");
            
            x_exact_ = values(n_, 0.0); 
            x_init_ = values(n_, 0.0);     

            {
                const Write<Vector> write1(x_init_);

                each_write(x_init_, [](const SizeType i) -> double 
                {   
                    SizeType m = (i+1)%4; 
                    
                    if(m ==1){
                        return 3.0; 
                    }
                    else if(m==2){
                        return -1.0; 
                    }
                    else if(m==3){
                        return 0.0; 
                    }
                    else{
                        return 1.0; 
                    }

                }   );                
            }      
        }

        std::string name() const override
        {
            return "Extended Powell singular"; 
        }

        SizeType dim() const override
        {
            return n_; 
        }

        bool exact_sol_known() const override
        {
            return false;  // just because we can not fit into precision
        }
        

        bool value(const Vector &x, Scalar &result) const override 
        {
            if( mpi_world_size() > 1){
                utopia_error("Function is not supported in parallel... \n"); 
                return false; 
            }

            const SizeType n = this->dim(); 
            assert(size(x).get(0) == n);
        
            result = 0.0; 
            Scalar xjp1, xjp2, xjp3,f1, f2, f3, f4; 

            {
                Read<Vector> read(x); 

                for(SizeType j=1; j <=n; j+=4)
                {
                    if (j + 1 <= n )
                      xjp1 = x.get(j);
                    else
                      xjp1 = 0.0;

                    if (j+2 <= n )
                      xjp2 = x.get(j+1);
                    else
                      xjp2 = 0.0;

                    if (j+3 <= n )
                      xjp3 = x.get(j+2);
                    else
                      xjp3 = 0.0;

                    f1 = x.get(j-1) + (10.0 * xjp1);

                    if (j+1 <= n )
                      f2 = xjp2 - xjp3;
                    else
                      f2 = 0.0;

                    if ( j + 2 <= n )
                      f3 = xjp1 - (2.0 * xjp2);
                    else
                      f3 = 0.0;

                    if ( j + 3 <= n )
                      f4 = x.get(j-1) - xjp3;
                    else
                      f4 = 0.0;

                    result += (f1 * f1) +  (5.0 * f2 * f2) + (f3 * f3 * f3 * f3) + (10.0 * f4 * f4 * f4 * f4);
                }
            }

            return true;
        }

        bool gradient(const Vector &x, Vector &g) const override 
        {   
            if( mpi_world_size() > 1){
                utopia_error("Function is not supported in parallel... \n"); 
                return false; 
            }

            const SizeType n = this->dim(); 
            assert(size(x).get(0) == n);
            
            if(empty(g))
            {
                g = zeros(n); 
            }

            Scalar xjp1, xjp2, xjp3,f1, f2, f3, f4, f43, f33; 
            Scalar df1dxj, df1dxjp1, df2dxjp3, df2dxjp2, df3dxjp1, df3dxjp2, df4dxj, df4dxjp3; 

            {
                Read<Vector> read(x); 
                Write<Vector> write(g); 

                for(SizeType j=1; j <=n; j+=4)
                {
                    if (j + 1 <= n )
                      xjp1 = x.get(j);
                    else
                      xjp1 = 0.0;

                    if (j+2 <= n )
                      xjp2 = x.get(j+1);
                    else
                      xjp2 = 0.0;

                    if (j+3 <= n )
                      xjp3 = x.get(j+2);
                    else
                      xjp3 = 0.0;

                    f1 = x.get(j-1) + (10.0 * xjp1);
                    df1dxj = 1.0;
                    df1dxjp1 = 10.0;

                    if (j+1 <= n ){
                        f2 = xjp2 - xjp3;
                        df2dxjp2 = 1.0;
                        df2dxjp3 = -1.0;                      
                    }
                    else{
                        f2 = 0.0;
                        df2dxjp2 = 0.0;
                        df2dxjp3 = 0.0;
                    }

                    if ( j + 2 <= n ){
                        f3 = xjp1 - (2.0 * xjp2);
                        df3dxjp1 = 1.0;
                        df3dxjp2 = -2.0;
                    }
                    else{
                        f3 = 0.0;
                        df3dxjp1 = 1.0;
                        df3dxjp2 = -2.0;
                    }

                    if ( j + 3 <= n ){
                        f4 = x.get(j-1) - xjp3;
                        df4dxj = 1.0;
                        df4dxjp3 = -1.0;
                    }
                    else{
                        f4 = 0.0;
                        df4dxj = 0.0;
                        df4dxjp3 = 0.0;
                    }

                    f43 = f4*f4*f4; 
                    f33 = f3*f3*f3;

                    g.set(j-1,  (2.0 * f1 * df1dxj) + (10.0 * 4.0 * f43 * df4dxj));


                    if ( j+1 <= n ){
                      g.set(j, (2.0 * f1 * df1dxjp1) + (4.0 * f33 * df3dxjp1));
                    }

                    if ( j+2 <= n ){
                      g.set(j+1, (2.0 * 5.0 * f2 * df2dxjp2) + (4.0 * f33 * df3dxjp2));
                    }

                    if ( j+3 <= n ){
                      g.set(j+2, (2.0 * 5.0 * f2 * df2dxjp3) + (10.0 * 4.0 * f43 * df4dxjp3));
                    }
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
                        
            const SizeType n = this->dim(); 
            assert(size(x).get(0) == n);

            H = sparse(n, n, 3.0); 

            Scalar xjp1, xjp2, xjp3,f1, f2, f3, f4, f43, f33; 
            Scalar df1dxj, df1dxjp1, df2dxjp3, df2dxjp2, df3dxjp1, df3dxjp2, df4dxj, df4dxjp3; 

            {
                Read<Vector> read(x); 
                Write<Matrix> write(H); 

                for(SizeType j=1; j <=n; j+=4)
                {
                    if (j + 1 <= n )
                      xjp1 = x.get(j);
                    else
                      xjp1 = 0.0;

                    if (j+2 <= n )
                      xjp2 = x.get(j+1);
                    else
                      xjp2 = 0.0;

                    if (j+3 <= n )
                      xjp3 = x.get(j+2);
                    else
                      xjp3 = 0.0;

                    f1 = x.get(j-1) + (10.0 * xjp1);
                    df1dxj = 1.0;
                    df1dxjp1 = 10.0;

                    if (j+1 <= n ){
                        f2 = xjp2 - xjp3;
                        df2dxjp2 = 1.0;
                        df2dxjp3 = -1.0;                      
                    }
                    else{
                        f2 = 0.0;
                        df2dxjp2 = 0.0;
                        df2dxjp3 = 0.0;
                    }

                    if ( j + 2 <= n ){
                        f3 = xjp1 - (2.0 * xjp2);
                        df3dxjp1 = 1.0;
                        df3dxjp2 = -2.0;
                    }
                    else{
                        f3 = 0.0;
                        df3dxjp1 = 1.0;
                        df3dxjp2 = -2.0;
                    }

                    if ( j + 3 <= n ){
                        f4 = x.get(j-1) - xjp3;
                        df4dxj = 1.0;
                        df4dxjp3 = -1.0;
                    }
                    else{
                        f4 = 0.0;
                        df4dxj = 0.0;
                        df4dxjp3 = 0.0;
                    }

                    f43 = f4*f4*f4; 
                    f33 = f3*f3*f3;

                    H.set(j-1, j-1, (2.0 * df1dxj * df1dxj) + (120.0 * f4 * f4 * df4dxj * df4dxj));


                    if ( j+1 <= n )
                    {
                        H.set(j-1,j, 2.0 * df1dxj * df1dxjp1);
                        H.set(j,j-1,  2.0 * df1dxj * df1dxjp1);
                        H.set(j,j , (2.0 * df1dxjp1 * df1dxjp1) + (12.0 * f3 * f3 * df3dxjp1 * df3dxjp1));
                    }

                    if ( j+2 <= n )
                    {
                        H.set(j-1,j+1, 0.0);
                        H.set(j+1,j-1, 0.0);
                        H.set(j,j+1,  12.0 * f3 * f3 * df3dxjp2 * df3dxjp1);
                        H.set(j+1,j,  12.0 * f3 * f3 * df3dxjp1 * df3dxjp2);
                        H.set(j+1,j+1,  (10.0 * df2dxjp2 * df2dxjp2) + (12.0 * f3 * f3 * df3dxjp2 * df3dxjp2));
                    }

                    if ( j+3 <= n ){
                        H.set(j-1,j+2, 120.0 * f4 * f4 * df4dxj * df4dxjp3);
                        H.set(j+2,j-1, 120.0 * f4 * f4 * df4dxj * df4dxjp3);
                        H.set(j,j+2,  0.0);
                        H.set(j+2,j,  0.0);
                        H.set(j+1,j+2, 10.0 * df2dxjp3 * df2dxjp2);
                        H.set(j+2,j+1, 10.0 * df2dxjp2 * df2dxjp3);
                        H.set(j+2,j+2, (10.0 * df2dxjp3 * df2dxjp3) + (120.0 * f4 * f4 * df4dxjp3 * df4dxjp3)); 
                    }
                }
            }

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
        SizeType n_; 
        Vector x_init_; 
        Vector x_exact_; 

    };    
    
}

#endif //UTOPIA_SOLVER_EXTENDED_POWELL_SINGULAR_22
