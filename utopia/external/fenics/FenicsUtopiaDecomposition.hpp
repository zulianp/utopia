/*
* @Author: alenakopanicakova
* @Date:   2016-05-26
* @Last Modified by:   alenakopanicakova
* @Last Modified time: 2016-06-06
*/

// based on cutting staff... 
#ifndef FENICS_UTOPIA_DOMAIN_DECOMPOSITION_HPP
#define FENICS_UTOPIA_DOMAIN_DECOMPOSITION_HPP


#include <utopia.hpp>


namespace utopia 
{
  template<class GlobalMatrix, class GlobalVector, class LocalMatrix, class LocalVector>
  class FenicsDecomposition :  public HorizontalDecomposition<GlobalMatrix, GlobalVector, LocalMatrix, LocalVector> 
  {

    public:

      FenicsDecomposition():   
                            HorizontalDecomposition<GlobalMatrix, GlobalVector, LocalMatrix, LocalVector>()
      {

      }


      ~FenicsDecomposition(){} 


        bool init(const GlobalVector &x) const override
        {
            VecGetLocalSize(utopia::raw_type(x), &local_dim);
            return true; 
        }

        // projecting local vector into global 
        bool interpolate(const LocalVector & x_k, GlobalVector & x) const override 
        {

            x = localZeros(local_dim);
            Range xr = range(x);
            const SizeType x_begin = xr.begin();
            const SizeType local_extent = xr.extent();

            {
                Write<GlobalVector> write(x);
                for (SizeType i = 0; i < local_extent; ++i)
                {
                    x.set(i + x_begin, x_k.get(i ));
                }
            }

            return true;
        }

        // interpolating local sub-matrixes into global matrix 
        bool interpolate(const LocalMatrix & M_k, GlobalMatrix & M) const override
        {
            GlobalMatrix I = localZeros(local_dim , local_dim );
            Range rr = rowRange(I);
            const SizeType I_begin = rr.begin();
            const SizeType I_extent = rr.extent();

            {
                Write<GlobalMatrix> write(I);

                for (SizeType i = 0; i < I_extent - 1; ++i)
                {
                    for (SizeType j = 0; j < I_extent - 1; ++j)
                    {
                        I.set(i + I_begin, j + I_begin, M_k.get(i,j));
                    }
                }
            }

            M = I ;

            return true; 

        }


        // projecting global vector into local 
        bool restrict(const GlobalVector &x, LocalVector &x_k) const override 
        {

            Range xr = range(x);
            const SizeType x_begin = xr.begin();
            const SizeType local_extent = xr.extent();

            x_k = zeros(local_extent);

            {
                Read<GlobalVector> read(x);
                for (SizeType i = 0; i < local_extent; ++i)
                {
                    x_k.set(i, x.get(x_begin + i));
                }
            }

            return true;
        }

        

        // restricting global matrix into local submatrixes 
        bool restrict(const GlobalMatrix & M, LocalMatrix & M_k) const override
        { 

            GlobalMatrix R = M;

            Range rr = rowRange(R);
            const SizeType R_begin = rr.begin();
            const SizeType R_extent = rr.extent();

            M_k = zeros(R_extent, R_extent);

            {
                Read<GlobalMatrix> w(R);

                for (SizeType i = 0; i < R_extent - 1 ; ++i)
                {
                    for (SizeType j = 0; j < R_extent - 1 ; ++j)
                    {
                        M_k.set(i , j , R.get(R_begin + i, R_begin + j));
                    }
                }
            }
            return true;
        }



  private:
    mutable PetscInt local_dim; 

  };

}

#endif //FENICS_UTOPIA_DOMAIN_DECOMPOSITION_HPP

