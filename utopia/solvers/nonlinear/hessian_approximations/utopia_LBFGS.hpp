#ifndef UTOPIA_LBFGS_HPP
#define UTOPIA_LBFGS_HPP

#include "utopia_Core.hpp"


namespace utopia
{
    // TO BE DONE.... 
    template<class SparseMatrix, class DenseMatrix, class Vector>
    class LBFGS : public HessianApproximation<SparseMatrix, Vector>
    {

        typedef UTOPIA_SCALAR(Vector)    Scalar;
        typedef UTOPIA_SIZE_TYPE(Vector) SizeType;

        public:

            LBFGS(const SizeType & m): m_(m), current_m_(0)
            {

            }

            virtual bool initialize(Function<SparseMatrix, Vector> &fun, const Vector &x) override
            {
                SizeType n = local_size(x).get(0);
                H0_ = local_identity(n, n);

                // initialize some other vectors and matrices

                this->initialized(true); 

                return true;
            }   

            virtual bool update(const Vector & /* s  */, const Vector &  /* y */ ) override
            {
                // TODO 
                return true; 
            }

            virtual bool apply_Hinv(const Vector & /* g */, Vector & /*s */) const override
            {
                // TODO 
                return true; 
            }


            virtual SparseMatrix & get_Hessian() override
            {
                std::cerr<<"--- not implemented yet---- \n"; 
                return H0_;
            } 

            virtual SparseMatrix & get_Hessian_inv() override
            {
                std::cerr<<"--- not implemented yet---- \n"; 
                return H0_;
            } 


            void set_memory_size(const SizeType & m)
            {
                m_ = m; 
            }

            SizeType get_memory_size() const 
            {
                return m_; 
            }



        private:
            // static_assert(utopia::is_sparse<SparseMatrix>::value, "BFGS does not support sparse matrices."); 
            // static_assert(!utopia::is_sparse<DenseMatrix>::value, "BFGS does not support sparse matrices."); 

            SizeType m_; // memory size 
            SizeType current_m_; // current amount of vectors in the memory 
            SparseMatrix H0_;  // identity 

        };



}

#endif //UTOPIA_LBFGS_HPP
