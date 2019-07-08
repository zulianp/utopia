#ifndef UTOPIA_STABILIZE_TRANSPORT_HPP
#define UTOPIA_STABILIZE_TRANSPORT_HPP

#include "utopia_fe_base.hpp"
#include "utopia_Traits.hpp"

#include <cmath>

namespace utopia {

    //A0 is the coupled transport term
    inline static void transport_stabilization(
        const USparseMatrix &A, 
        USparseMatrix &S)
    { 
        using SizeType = Traits<USparseMatrix>::Scalar;
        using Scalar   = Traits<USparseMatrix>::SizeType;

        USparseMatrix A_t = transpose(A);
        auto ls = local_size(A);
        S = A;
        S *= 0;

        {
            Read<USparseMatrix>  r(A_t);
            Write<USparseMatrix> w(S);

            each_read(A, [&](const SizeType i, const SizeType j, Scalar value){
                if(i != j) {
                    const Scalar value_t = A_t.get(i, j);
                    Scalar max_val = std::max(value, value_t);

                    if(max_val > 0.0){
                        max_val *= -1.0;

                        S.set(i, j, max_val);
                    }
                }
            });
        }


        UVector diag_elem = -1.0 * sum(S, 1);
        USparseMatrix S_diag = diag(diag_elem);    
        S += S_diag;
    }


    // auto op = std::make_shared<Factorization<USparseMatrix, UVector> >(MATSOLVERMUMPS,PCLU);

    // USparseMatrix D_b = diag(sum(B, 1));

    // UVector diag_elem_b = 1./sum((B),1);

    // USparseMatrix Dinv = diag(diag_elem_b);

    // USparseMatrix T =  Dinv * B;

    // const_cast<StoreTransferOperators&>(_fe_problem.getUserObject<StoreTransferOperators>(_userobject_name_T)).setTransferOperator()= std::make_shared<USparseMatrix>(T);

    // USparseMatrix A_m_tot_cons = transpose(T) * A_f_t * T + A_m_t;

    // USparseMatrix S_m_cons = A_m_t;

    // S_m_cons*=0;

    // stabilize_A_matrix_cons(A_m_tot_cons, S_m_cons);


    // const_cast<StoreTransferOperators&>(_fe_problem.getUserObject<StoreTransferOperators>(_userobject_name_1)).setTransferOperator()= std::make_shared<USparseMatrix>(S_m_cons);


    // USparseMatrix A_sc = mass_lumped_mp * inv_dt + transpose(T) * mass_lumped_fp * inv_dt * T + A_m_tot_cons + S_m_cons;

    // UVector rhs_sc     = rhs_m_t + transpose(T) * rhs_f_t;

    // constraint_concentration_mat(rhs_m_c,  A_sc, _constraint_m);

    // constraint_concentration_vec(rhs_m_c,  rhs_sc, _constraint_m);

}

#endif //UTOPIA_STABILIZE_TRANSPORT_HPP