#ifndef UTOPIA_LINEAR_ELASTICITY_HPP
#define UTOPIA_LINEAR_ELASTICITY_HPP

#include "utopia_ElasticMaterial.hpp"
#include "ui/utopia_UIScalarSampler.hpp"
#include "utopia_ElementWisePseudoInverse.hpp"
#include "utopia_fe_core.hpp"
#include "utopia_libmesh.hpp"
#include "utopia_FEFilter.hpp"
#include "utopia_FEEval_Filter.hpp"
#include "utopia_VonMisesStress.hpp"


#include <libmesh/tensor_value.h>
#include <libmesh/fe.h>

namespace utopia {

    template<class FunctionSpaceT, class Matrix, class Vector>
    class LinearElasticity final : public ElasticMaterial<Matrix, Vector> {
    public:
        LinearElasticity(FunctionSpaceT &V, const LameeParameters &params)
        : V_(V), params_(params), initialized_(false)
        {}

        bool init(Matrix &hessian)
        {
            if(initialized_) return true;
            initialized_ = assemble_hessian(hessian);
            return initialized_;
        }

        void clear() override
        {
            initialized_ = false;
        }

        bool is_linear() const override { return true; }

        
        // bool assemble_hessian_and_gradient(const Vector &x, Matrix &hessian, Vector &gradient) override
        bool assemble_hessian_and_gradient(const Vector &x, Matrix &hessian, Vector &gradient) override
        {
            if(!init(hessian)) {
                return false;
            }

            gradient = hessian * x;
            return true;
        }

        bool stress(const Vector &x, Vector &result) override {
            Matrix hessian;

            if(!assemble_hessian(hessian)) {
                return false;
            }

            result = hessian * x;
            return true;
        }

        bool normal_stress(const UVector &x, UVector &out)
        {
            auto u = trial(V_);
            auto vx = test(V_[0]);
            
            auto mu     = params_.var_mu();
            auto lambda = params_.var_lambda();

            auto uk = interpolate(x, u);

            auto strain = transpose(grad(uk)) + grad(uk);
            auto stress = mu * strain + lambda * trace(strain) * identity();
            auto normal_stress = dot(normal(), stress * normal());
            
            UVector mass_vector;
            bool ok = utopia::assemble(surface_integral(dot(normal_stress, vx)), out); assert(ok);
            if(!ok) return false;

            utopia::assemble(surface_integral(dot(coeff(1.0), vx)), mass_vector);

            e_pseudo_inv(mass_vector, mass_vector, 1e-14);
            out = e_mul(mass_vector, out);
            return true;
        }  

        bool von_mises_stress(const UVector &x, UVector &out)
        {
            auto u = trial(V_);
            auto vx = test(V_[0]);
            
            auto mu     = params_.var_mu();
            auto lambda = params_.var_lambda();

            auto uk = interpolate(x, u);

            auto strain = transpose(grad(uk)) + grad(uk);
            auto stress = mu * strain + lambda * trace(strain) * identity();

            auto vm = filter(stress, VonMisesStress::apply);

            UVector mass_vector;
            bool ok = utopia::assemble(integral(dot(vm, vx)), out); assert(ok);
            if(!ok) return false;

            utopia::assemble(integral(dot(coeff(1.0), vx)), mass_vector);

            e_pseudo_inv(mass_vector, mass_vector, 1e-14);
            out = e_mul(mass_vector, out);
            return true;
        }

    private:
        FunctionSpaceT &V_;
        LameeParameters params_;
        bool initialized_;

        bool assemble_hessian(Matrix &hessian)
        {
            auto u = trial(V_);
            auto v = test(V_);

            auto mu     = params_.var_mu();
            auto lambda = params_.var_lambda();

            auto e_u = 0.5 * ( transpose(grad(u)) + grad(u) );
            auto e_v = 0.5 * ( transpose(grad(v)) + grad(v) );

            auto b_form = integral((2. * mu) * inner(e_u, e_v) + lambda * inner(div(u), div(v)));

            return assemble(b_form, hessian);
        }

        // static double von_mises_stress_2(const double *stress)
        // {
        //     using std::sqrt;

        //     double result =  0.5 * ( stress[0] - stress[3] ) *
        //     ( stress[0] - stress[3] ) +
        //     3.0  *  stress[1] * stress[1];
            
        //     result = sqrt( fabs(result) );
        //     assert(result == result && "von_mises_stress_2: result is nan");
        //     return result;
        // }

        // static double von_mises_stress_3(const double *stress)
        // {
        //     using std::sqrt;

        //     double result =  0.5 * ( stress[0] - stress[4] ) *
        //     ( stress[0] - stress[4] ) +
        //     3.0  *  stress[1] * stress[1];
            
        //     result += 0.5 * (stress[8] - stress[4]) * (stress[8] - stress[4]) + 3.0  * stress[7] * stress[7];
        //     result += 0.5 * (stress[8] - stress[0]) * (stress[8] - stress[0]) + 3.0  * stress[6] * stress[6];
            
        //     result = sqrt( fabs(result) );
            
        //     assert(result == result && "von_mises_stress_3: result is nan");
        //     return result;
        // }

        // static double von_mises_stress(const int n_dims, const double * stress)
        // {
        //     switch(n_dims) {
        //         case 2: { return von_mises_stress_2(stress); }
        //         case 3: { return von_mises_stress_3(stress); }
        //         default : { assert(false && "von_mises_stress: not supported for dim."); return 0.0; }
        //     }
        //     return 0.0;
        // }

        // static void stress_linear_elasticity(const double mu, const double lambda, const libMesh::DenseMatrix<double> &grad_u, libMesh::DenseMatrix<double> &stress)
        // {
        //     stress = grad_u;
        //     const int n = stress.m();
            
        //     double trace_grad_u = 0;
            
        //     for(int i = 0; i < n; ++i) {
        //         trace_grad_u += grad_u(i, i);
                
        //         for(int j = 0; j < n; ++j) {
        //             stress(i, j) += grad_u(j, i);
        //         }
        //     }
            
        //     stress *= mu;
        //     double temp = (lambda * trace_grad_u);
        //     for(int i = 0; i < n; ++i) {
        //         stress(i, i) += temp;
        //     }
        // }

        // template<class FE>
        // static void von_mises_stress_linear_elasticity(FE &fe, const int dims, const double mu, const double lambda, const libMesh::DenseVector<double> &u,
        //                                         libMesh::DenseVector<double> &von_mises_stress_vec,
        //                                         libMesh::DenseVector<double> &mass_vec)
        // {
        //     auto &fun = fe.get_fe().get_phi();
        //     auto &grad_fun = fe.get_fe().get_dphi();
        //     auto &JxW   = fe.get_fe().get_JxW();
            
        //     von_mises_stress_vec.resize(fun.size());
        //     von_mises_stress_vec.zero();
            
        //     mass_vec.resize(fun.size());
        //     mass_vec.zero();
            
        //     double vm_stress = 0.0;
        //     double mass = 0.0;
            
        //     libMesh::DenseMatrix<double> grad_u;
        //     libMesh::DenseMatrix<double> stress;
        //     for(SizeType qp = 0; qp < JxW.size(); ++qp) {
        //         libMesh::DenseMatrix<double> grad_u(dims, dims);
        //         grad_u.zero();
                
        //         for(SizeType i = 0; i < grad_fun.size(); ++i) {
        //             for(SizeType di = 0; di < dims;  ++di) {
        //                 for(SizeType dj = 0; dj < dims; ++dj) {
        //                     grad_u(di, dj) += grad_fun[i][qp](di, dj) * u(i);
        //                 }
        //             }
                    
        //         }
                
        //         stress_linear_elasticity(mu, lambda, grad_u, stress);
        //         //von mises
                
        //         for(SizeType i = 0; i < fun.size(); ++i) {
        //             auto FxJxW = fun[i][qp] * von_mises_stress(dims, &stress.get_values()[0]) * JxW[qp];
        //             auto MxJxW = fun[i][qp] * JxW[qp];
                    
        //             for(SizeType d = 0; d < dims; ++d) {
        //                 von_mises_stress_vec(i) += FxJxW(d);
        //                 mass_vec(i) += MxJxW(d);
        //             }
        //         }
        //     }
        // }


        // static void assemble_linear_elasticity(libMesh::FEBase &fe, const double lambda, const double mu, libMesh::DenseMatrix<double> &mat)
        // {
        //     typedef libMesh::DenseMatrix<double> DenseMatrixT;
        //     typedef libMesh::TensorValue<double> TensorValueT;
        //     typedef libMesh::DenseVector<double> DenseVectorT;
        //     typedef unsigned int uint;
            
        //     auto &grad  = fe.get_dphi();
        //     auto &JxW   = fe.get_JxW();
        //     auto &div   = fe.get_div_phi();
            
        //     uint n_quad_points = grad[0].size();
            
        //     mat.resize(grad.size(), grad.size());
        //     mat.zero();
            
        //     std::vector<TensorValueT> strain(grad.size());
            
        //     for (uint qp = 0; qp < n_quad_points; qp++) {
        //         //precompute straint tensor for each quadrature point
        //         for(uint i = 0; i < grad.size(); ++i) {
        //             strain[i]  = grad[i][qp];
        //             strain[i] += grad[i][qp].transpose();
        //         }
                
        //         for (uint i = 0; i < strain.size(); i++) {
        //             for (uint j = i; j < strain.size(); j++) {
        //                 mat(i, j) += ( mu * 0.5 * strain[i].contract(strain[j]) + lambda * div[i][qp] * div[j][qp] ) * JxW[qp];
        //             }
        //         }
        //     }
            
        //     //exploit symmetry
        //     for(uint i = 0; i < mat.n(); ++i) {
        //         for(uint j = i+1; j < mat.n(); ++j) {
        //             mat(j, i) = mat(i, j);
        //         }
        //     }
        // }

        // template<class FE>
        // static void assemble_linear_elasticity_matrix(
        //                                        const LameeParameters &params,
        //                                        FE &fe,
        //                                        Matrix &mat)
        // {
        //     using namespace libMesh;
            
        //     auto e_begin = fe.mesh().active_local_elements_begin();
        //     auto e_end   = fe.mesh().active_local_elements_end();
            
        //     std::vector<dof_id_type> dof_indices;
            
        //     DenseMatrix<Real> el_mat;
        //     for(auto e_it = e_begin; e_it != e_end; ++e_it) {
        //         fe.set_element(**e_it);
                
        //         const int block_id = (*e_it)->subdomain_id();
                
        //         const double mu     = params.mu(block_id);
        //         const double lambda = params.lambda(block_id);
                
        //         assemble_linear_elasticity(fe, lambda, mu, el_mat);
                
        //         fe.dof_map().dof_indices(*e_it, dof_indices, fe.var_num());
        //         add_matrix(el_mat, dof_indices, dof_indices, mat);
        //     }
        // }

        // template<class FE>
        // static void assemble_von_mises_stress(
        //                                const LameeParameters &params,
        //                                FE &fe,
        //                                const Vector &u,
        //                                Vector &stress)
        // {
        //     using namespace libMesh;
            
            
        //     auto e_begin = fe.mesh().active_local_elements_begin();
        //     auto e_end   = fe.mesh().active_local_elements_end();
            
        //     std::vector<dof_id_type> dof_indices;
            
        //     Vector mass = local_zeros(local_size(u));
        //     stress = local_zeros(local_size(u));
            
        //     Read<Vector> r_u(u);
        //     {
        //         Write<Vector> w_s(stress), w_m(mass);
                
        //         DenseVector<Real> u_local, stress_local, mass_local;
        //         for(auto e_it = e_begin; e_it != e_end; ++e_it) {
        //             fe.set_element(**e_it);
                    
        //             const int block_id = (*e_it)->subdomain_id();
                    
        //             const double mu     = params.mu(block_id);
        //             const double lambda = params.lambda(block_id);
                    
        //             fe.dof_map().dof_indices(*e_it, dof_indices, fe.var_num());
                    
        //             get_vector(u, dof_indices, u_local);
                    
        //             von_mises_stress_linear_elasticity(fe, fe.mesh().mesh_dimension(), mu, lambda, u_local, stress_local, mass_local);
        //             add_vector(stress_local, dof_indices, stress);
        //             add_vector(mass_local, dof_indices, mass);
        //         }
        //     }
            
        //     stress = e_mul(stress, 1./mass);
        // }

    };
}

#endif //UTOPIA_LINEAR_ELASTICITY_HPP
