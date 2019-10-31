// #include "utopia_NewLinearElasticity.hpp"

// namespace utopia {

//     template<class FunctionSpaceT, class Matrix, class Vector>
//     NewLinearElasticity<FunctionSpaceT, Matrix, Vector>::NewLinearElasticity(FunctionSpaceT &V, const LameeParameters &params)
//     : V_(V), params_(params), initialized_(false), rescaling_(1.0)
//     {}

//     template<class FunctionSpaceT, class Matrix, class Vector>
//     bool NewLinearElasticity<FunctionSpaceT, Matrix, Vector>::normal_stress(const UVector &x, UVector &out, const int subspace)
//     {
//         auto u  = trial(V_);
//         auto vx = test(V_[subspace]);
        
//         auto mu     = params_.var_mu();
//         auto lambda = params_.var_lambda();

//         auto uk = interpolate(x, u);

//         auto strain = transpose(grad(uk)) + grad(uk);
//         auto stress = mu * strain + lambda * trace(strain) * identity();
//         auto normal_stress = dot(normal(), stress * normal());
        
//         UVector mass_vector;
//         bool ok = utopia::assemble(surface_integral(dot(normal_stress, vx)), out); assert(ok);
//         if(!ok) return false;

//         utopia::assemble(surface_integral(dot(coeff(1.0), vx)), mass_vector);

//         e_pseudo_inv(mass_vector, mass_vector, 1e-14);
//         out = e_mul(mass_vector, out);
//         return true;
//     }  

//     template<class FunctionSpaceT, class Matrix, class Vector>
//     bool NewLinearElasticity<FunctionSpaceT, Matrix, Vector>::von_mises_stress(const UVector &x, UVector &out, const int subspace)
//     {
//         auto u = trial(V_);
//         auto vx = test(V_[subspace]);
        
//         auto mu     = params_.var_mu();
//         auto lambda = params_.var_lambda();

//         auto uk = interpolate(x, u);

//         auto strain = transpose(grad(uk)) + grad(uk);
//         auto stress = mu * strain + lambda * trace(strain) * identity();

//         auto vm = filter(stress, VonMisesStress::apply);

//         UVector mass_vector;
//         bool ok = utopia::assemble(integral(dot(vm, vx)), out); assert(ok);
//         if(!ok) return false;

//         utopia::assemble(integral(dot(coeff(1.0), vx)), mass_vector);

//         e_pseudo_inv(mass_vector, mass_vector, 1e-14);
//         out = e_mul(mass_vector, out);
//         return true;
//     }

//     template<class FunctionSpaceT, class Matrix, class Vector>
//     NewLinearElasticity<FunctionSpaceT, Matrix, Vector>::assemble_hessian(Matrix &hessian)
//     {
//         // init_bilinear_integrator();
//         // assert( this->bilinear_integrator() );
//         auto u = trial(V_);
//         auto v = test(V_);

//         auto mu     = params_.var_mu();
//         auto lambda = params_.var_lambda();

//         auto e_u = 0.5 * ( transpose(grad(u)) + grad(u) );
//         auto e_v = 0.5 * ( transpose(grad(v)) + grad(v) );

//         auto b_form = integral(((2. * rescaling_) * mu) * inner(e_u, e_v) + (rescaling_ * lambda) * inner(div(u), div(v)));
        
//         return assemble(b_form, hessian);
//         // return assemble(*this->bilinear_integrator(), hessian);
//     }

//     template class NewLinearElasticity<LibMeshFunctionSpace, USparseMatrix, UVector>;

// }
