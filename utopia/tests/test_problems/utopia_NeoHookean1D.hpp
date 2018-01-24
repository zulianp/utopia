#ifndef UTOPIA_SOLVER_GLOBALANDLOCALNEOHOOKEAN1D_ASSEMBLE_HPP
#define UTOPIA_SOLVER_GLOBALANDLOCALNEOHOOKEAN1D_ASSEMBLE_HPP

#include <vector>
#include "utopia_GLFunction.hpp"


namespace utopia {

    // Currently works for Petsc for global and Blas for local.
    template<class GlobalMatrix, class GlobalVector, class LocalMatrix, class LocalVector>
    class GlobalAndLocalNeoHookean1D  : public GLFunction<GlobalMatrix, GlobalVector, LocalMatrix, LocalVector> {

        typedef typename utopia::Traits<LocalVector>::Scalar Scalar;
        typedef typename utopia::Traits<LocalVector>::SizeType SizeType;

    public:

        GlobalAndLocalNeoHookean1D(SizeType localN=10, Scalar mu = 1.0, Scalar lambda = 1.0) :
                localN(localN), mu(mu), lambda(lambda)
        {
            // std::cout<<"GlobalAndLocalNeoHookean1D:: constructor... \n"; 


            GlobalVector dummy = local_zeros(localN);

            is_first_domain = (range(dummy).begin() == 0);
            is_last_domain = (range(dummy).end() == dummy.size().get(0));

            P = local_identity(localN + (!is_first_domain) + (!is_last_domain), localN);
            I = local_identity(localN, localN + (!is_first_domain) + (!is_last_domain));

            SizeType N = P.size().get(1);

            b = local_values(localN, 1.0);

            if(N!=localN)
            {
                Write<GlobalMatrix> write(P);
                Range mr = row_range(P);
                Range vr = range(b);

                if (!is_first_domain && !is_last_domain)
                {
                    P.set(mr.end()-1, vr.begin()-1, 1.0);
                    P.set(mr.end()-2, vr.end(), 1.0);
                }
                if(is_first_domain)
                {
                    P.set(mr.end()-1, vr.end(), 1.0);
                }
                if(is_last_domain)
                {
                    P.set(mr.end()-1, vr.begin()-1, 1.0);
                }

            }

            h = 1.0/(N-1);

            GlobalMatrix M;
            M = local_sparse(localN, localN, 3);
            {
                Write<GlobalMatrix> write(M);

                Range row_range = rowRange(M);
                SizeType row_begin = row_range.begin();
                SizeType row_end = row_range.end();

                if(1 >= row_begin && 1 < row_end)
                    M.set(1, 0, h/6.0);
                if(N-1 >= row_begin && N-1 < row_end)
                    M.set(N-1, N-1, h/3.0);

                for (SizeType i=row_begin; i<std::min(row_end, N-1); ++i)
                {
                    if (i==0)
                    {
                        M.set(0, 0, h/3.0);
                        M.set(0, 1, h/6.0);
                    }
                    else
                    {
                        M.set(i, i, 2.0 * h / 3.0);
                        M.set(i, i+1, h/6.0);
                        M.set(i+1, i, h/6.0);
                    }
                    if(i==row_begin && i)
                    {
                        M.set(i, i-1, h/6.0);
                    }
                }
            }

            GlobalVector ones = local_values(localN, 1, -0.03);

            force = M*ones;

            restrict(force, local_force);
            if(!is_first_domain)
                local_force.set(0, local_force.get(0)/2.0);

            if(!is_last_domain)
                local_force.set(local_force.size().get(0)-1, local_force.get(local_force.size().get(0)-1)/2.0);
            // std::cout<<"GlobalAndLocalNeoHookean1D:: constructor... end \n"; 

        }


        bool init(const GlobalVector &x) const override
        {
            return true; 
        }

        bool  value(const GlobalVector &point, Scalar & energy) const override
        {
            // std::cout<<"value \n"; 
            GlobalVector W = local_zeros(1);
            Scalar force_contribution = dot(force, point); // = 0;

            LocalVector distributed_point;
            restrict(point, distributed_point);

            Scalar local_W = 0.0;
            {
                Range dp_range = range(distributed_point);

                assert(dp_range.extent()>1);

                for (SizeType i=dp_range.begin();
                     i<dp_range.end()-1-(!is_last_domain);
                     ++i)
                {
                    Scalar u_prime = distributed_point.get(i+1)-distributed_point.get(i);
                    Scalar f_local = 1.0 + (u_prime)/h;

                    if (f_local < 0.0)
                    {
                       std::cerr << "Deformation determinant is negative" << std::endl;
                    }

                    Scalar log_f = log(f_local);
                    Scalar element_energy = 0.5 * mu * (f_local*f_local - 1.0) - mu * log_f + 0.5 * lambda * (log_f*log_f);

                    local_W += element_energy*h;
                }
            }

            {
                Range W_range = range(W);
                Write<GlobalVector> write(W);
                W.set(W_range.begin(), local_W);
            }

            Scalar W_sum = sum(W);

            energy =  W_sum - force_contribution;
            return true; 
        }



        bool local_value(const LocalVector &point, Scalar & energy) const override
        {
            std::vector<bool> element_computed(point.size().get(0), false);

            Scalar W = 0.0;

            {
                for (SizeType i=!is_first_domain; i<point.size().get(0)-1-(!is_last_domain); ++i)
                {

                    if(i>0 && !element_computed.at(i-1))
                    {
                        Scalar u_prime_loc =  (point.get(i) - point.get(i-1)) / h;
                        Scalar f_local = 1.0 + u_prime_loc;

                        if (f_local < 0.0)
                        {
                           std::cerr << "deformation determinant is negative in local energy" << std::endl;
                        }

                        Scalar log_f = log(f_local);
                        Scalar local_energy = 0.5 * mu * (f_local*f_local - 1.0) - mu * log_f + 0.5 * lambda * (log_f*log_f);

                        W += local_energy * h;

//                      if(is_last_domain)
//                      std::cout << local_energy*h << std::endl;

                        element_computed[i-1] = true;
                    }
                    if(i<point.size().get(0)-1 && !element_computed[i])
                    {
                        Scalar u_prime_loc =  (point.get(i+1) - point.get(i)) / h;
                        Scalar f_local = 1.0 + u_prime_loc;

                        if (f_local < 0.0)
                        {
                           std::cerr << "deformation determinant is negative" << std::endl;
                        }

                        Scalar log_f = log(f_local);
                        Scalar local_energy = 0.5 * mu * (f_local*f_local - 1.0) - mu * log_f + 0.5 * lambda * (log_f*log_f);

                        W += local_energy * h;

                        element_computed.at(i) = true;
                    }
                }
            }

            Scalar force_contribution = dot(local_force, point);
            
            energy = W - force_contribution; 

            return true; 

        }


        bool gradient(const GlobalVector &point, GlobalVector &result) const override
        {
            // std::cout<<"gradient \n"; 
            LocalVector distributed_point;
            restrict(point, distributed_point);

            result = local_zeros(localN);

            LocalVector distributed_gradient;
            restrict(result, distributed_gradient);

            {
                Range dp_range = range(distributed_point);

                for (SizeType i=dp_range.begin();
                     i<dp_range.end()-1;

                     ++i)
                {
                    Scalar f_local = 1.0 + (distributed_point.get(i+1)-distributed_point.get(i))/h;
                    Scalar res_loc = mu * (f_local - 1.0 / f_local) + lambda * log(f_local) / f_local;
                    distributed_gradient.set(i, distributed_gradient.get(i) - res_loc);
                    distributed_gradient.set(i+1, distributed_gradient.get(i+1) + res_loc);
                }
            }

            interpolate(distributed_gradient, result);
            result -= force;

            {
                Write<GlobalVector> write(result);
                if(is_first_domain)
                    result.set(0, 0.0);
            }

            return true;
        }



        bool local_gradient(const LocalVector &point, LocalVector &result) const override 
        {
            SizeType local_size = point.size().get(0);
            result = zeros(local_size);


            for (SizeType i = 0; i < local_size-1; ++i)
            {
                Scalar u_prime_loc =  (point.get(i+1) - point.get(i)) / h;
                Scalar f_local = 1.0 + u_prime_loc;
                Scalar res_loc = mu * (f_local - 1.0 / f_local) + lambda * log(f_local) / f_local;

                result.set(i, result.get(i) - res_loc);
                result.set(i+1, result.get(i+1) + res_loc);
            }

            result -= local_force;

            //BC
            if(is_first_domain)
                result.set(0, 0.0);

            if(!is_first_domain)
                result.set(0, 0.0);
            if(!is_last_domain)
                result.set(result.size().get(0)-1, 0.0);


            return true;
        }



        bool hessian(const GlobalVector &point, GlobalMatrix &result) const override 
        {
            // std::cout<<"hessian \n"; 
            result = sparse(point.size().get(0), point.size().get(0), 1);

            LocalVector distributed_point;
            restrict(point, distributed_point);

            Range row_range = rowRange(result);
            Range dp_range = range(distributed_point);

            SizeType dp_i = dp_range.begin();
            {
                Write<GlobalMatrix> write(result);
                if(!is_first_domain)
                {
                    SizeType row = row_range.begin();

                    Scalar f_local = 1.0 + (distributed_point.get(dp_i+1)-distributed_point.get(dp_i))/h;

                    LocalMatrix loc_hessian(1, 2, {-1, 1});
                    Scalar inv_h = (1.0/h);
                    LocalMatrix loc_hessian_scaled;
                    loc_hessian_scaled = loc_hessian*inv_h;
                    Scalar coefficient = mu * (1.0 + 1.0 / (f_local*f_local)) + lambda * (1 - log(f_local)) / (f_local*f_local);


                    LocalMatrix weighted_loc_hessian = coefficient * loc_hessian_scaled;


                    for (int j = -1; j<1; ++j)
                        result.set(row, row+j, weighted_loc_hessian.get(0, 1+j));

                    ++dp_i;
                }
            }

            for (SizeType row = row_range.begin(); row<row_range.end()-1; ++row, ++dp_i)
            {
                Scalar f_local = 1.0 + (distributed_point.get(dp_i+1)-distributed_point.get(dp_i))/h;

                LocalMatrix loc_hessian = identity(2, 2);
                loc_hessian.set(0,1,-1.0);
                loc_hessian.set(1,0,-1.0);
                Scalar inv_h = (1.0/h);
                LocalMatrix loc_hessian_scaled;
                loc_hessian_scaled = loc_hessian*inv_h;
                Scalar coefficient = mu * (1.0 + 1.0 / (f_local*f_local)) + lambda * (1 - log(f_local)) / (f_local*f_local);

                LocalMatrix weighted_loc_hessian = coefficient * loc_hessian_scaled;

                LocalMatrix old_loc = zeros(2,2);
                {
                    Read<GlobalMatrix> read(result);
                    for (SizeType j=0; j<2; ++j)
                        for (SizeType k=0; k<2; ++k)
                        {
                            old_loc.set(j,k, result.get(row+j, row+k));
                        }
                }

                {
                    Write<GlobalMatrix> write(result);
                    for (SizeType j=0; j<2; ++j)
                        for (SizeType k=0; k<2; ++k)
                            result.set(row+j, row+k, weighted_loc_hessian.get(j, k)+old_loc.get(j,k));
                }

            }

            LocalMatrix old_loc = zeros(1, 2);
            if(row_range.end() != result.size().get(0))
            {
                Read<GlobalMatrix> read(result);
                for (SizeType i=0; i<2; ++i)
                    old_loc.set(0, i, result.get(row_range.end()-1, row_range.end()-1+i));
            }

            {
                Write<GlobalMatrix> write(result);
                if(!is_last_domain)
                {
                    SizeType row = row_range.end()-1;

                    Scalar f_local = 1.0 + (distributed_point.get(dp_i+1)-distributed_point.get(dp_i))/h;

                    LocalMatrix loc_hessian(1, 2, {1, -1});
                    Scalar inv_h = (1.0/h);
                    LocalMatrix loc_hessian_scaled;
                    loc_hessian_scaled = loc_hessian*inv_h;
                    Scalar coefficient = mu * (1.0 + 1.0 / (f_local*f_local)) + lambda * (1 - log(f_local)) / (f_local*f_local);

                    LocalMatrix weighted_loc_hessian = coefficient * loc_hessian_scaled;

                    for (int j = 0; j<2; ++j)
                        result.set(row, row+j, weighted_loc_hessian.get(0, j)+old_loc.get(0, j));

                    ++dp_i;
                }
            }

            {
                Write<GlobalMatrix> write(result);
                if(row_range.begin() == 0)
                {
                    result.set(0, 0, 1.0);
                    for(SizeType i=1; i<result.size().get(1); ++i)
                    {
                        result.set(0, i, 0.0);
                    }
                }
            }

            return true;
        }




        bool local_hessian(const LocalVector &point, LocalMatrix &result) const override
        {
            SizeType local_size = point.size().get(0);
            result = zeros(point.size().get(0), point.size().get(0));

            for (SizeType row = 0; row<local_size-1; ++row)
            {
                Scalar f_local = 1.0 + (point.get(row+1)-point.get(row))/h;

                LocalMatrix loc_hessian = identity(2, 2);
                loc_hessian.set(0,1,-1.0);
                loc_hessian.set(1,0,-1.0);
                Scalar inv_h = (1.0/h);
                LocalMatrix loc_hessian_scaled;
                loc_hessian_scaled = loc_hessian*inv_h;
                Scalar coefficient = mu * (1.0 + 1.0 / (f_local*f_local)) + lambda * (1 - log(f_local)) / (f_local*f_local);

                auto loc_old = result.range(row, row+2, row, row+2);


                LocalMatrix tmp = loc_old;
                LocalMatrix weighted_loc_hessian = coefficient * loc_hessian_scaled;

                //TODO understand why
                LocalMatrix tmp2 = tmp+weighted_loc_hessian;

                result.range(row, row+2, row, row+2) = tmp2;

            }

            LocalMatrix zero_line = values(1, result.size().get(1), 0.0);

            // BC and overlap
            result.range(0, 1, 0, result.size().get(1)) = zero_line;
            result.set(0, 0, 1.0);

            if(!is_last_domain)
            {
                result.range(result.size().get(0)-1, result.size().get(0), 0, result.size().get(1)) = zero_line;
                result.set(result.size().get(0)-1, result.size().get(1)-1, 1.0);
            }

            return true;
        }




        /// As a very simple restriction, we simply project the global vector on a local vector, rearranging the components according to our
        /// distribution.
        bool restrict(const GlobalVector &global_vector, LocalVector &local_vector) const override 
        {

            // std::cout<<"restrict 1 ... \n"; 
            GlobalVector distributed_vector = P * global_vector;

            // std::cout<<"restrict 2 ... \n"; 
            Range xr = range(distributed_vector);
            const SizeType xb = xr.begin();
            const SizeType local_extent = xr.extent();

            // std::cout<<"restrict 3 ... \n"; 
            local_vector = zeros(local_extent);
            // std::cout<<"restrict 4 ... \n"; 
            {
                Read<GlobalVector> read(distributed_vector);
                for (SizeType i = 0; i < local_extent-1+is_first_domain; ++i)
                    local_vector.set(i+(!is_first_domain), distributed_vector.get(xb+i));
                if(!is_first_domain)
                    local_vector.set(0, distributed_vector.get(xr.end()-1));
            }
            // std::cout<<"restrict end... \n"; 

            return true;
        }


        /// As a very simple interpolation, we simply project the local vector on a global vector, rearranging.
        bool interpolate(const LocalVector &local_vector,  GlobalVector &global_vector) const override 
        {
            // std::cout<<"interpolate ... \n"; 
            GlobalVector interpolated = local_zeros(localN+(!is_first_domain)+(!is_last_domain));
            Range xr = range(interpolated);
            const SizeType xb = xr.begin();
            const SizeType local_extent = xr.extent();

            {
                Write<GlobalVector> write(interpolated);
                for (SizeType i = 0; i < local_extent-1+is_first_domain; ++i)
                    interpolated.set(i+xb, local_vector.get(i+(!is_first_domain)));
                if(!is_first_domain)
                    interpolated.set(xb+local_extent-1, local_vector.get(0));
            }

            global_vector = I*interpolated;

            // std::cout<<"interpolate end... \n"; 
            return true;
        }


        bool restrict(const GlobalMatrix &, LocalMatrix &) const override 
        {
            std::cout<<"not yet... \n"; 
            return false; 
        }


        bool interpolate(const LocalMatrix &, GlobalMatrix &) const override 
        {
            std::cout<<"not yet... \n"; 
            return false; 
        }


    private:
        const SizeType localN;
        Scalar mu;
        Scalar lambda;
        GlobalVector force;
        GlobalVector b;
        Scalar h;
        GlobalMatrix P;
        GlobalMatrix I;
        LocalVector local_force;
        bool is_first_domain;
        bool is_last_domain;
    };
}

#endif //UTOPIA_SOLVER_GLOBALANDLOCALNEOHOOKEAN1D_ASSEMBLE_HPP
