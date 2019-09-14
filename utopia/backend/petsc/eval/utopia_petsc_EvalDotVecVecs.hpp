#ifndef UTOPIA_EVAL_PETSC_DOT_VEC_VECS_HPP
#define UTOPIA_EVAL_PETSC_DOT_VEC_VECS_HPP

#include "utopia_Eval_Empty.hpp"
#include "utopia_ForwardDeclarations.hpp"
#include <vector>

namespace utopia
{

    template<class Vector>
    class EvalDots <Vector, PETSC>
    {
        public:
            static void apply(const Vector &v1, const std::vector<std::shared_ptr<Vector> > &vectors, std::vector<typename utopia::Traits<Vector>::Scalar> & results)
            {
                typename utopia::Traits<Vector>::SizeType n =  vectors.size();

                if(n!=static_cast<SizeType>(results.size()))
                    results.resize(n);

                std::vector<Vec> vecs(n);

                for(auto i=0; i < static_cast<SizeType>(vectors.size()); i++)
                    vecs[i]=(raw_type(*vectors[i]));

                 VecMDot(raw_type(v1), n, vecs.data(), results.data());
            }

            static void apply(const Vector &v1, const std::vector<Vector> &vectors, std::vector<typename utopia::Traits<Vector>::Scalar> & results)
            {
                std::vector<Vec> vecs;
                for(auto i=0; i < static_cast<SizeType>(vectors.size()); i++)
                {
                    if(!empty(vectors[i]))
                        vecs.push_back(raw_type(vectors[i]));
                }

                typename utopia::Traits<Vector>::SizeType n =  vecs.size();

                if(n != vectors.size() || n != results.size())
                {
                    std::vector<typename utopia::Traits<Vector>::Scalar> result_new(n);
                    VecMDot(raw_type(v1), n, vecs.data(), result_new.data());

                    for(auto i=0; i < n; i++)
                        results[i] = result_new[i];
                }
                else
                {
                    if(n!=results.size())
                        results.resize(n);

                    VecMDot(raw_type(v1), n, vecs.data(), results.data());
                }
            }

            static void apply(const Vector &v11, const Vector &v12, typename utopia::Traits<Vector>::Scalar & result1, const Vector &v21, const Vector &v22, typename utopia::Traits<Vector>::Scalar & result2)
            {
                VecDotBegin(raw_type(v11), raw_type(v12), &result1);
                VecDotBegin(raw_type(v21), raw_type(v22), &result2);

                VecDotEnd(raw_type(v11), raw_type(v12), &result1);
                VecDotEnd(raw_type(v21), raw_type(v22), &result2);
            }


            static void apply(const Vector &v11, const Vector &v12, typename utopia::Traits<Vector>::Scalar & result1, const Vector &v21, const Vector &v22, typename utopia::Traits<Vector>::Scalar & result2, const Vector &v31,const Vector &v32, typename utopia::Traits<Vector>::Scalar & result3)
            {
                VecDotBegin(raw_type(v11), raw_type(v12), &result1);
                VecDotBegin(raw_type(v21), raw_type(v22), &result2);
                VecDotBegin(raw_type(v31), raw_type(v32), &result3);

                VecDotEnd(raw_type(v11), raw_type(v12), &result1);
                VecDotEnd(raw_type(v21), raw_type(v22), &result2);
                VecDotEnd(raw_type(v31), raw_type(v32), &result3);
            }      

    };




    template<class Vector>
    class EvalNorm2s<Vector, PETSC>
    {
        public:
            static void apply(const Vector &v1, const Vector &v2, typename utopia::Traits<Vector>::Scalar & result1, typename utopia::Traits<Vector>::Scalar & result2)
            {
                VecNormBegin(raw_type(v1), NORM_2, &result1);
                VecNormBegin(raw_type(v2), NORM_2, &result2);
                PetscCommSplitReductionBegin(PetscObjectComm((PetscObject)raw_type(v1)));
                VecNormEnd(raw_type(v1), NORM_2, &result1);
                VecNormEnd(raw_type(v2), NORM_2, &result2);
            }

            static void apply(const Vector &v1, const Vector &v2, const Vector &v3,  typename utopia::Traits<Vector>::Scalar & result1, typename utopia::Traits<Vector>::Scalar & result2, typename utopia::Traits<Vector>::Scalar & result3)
            {
                VecNormBegin(raw_type(v1), NORM_2, &result1);
                VecNormBegin(raw_type(v2), NORM_2, &result2);
                VecNormBegin(raw_type(v3), NORM_2, &result3);
                PetscCommSplitReductionBegin(PetscObjectComm((PetscObject)raw_type(v1)));
                VecNormEnd(raw_type(v1), NORM_2, &result1);
                VecNormEnd(raw_type(v2), NORM_2, &result2);
                VecNormEnd(raw_type(v3), NORM_2, &result3);
            }
    };





}


#endif //UTOPIA_EVAL_PETSC_DOT_VEC_VECS_HPP
