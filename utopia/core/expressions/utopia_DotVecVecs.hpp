#ifndef UTOPIA_EVAL_DOT_VEC_VECS_HPP
#define UTOPIA_EVAL_DOT_VEC_VECS_HPP

#include "utopia_Eval_Empty.hpp"
#include "utopia_ForwardDeclarations.hpp"
#include "utopia_Tensor.hpp"

namespace utopia
{
    template<class Vector, int Backend = Traits<Vector>::Backend>
    class EvalDots
    {
        public:
            static void apply(const Vector &v1, const std::vector<std::shared_ptr<Vector> > &vectors, std::vector<typename utopia::Traits<Vector>::Scalar> & results)
            {
                typename utopia::Traits<Vector>::SizeType n =  vectors.size();

                if(n!=results.size()){
                    results.resize(n);
                }

                for(auto i = 0; i < n; i++)
                    results[i] = dot(v1, *(vectors[i]));
            }

            static void apply(const Vector &v1, const std::vector<Vector> &vectors, std::vector<typename utopia::Traits<Vector>::Scalar> & results)
            {
                typename utopia::Traits<Vector>::SizeType n =  vectors.size();

                if(n!=results.size()){
                    results.resize(n);
                }

                for(auto i = 0; i < n; i++)
                    results[i] = dot(v1, (vectors[i]));
            }

            static void apply(const Vector &v11, const Vector &v12, typename utopia::Traits<Vector>::Scalar & result1, const Vector &v21, const Vector &v22, typename utopia::Traits<Vector>::Scalar & result2)
            {
                result1 = dot(v11, v12);
                result2 = dot(v21, v22);
            }  

            static void apply(const Vector &v11, const Vector &v12, typename utopia::Traits<Vector>::Scalar & result1, const Vector &v21, const Vector &v22, typename utopia::Traits<Vector>::Scalar & result2, const Vector &v31,const Vector &v32, typename utopia::Traits<Vector>::Scalar & result3)
            {
                result1 = dot(v11, v12);
                result2 = dot(v21, v22);
                result3 = dot(v31, v32);
            }                        
    };


    template<class Vector, int Backend = Traits<Vector>::Backend>
    class EvalNorm2s
    {
        public:
            static void apply(const Tensor<Vector, 1> &v1, const Tensor<Vector, 1> &v2, typename utopia::Traits<Vector>::Scalar & result1, typename utopia::Traits<Vector>::Scalar & result2)
            {
                result1 = norm2(v1);
                result2 = norm2(v2);
            }

            static void apply(const Tensor<Vector, 1> &v1, const Tensor<Vector, 1> &v2, const Tensor<Vector, 1> &v3, typename utopia::Traits<Vector>::Scalar & result1, typename utopia::Traits<Vector>::Scalar & result2, typename utopia::Traits<Vector>::Scalar & result3)
            {
                result1 = norm2(v1);
                result2 = norm2(v2);
                result3 = norm2(v3);
            }
    };

    template<class Vector>
    void dots(const Tensor<Vector, 1> &v1, const std::vector<std::shared_ptr<Vector> > &vectors, std::vector<typename utopia::Traits<Vector>::Scalar> & results)
    {
        EvalDots<Vector>::apply(v1.derived(), vectors, results);
    }

    template< class Vector>
    void dots(const Tensor<Vector, 1> &v1, const std::vector<Vector> &vectors, std::vector<typename utopia::Traits<Vector>::Scalar> & results)
    {
        EvalDots<Vector>::apply(v1.derived(), vectors, results);
    }


    template< class Vector>
    void dots(const Tensor<Vector, 1> &v11,const Tensor<Vector, 1> &v12, typename utopia::Traits<Vector>::Scalar & result1, const Tensor<Vector, 1> &v21,const Tensor<Vector, 1> &v22, typename utopia::Traits<Vector>::Scalar & result2)
    {
        EvalDots<Vector>::apply(v11.derived(), v12.derived(), result1, v21.derived(), v22.derived(), result2);
    }


    template< class Vector>
    void dots(const Tensor<Vector, 1> &v11,const Tensor<Vector, 1> &v12, typename utopia::Traits<Vector>::Scalar & result1, const Tensor<Vector, 1> &v21,const Tensor<Vector, 1> &v22, typename utopia::Traits<Vector>::Scalar & result2, const Tensor<Vector, 1> &v31,const Tensor<Vector, 1> &v32, typename utopia::Traits<Vector>::Scalar & result3)
    {
        EvalDots<Vector>::apply(v11.derived(), v12.derived(), result1, v21.derived(), v22.derived(), result2, v31.derived(), v32.derived(), result3);
    }


    template< class Vector>
    void norms2(const Tensor<Vector, 1> &v1,const Tensor<Vector, 1> &v2, typename utopia::Traits<Vector>::Scalar & result1, typename utopia::Traits<Vector>::Scalar & result2)
    {
        EvalNorm2s<Vector>::apply(v1.derived(), v2.derived(), result1, result2);
    }


    template< class Vector>
    void norms2(const Tensor<Vector, 1> &v1,const Tensor<Vector, 1> &v2, const Tensor<Vector, 1> &v3, typename utopia::Traits<Vector>::Scalar & result1, typename utopia::Traits<Vector>::Scalar & result2, typename utopia::Traits<Vector>::Scalar & result3)
    {
        EvalNorm2s<Vector>::apply(v1.derived(), v2.derived(), v3.derived(),  result1, result2, result3);
    }

}


#endif //UTOPIA_EVAL_DOT_VEC_VECS_HPP
