#ifndef UTOPIA_EVAL_DOT_VEC_VECS_HPP
#define UTOPIA_EVAL_DOT_VEC_VECS_HPP

#include "utopia_Eval_Empty.hpp"
#include "utopia_ForwardDeclarations.hpp"

namespace utopia
{
    template<class Vector, int Backend = Traits<Vector>::Backend>
    class EvalDots
    {
        public:
            static void apply(const Wrapper<Vector, 1> &v1, const std::vector<std::shared_ptr<Wrapper<Vector, 1> > > &vectors, std::vector<typename utopia::Traits<Vector>::Scalar> & results)
            {
                typename utopia::Traits<Vector>::SizeType n =  vectors.size();

                if(n!=results.size())
                    results.resize(n);

                for(auto i = 0; i < n; i++)
                    results[i] = dot(v1, *(vectors[i]));
            }

            static void apply(const Wrapper<Vector, 1> &v1, const std::vector<Wrapper<Vector, 1> > &vectors, std::vector<typename utopia::Traits<Vector>::Scalar> & results)
            {
                typename utopia::Traits<Vector>::SizeType n =  vectors.size();

                if(n!=results.size())
                    results.resize(n);

                for(auto i = 0; i < n; i++)
                    results[i] = dot(v1, (vectors[i]));
            }
    };


    template<class Vector, int Backend = Traits<Vector>::Backend>
    class EvalNorm2s
    {
        public:
            static void apply(const Wrapper<Vector, 1> &v1, const Wrapper<Vector, 1> &v2, typename utopia::Traits<Vector>::Scalar & result1, typename utopia::Traits<Vector>::Scalar & result2)
            {
                result1 = norm2(v1);
                result2 = norm2(v2);
            }

            static void apply(const Wrapper<Vector, 1> &v1, const Wrapper<Vector, 1> &v2, const Wrapper<Vector, 1> &v3, typename utopia::Traits<Vector>::Scalar & result1, typename utopia::Traits<Vector>::Scalar & result2, typename utopia::Traits<Vector>::Scalar & result3)
            {
                result1 = norm2(v1);
                result2 = norm2(v2);
                result3 = norm2(v3);
            }
    };




    template< class Vector>
    void dots(const Wrapper<Vector, 1> &v1, const std::vector<std::shared_ptr<Wrapper<Vector, 1> > > &vectors, std::vector<typename utopia::Traits<Vector>::Scalar> & results)
    {
        EvalDots<Vector>::apply(v1, vectors, results);
    }


    template< class Vector>
    void dots(const Wrapper<Vector, 1> &v1, const std::vector<Wrapper<Vector, 1> > &vectors, std::vector<typename utopia::Traits<Vector>::Scalar> & results)
    {
        EvalDots<Vector>::apply(v1, vectors, results);
    }



    template< class Vector>
    void norms2(const Wrapper<Vector, 1> &v1,const Wrapper<Vector, 1> &v2, typename utopia::Traits<Vector>::Scalar & result1, typename utopia::Traits<Vector>::Scalar & result2)
    {
        EvalNorm2s<Vector>::apply(v1, v2, result1, result2);
    }


    template< class Vector>
    void norms2(const Wrapper<Vector, 1> &v1,const Wrapper<Vector, 1> &v2, const Wrapper<Vector, 1> &v3, typename utopia::Traits<Vector>::Scalar & result1, typename utopia::Traits<Vector>::Scalar & result2, typename utopia::Traits<Vector>::Scalar & result3)
    {
        EvalNorm2s<Vector>::apply(v1, v2, v3,  result1, result2, result3);
    }



}


#endif //UTOPIA_EVAL_DOT_VEC_VECS_HPP
