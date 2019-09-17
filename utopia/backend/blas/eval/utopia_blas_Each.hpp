#ifndef UTOPIA_BLAS_EACH_HPP
#define UTOPIA_BLAS_EACH_HPP

#include "utopia_Each.hpp"
#include "utopia_blas_Types.hpp"
#include "utopia_blas_ForwardDeclarations.hpp"

namespace utopia {
        template<typename T, int FILL_TYPE>
        class Each<BlasVector<T>, 1, FILL_TYPE> {
        public:
            using BlasVector = utopia::BlasVector<T>;
            
            template<class Fun>
            inline static void apply_read(const BlasVector &v, Fun fun)
            {
                For<>::apply(
                    0,
                    v.size(),
                    [&v, &fun](const std::size_t i) {
                        fun(i, v.get(i));
                    }
                );
            }

            template<class Fun>
            inline static void apply_write(BlasVector &v, Fun fun)
            {
                For<>::apply(
                    0,
                    v.size(),
                    [&v, &fun](const std::size_t i) {
                        v.set(i, fun(i));
                    }
                );
            }

            template<class Fun>
            inline static void apply_transform(const BlasVector &in, BlasVector &out, Fun fun)
            {
                auto s = in.size();
                if(s != out.size()) {
                    out.resize(s);
                }

                For<>::apply(
                    0,
                    s,
                    [&in, &out, &fun](const std::size_t i) {
                        out.set(i, fun(i, in.get(i)));
                    }
                );
            }
        };
}

#endif //UTOPIA_BLAS_EACH_HPP
