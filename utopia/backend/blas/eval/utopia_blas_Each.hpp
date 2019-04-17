#ifndef UTOPIA_BLAS_EACH_HPP
#define UTOPIA_BLAS_EACH_HPP

#include "utopia_Each.hpp"
#include "utopia_blas_Types.hpp"

namespace utopia {
        template<int FILL_TYPE>
        class Each<Vectord, 1, FILL_TYPE> {
        public:
            template<class Fun>
            inline static void apply_read(const Vectord &v, Fun fun)
            {
                const auto &impl = raw_type(v);

                For<>::apply(
                    0,
                    impl.size(),
                    [&impl, &fun](const std::size_t i) {
                        fun(i, impl[i]);
                    }
                );
            }

            template<class Fun>
            inline static void apply_write(Vectord &v, Fun fun)
            {
                auto &impl = raw_type(v);

                For<>::apply(
                    0,
                    impl.size(),
                    [&impl, &fun](const std::size_t i) {
                        impl[i] = fun(i);
                    }
                );
            }

            template<class Fun>
            inline static void apply_transform(const Vectord &in, Vectord &out, Fun fun)
            {
                const auto &impl_in = raw_type(in);
                auto &impl_out = raw_type(out);
                auto s = impl_in.size();

                if(s != impl_out.size()) {
                    impl_out.resize(s);
                }

                For<>::apply(
                    0,
                    s,
                    [&impl_in, &impl_out, &fun](const std::size_t i) {
                        impl_out[i] = fun(i, impl_in[i]);
                    }
                );
            }
        };
}

#endif //UTOPIA_BLAS_EACH_HPP
