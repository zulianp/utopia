#ifndef UTOPIA_STATIC_MATH_HPP
#define UTOPIA_STATIC_MATH_HPP

#include "utopia_Base.hpp"

namespace utopia {

    template <int N>
    class Factorial {
    public:
        static const int value = Factorial<N - 1>::value * N;
    };

    template <>
    class Factorial<1> {
    public:
        static const int value = 1;
    };

    template <>
    class Factorial<0> {
    public:
        static const int value = 1;
    };

    template <>
    class Factorial<-1> {
    public:
        static const int value = 1;
    };

    template <int Base, int Expon>
    class Power {
    public:
        static const int value = Base * Power<Base, Expon - 1>::value;
    };

    template <int Base>
    class Power<Base, 1> {
    public:
        static const int value = Base;
    };

    template <int Base>
    class Power<Base, 0> {
    public:
        static const int value = 1;
    };

    template <int N, int K>
    class CombinationsAux {
    public:
        static const int NChooseK = Factorial<N>::value / (Factorial<K>::value * Factorial<N - K>::value);

        static void apply(std::array<std::array<int, K>, NChooseK> &combs) {
            std::array<int, K> data;
            int comb_index = 0;
            apply(data, 0, 0, combs, comb_index);
        }

    private:
        static void apply(std::array<int, K> &data,
                          const int index,
                          const int i,
                          std::array<std::array<int, K>, NChooseK> &combs,
                          int &comb_index) {
            if (index == K) {
                std::copy(std::begin(data), std::end(data), std::begin(combs[comb_index++]));
                return;
            }

            if (i >= N) {
                return;
            }

            data[index] = i;

            apply(data, index + 1, i + 1, combs, comb_index);

            // current is excluded, replace it with next (Note that
            // i+1 is passed, but index is not changed)
            apply(data, index, i + 1, combs, comb_index);
        }
    };

    template <int N, int ChooseM>
    class Combinations {
    public:
        static const int value = Factorial<N>::value / (Factorial<ChooseM>::value * Factorial<N - ChooseM>::value);
        std::array<std::array<int, ChooseM>, value> combs;

        static void print_all() {
            for (auto const &c : instance().combs) {
                for (auto n : c) {
                    std::cout << n << " ";
                }

                std::cout << std::endl;
            }
        }

        template <typename T>
        static void choose(const int k, const std::array<T, N> &in, std::array<T, ChooseM> &out) {
            assert(k < value);
            const auto &comb = instance().combs[k];

            for (int i = 0; i < ChooseM; ++i) {
                assert(comb[i] < N);
                assert(comb[i] >= 0);

                out[i] = in[comb[i]];
            }
        }

        static void generate(const int k, int comb[ChooseM]) {
            std::copy(instance().combs[k].begin(), instance().combs[k].end(), comb);
        }

    private:
        Combinations() { CombinationsAux<N, ChooseM>::apply(combs); }

        inline static const Combinations &instance() {
            static const Combinations instance_;
            return instance_;
        }
    };

    template <int N>
    class Combinations<N, 1> {
    public:
        static const int value = N;
        static void generate(const int k, int comb[1]) { comb[0] = k; }

        template <typename T>
        static void choose(const int k, const std::array<T, N> &in, std::array<T, 1> &out) {
            assert(k < N);
            out[0] = in[k];
        }
    };

}  // namespace utopia

#endif  // UTOPIA_STATIC_MATH_HPP
