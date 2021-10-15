#ifndef UTOPIA_KOKKOS_SUBDOMAIN_VALUE_HPP
#define UTOPIA_KOKKOS_SUBDOMAIN_VALUE_HPP

#include "utopia_Input.hpp"
#include "utopia_kokkos_FE.hpp"

#include <Kokkos_UnorderedMap.hpp>

namespace utopia {
    namespace kokkos {

        template <class FE>
        class SubdomainValue /*: public Configurable*/ {
        public:
            using Scalar = typename FE::Scalar;
            using ExecutionSpace = typename FE::ExecutionSpace;
            using UnorderedMap = ::Kokkos::UnorderedMap<int, Scalar>;
            using HostUnorderedMap = typename UnorderedMap::HostMirror;
            static const int DEFAULT_MAP_CAPACITY = 10;

            void read(Input &in) /*override*/ {
                in.get_all([this](Input &sub_is) {
                    int block_id = -1;
                    std::string value;
                    sub_is.get("value", value);
                    sub_is.get("block", block_id);

                    if (!value.empty()) {
                        host_map.insert(block_id, std::atof(value.c_str()));
                        empty = false;
                    }
                });

                map = UnorderedMap(host_map.size());
                ::Kokkos::deep_copy(map, host_map);
            }

            UTOPIA_INLINE_FUNCTION Scalar value(const int block_id) const {
                auto idx = map.find(block_id);
                if (map.valid_at(idx)) {
                    return map.value_at(idx);
                } else {
                    return default_value;
                }
            }

            SubdomainValue(Scalar default_value, int map_capacity = DEFAULT_MAP_CAPACITY)
                : default_value(default_value), host_map(map_capacity) {}

            Scalar default_value;
            UnorderedMap map;
            HostUnorderedMap host_map;
            bool empty{true};
        };

        template <class FE, class Kernel>
        class SubdomainValueTimesKernel {
        public:
            using Scalar = typename FE::Scalar;

            template <class... Args>
            UTOPIA_INLINE_FUNCTION Scalar operator()(const int block_id, Args... args) {
                auto value = subdomain_value.value(block_id);
                return value * kernel(args...);
            }

            UTOPIA_INLINE_FUNCTION SubdomainValueTimesKernel(const SubdomainValue<FE> &subdomain_value,
                                                             const Kernel &kernel)
                : subdomain_value(subdomain_value), kernel(kernel) {}

            SubdomainValue<FE> subdomain_value;
            Kernel kernel;
        };

    }  // namespace kokkos
}  // namespace utopia

#endif  // UTOPIA_KOKKOS_SUBDOMAIN_VALUE_HPP
