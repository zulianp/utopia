#ifndef UTOPIA_INTREPID2_SUBDOMAIN_FUNCTION_HPP
#define UTOPIA_INTREPID2_SUBDOMAIN_FUNCTION_HPP

#include "utopia_Input.hpp"
#include "utopia_intrepid2_FE.hpp"

#include <Kokkos_UnorderedMap.hpp>

namespace utopia {
    namespace intrepid2 {
        template <typename Scalar>
        class SubdomainValue : public Configurable {
        public:
            using ExecutionSpace = typename FE<Scalar>::ExecutionSpace;
            using UnorderedMap = ::Kokkos::UnorderedMap<int, Scalar>;
            using HostUnorderedMap = typename UnorderedMap::HostMirror;
            static const int DEFAULT_MAP_CAPACITY = 10;

            void read(Input &in) override {
                in.get_all([this](Input &sub_is) {
                    int block_id = -1;
                    std::string value;
                    sub_is.get("value", value);
                    sub_is.get("block", block_id);

                    if (!value.empty()) {
                        host_map.insert(block_id, std::atof(value.c_str()));
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
        };
    }  // namespace intrepid2
}  // namespace utopia

#endif  // UTOPIA_INTREPID2_SUBDOMAIN_FUNCTION_HPP
