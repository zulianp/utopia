#include "utopia_InputParameters.hpp"
#include "utopia_MPI.hpp"
#include "utopia_ui.hpp"

namespace utopia {

    void InputParameters::init(const int argc, char *argv[]) {
        std::string key;
        std::string value;
        std::size_t key_len = 0;

        for (int i = 1; i < argc; i++) {
            key = argv[i];
            key_len = key.size();

            if (key.compare(0, 5, "@file") == 0 && i + 1 < argc) {
                auto pos = key.find_first_of(".", 5, 1);

                if (pos == std::string::npos) {
                    aux_root_ = open_istream(argv[i + 1]);
                    ++i;  // INCREMENT
                } else {
                    std::string node_key = key.substr(pos + 1, key.size() - pos - 1);
                    auto istr = open_istream(argv[i + 1]);

                    if (istr) {
                        nodes_[node_key] = std::move(istr);
                    }

                    ++i;  // INCREMENT
                }

            } else {
                if (key[0] == '-') {
                    if (key_len > 2 && key[1] == '-') {
                        // boolean option
                        set(key.substr(2, key_len - 2), true);
                    } else {
                        if (i + 1 < argc) {
                            value = argv[i + 1];
                            set(key.substr(1, key_len - 1), value);
                            ++i;  // INCREMENT
                        } else {
                            if (mpi_world_rank() == 0) {
                                std::cerr << "no value for key: " << key << std::endl;
                            }
                        }
                    }
                }
            }
        }

        if (mpi_world_rank() == 0) {
            describe(std::cout);
        }
    }

}  // namespace utopia
