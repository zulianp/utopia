#include "utopia_InputParameters.hpp"
#include "utopia_MPI.hpp"
#include "utopia_ui.hpp"

#include "utopia_Path.hpp"

#include "utopia_IOStream.hpp"

namespace utopia {

    static std::unique_ptr<Input> include(const Path &script_dir,
                                          std::unique_ptr<Input> &&in,
                                          int recursion = 0,
                                          int max_recursions = 4) {
        std::unique_ptr<Input> ret;

        if (recursion < max_recursions) {
            in->get("include", [&](Input &array_node) {
                array_node.get_all([&](Input &node) {
                    Path relative_path;
                    node.get("path", relative_path);
                    if (relative_path.empty()) {
                        ret = std::move(in);
                        utopia::err() << "relative_path is empty\n";
                    } else {
                        std::string insert_node;
                        node.get("node", insert_node);

                        auto temp = std::make_unique<InputParameters>();
                        Path absolute_path = script_dir / relative_path;
                        Path sub_script_dir = absolute_path.parent();

                        temp->add_root(std::move(in));

                        auto included_file =
                            include(sub_script_dir, open_istream(absolute_path), recursion + 1, max_recursions);

                        if (insert_node.empty()) {
                            temp->add_root(std::move(included_file));
                        } else {
                            temp->add_node(insert_node, std::move(included_file));
                        }

                        ret = std::move(temp);
                    }
                });
            });
        }

        if (!ret) {
            ret = std::move(in);
        }

        return ret;
    }

    void InputParameters::init(const int argc, char *argv[], const bool verbose) {
        std::string key;
        std::string value;
        std::size_t key_len = 0;
        Path script_dir = ".";

        for (int i = 1; i < argc; i++) {
            key = argv[i];
            key_len = key.size();

            if (key.compare(0, 5, "@file") == 0 && i + 1 < argc) {
                auto pos = key.find_first_of(".", 5, 1);

                script_dir = Path(argv[i + 1]).parent();

                if (pos == std::string::npos) {
                    add_root(include(script_dir, open_istream(argv[i + 1])));

                    ++i;  // INCREMENT
                } else {
                    std::string node_key = key.substr(pos + 1, key.size() - pos - 1);
                    auto istr = include(script_dir, open_istream(argv[i + 1]));

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

        if (verbose && mpi_world_rank() == 0) {
            describe(std::cout);
        }
    }

}  // namespace utopia
