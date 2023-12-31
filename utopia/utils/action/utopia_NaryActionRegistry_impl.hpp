#ifndef UTOPIA_NARY_ACTION_REGISTRY_IMPL_HPP
#define UTOPIA_NARY_ACTION_REGISTRY_IMPL_HPP

#include "utopia_ActionRegistry.hpp"
#include "utopia_Chrono.hpp"
#include "utopia_MPI.hpp"

#include <iostream>

namespace utopia {

    template <class... Args>
    NaryActionRegistry<Args...>::NaryActionRegistry() : type_("action") {
        // rank_ = utopia::mpi_world_rank();
    }

    template <class... Args>
    NaryActionRegistry<Args...>::~NaryActionRegistry() {
        if (rank_ == 0 && n_action_applied_ > 0) {
            std::cout << "[Status] " << type() << " executed " << n_action_applied_ << std::endl;
        }
    }

    template <class... Args>
    char NaryActionRegistry<Args...>::add_action(const std::string &action_name,
                                                 NaryActionRegistry::ExecuteAction action) {
        actions_[action_name] = action;
        return 0;
    }

    template <class... Args>
    void NaryActionRegistry<Args...>::describe(std::ostream &os) const {
        for (const auto &u : actions_) {
            os << "\t" << u.first << "\n";
        }

        os << std::flush;
    }

    template <class... Args>
    int NaryActionRegistry<Args...>::apply_with_backend(const std::string &action_name,
                                                        const std::string &backend,
                                                        Args &&... args) {
        rank_ = mpi_world_rank();

        std::string search_string;
        if (!backend.empty()) {
            search_string = backend + "_";
        }

        search_string += action_name;

        auto it = actions_.find(search_string);
        if (it == actions_.end()) {
            // if(rank_ == 0) {
            //     std::cerr << "[Error] no " << type() << " with name " << action_name << std::endl;
            // }

            return -1;
        }

        try {
            if (rank_ == 0 && verbose()) {
                std::cout << "--------------------------------------------------------\n";
                std::cout << "begin:\t[" << action_name << "] " << type() << std::endl;
                std::cout << "--------------------------------------------------------\n";
            }

            it->second(std::forward<Args>(args)...);

            if (rank_ == 0 && verbose()) {
                std::cout << "--------------------------------------------------------\n";
                std::cout << "end:\t[" << action_name << "] " << type() << std::endl;
                std::cout << "--------------------------------------------------------\n";
            }

            action_applied();

        } catch (std::exception &ex) {
            std::cerr << "[Failure] in " << action_name << " " << ex.what() << std::endl;
            error_code_ = 1;
        }

        return error_code_;
    }

}  // namespace utopia

#endif
