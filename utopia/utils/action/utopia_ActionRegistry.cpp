#include "utopia_ActionRegistry.hpp"
#include "utopia_MPI.hpp"
#include "utopia_Chrono.hpp"

#include <iostream>

namespace utopia {

    ActionRegistry::ActionRegistry() : type_("action") {
        // rank_ = utopia::mpi_world_rank();
    }

    ActionRegistry::~ActionRegistry()
    {
        if(rank_ == 0 && n_action_applied_ > 0) {
            std::cout << "[Status] " << type() << " executed " << n_action_applied_ << std::endl;
        }
    }

    char ActionRegistry::add_action(const std::string &action_name, ActionRegistry::ExecuteAction action)
    {
        actions_[action_name] = action;
        return 0;
    }

    void ActionRegistry::describe(std::ostream &os) const
    {
        for(const auto &u : actions_) {
            os << "\t" << u.first << "\n";
        }

        os << std::flush;
    }

    void ActionRegistry::apply(const std::string &action_name, ExecuteAction action)
    {
        try {
            if(rank_ == 0 && verbose()) {                                     \
                std::cout << "--------------------------------------------------------\n";       \
                std::cout << "begin:\t[" << action_name << "] " << type() << std::endl;  \
                std::cout << "--------------------------------------------------------\n";       \
            }

            action();

            if(rank_ == 0 && verbose()) {                                     \
                std::cout << "--------------------------------------------------------\n";       \
                std::cout << "end:\t[" << action_name << "] " << type()  << std::endl;  \
                std::cout << "--------------------------------------------------------\n";       \
            }

            action_applied();

        } catch(std::exception &ex) {
                std::cerr << "[Failure] in " << action_name << " " << ex.what() << std::endl;
                error_code_ = 1;
        }
    }

    int ActionRegistry::apply_all()
    {
        for(std::map<std::string, ExecuteAction>::const_iterator it = actions_.begin(); it != actions_.end(); ++it) {
            apply(it->first, it->second);
        }

        return error_code_;

    }

    int ActionRegistry::apply_aux(const std::map<std::string, ExecuteAction> &units, const std::string &action_name)
    {
        auto it = units.find(action_name);
        if(it == units.end()) {
            return -1;
        }
        apply(it->first, it->second);
        return error_code_;
    }

    int ActionRegistry::apply(const std::string &action_name)
    {
        rank_ = mpi_world_rank();

        int ret = apply_aux(actions_, action_name);

        // if(ret == -1 && rank_ == 0) {
        //     std::cerr << "[Error] no " << type() << " with name " << action_name << std::endl;
        // }

        return ret;
    }

}
