#ifndef UTOPIA_GENERIC_ACTION_REGISTRY_H
#define UTOPIA_GENERIC_ACTION_REGISTRY_H

#include <map>
#include <string>
#include <functional>
#include <iostream>

namespace utopia {

    template<class... Args>
    class NaryActionRegistry {
    public:
        using Count = long;
        typedef void (*ExecuteAction)(Args &&...);

        char add_action(const std::string &action_name, ExecuteAction apply_test);
        int apply(const std::string &action_name, Args &&...args);

        void describe(std::ostream &os = std::cout) const;

        inline bool verbose() const { return verbose_; }
        inline void verbose(const bool val) { verbose_ = val; }
        inline void action_applied() { ++n_action_applied_; }
        inline Count n_action_applied() { return n_action_applied_; }
        inline bool empty() const { return actions_.empty(); }
        inline std::size_t size() const { return actions_.size(); }


        inline static NaryActionRegistry &instance()
        {
            static NaryActionRegistry instance_;
            return instance_;
        }

        inline void set_type(const std::string &type)
        {
            type_ = type;
        }

        inline const std::string &type() const
        {
            return type_;
        }

        NaryActionRegistry();
        virtual ~NaryActionRegistry();

    private:
        std::map<std::string, ExecuteAction> actions_;
        bool verbose_;
        int error_code_ = 0;
        Count n_action_applied_;
        std::string type_;
        int rank_;
    };

}

#endif //UTOPIA_GENERIC_ACTION_REGISTRY_H
