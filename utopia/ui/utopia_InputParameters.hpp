#ifndef UTOPIA_INPUT_PARAMETERS_HPP
#define UTOPIA_INPUT_PARAMETERS_HPP

#include "utopia_Convertible.hpp"
#include "utopia_Input.hpp"
#include "utopia_make_unique.hpp"

#include <cassert>
#include <map>
#include <memory>
#include <string>
#include <utility>

namespace utopia {

    class InputParameters final : public Input {
    public:
        using Input::get;

        void init(const int argc, char *argv[], const bool verbose = false);

        inline bool empty() const { return nodes_.empty() && values_.empty(); }

        inline SizeType size() const override { return nodes_.size() + values_.size(); }

        inline void get(const std::string &key, bool &val) override { aux_get(key, val); }

        inline void get(const std::string &key, double &val) override { aux_get(key, val); }

        inline void get(const std::string &key, int &val) override { aux_get(key, val); }

        inline void get(const std::string &key, std::string &val) override { aux_get(key, val); }

        inline void get(const std::string &key, long &val) override { aux_get(key, val); }

        inline void get(const std::string &key, unsigned long &val) override { aux_get(key, val); }

        inline void get(const std::string &key, long long int &val) override { aux_get(key, val); }

        inline void get(const std::string &key, Configurable &val) override {
            auto node_ptr = node(key);

            if (node_ptr) {
                val.read(*node_ptr);
            } else {
                // if (aux_root_) {
                //     aux_root_->get(key, val);
                // }

                for (auto ar = aux_roots_.rbegin(); ar != aux_roots_.rend(); ++ar) {
                    (*ar)->get(key, val);
                }
            }
        }

        inline void get(const std::string &key, std::function<void(Input &)> lambda) override {
            auto node_ptr = node(key);

            if (node_ptr) {
                lambda(*node_ptr);
            } else {
                // if (aux_root_) {
                //     aux_root_->get(key, lambda);
                // }

                for (auto ar = aux_roots_.rbegin(); ar != aux_roots_.rend(); ++ar) {
                    (*ar)->get(key, lambda);
                }
            }
        }

        void get_all(std::function<void(Input &)> lambda) override {
            for (const auto &n : nodes_) {
                lambda(*n.second);
            }

            for (auto ar = aux_roots_.rbegin(); ar != aux_roots_.rend(); ++ar) {
                (*ar)->get_all(lambda);
            }

            // if (aux_root_) {
            //     aux_root_->get_all(lambda);
            // }
        }

        void get(std::vector<std::shared_ptr<IConvertible>> &values) override {
            if (values_.empty()) return;

            values.reserve(values_.size());

            for (const auto &v : values_) {
                values.push_back(std::shared_ptr<IConvertible>(v.second->clone()));
            }
        }

        inline bool good() const override { return true; }

        std::shared_ptr<Input> node(const std::string &key) const {
            auto it = nodes_.find(key);
            if (it != nodes_.end()) {
                return it->second;
            }

            // IS this assert necessary???
            // assert(false);
            return nullptr;
        }

        inline void set(const std::string &key, const bool &val) { aux_set(key, val); }

        inline void set(const std::string &key, const double &val) { aux_set(key, val); }

        inline void set(const std::string &key, const int &val) { aux_set(key, val); }

        inline void set(const std::string &key, const long &val) { aux_set(key, val); }

        inline void set(const std::string &key, const long long &val) { aux_set(key, val); }

        inline void set(const std::string &key, const char *val) { aux_set(key, std::string(val)); }

        inline void set(const std::string &key, const std::string &val) { aux_set(key, val); }

        inline void set(const std::string &key, const std::shared_ptr<Input> &in) { nodes_[key] = in; }

        void describe(std::ostream &os) const { aux_describe(os, 0); }

        inline InputParameters() {}
        InputParameters(InputParameters &&other)
            : Input(std::move(other)),
              values_(std::move(other.values_)),
              nodes_(std::move(other.nodes_)),
              aux_roots_(std::move(other.aux_roots_)) {}

        inline void add_root(std::unique_ptr<Input> &&root) { aux_roots_.push_back(std::move(root)); }
        inline void add_node(const std::string &name, std::shared_ptr<Input> &&node) { nodes_[name] = std::move(node); }

        bool key_exists(const std::string &key) const override;

    private:
        std::map<std::string, std::unique_ptr<IConvertible>> values_;
        std::map<std::string, std::shared_ptr<Input>> nodes_;
        std::vector<std::unique_ptr<Input>> aux_roots_;

        template <typename Out>
        void aux_get(const std::string &key, Out &out) const {
            auto it = values_.find(key);

            if (it != values_.end()) {
                it->second->get(out);
            } else
            // if (aux_root_)
            {
                for (auto ar = aux_roots_.rbegin(); ar != aux_roots_.rend(); ++ar) {
                    (*ar)->get(key, out);
                }
            }
        }

        template <typename In>
        void aux_set(const std::string &key, const In &in) {
            values_[key] = utopia::make_unique<Convertible<In>>(in);
        }

        void aux_describe(std::ostream &os, const int level) const {
            std::string indent(""), str;
            indent.resize(2 * level, ' ');

            for (const auto &kv : values_) {
                kv.second->get(str);
                os << indent << kv.first << " : " << str << "\n";
            }

            for (const auto &n : nodes_) {
                os << indent << n.first << " : {node}\n";
            }
        }
    };

    template <class Arg>
    inline std::pair<std::string, Arg> param(const std::string &key, Arg &&value) {
        return {key, std::forward<Arg>(value)};
    }

    inline std::pair<std::string, std::string> param(const std::string &key, const char *value) { return {key, value}; }

    template <class First>
    inline void param_set(InputParameters &params, std::pair<std::string, First> &&p) {
        params.set(std::forward<std::string>(p.first), std::forward<First>(p.second));
    }

    inline void param_set(InputParameters &params, std::pair<std::string, InputParameters> &&p) {
        params.set(std::move(p.first), std::make_shared<InputParameters>(std::move(p.second)));
    }

    template <class First>
    inline void param_append(InputParameters &params, std::pair<std::string, First> &&first) {
        param_set(params, std::forward<std::pair<std::string, First>>(first));
    }

    template <class First, class... Args>
    inline void param_append(InputParameters &params,
                             std::pair<std::string, First> &&first,
                             std::pair<std::string, Args> &&... list) {
        params.set(std::forward<std::string>(first.first), std::forward<First>(first.second));
        param_append(params, std::forward<std::pair<std::string, Args>>(list)...);
    }

    template <class... Args>
    inline InputParameters param_list(std::pair<std::string, Args> &&... list) {
        InputParameters ret;
        param_append(ret, std::forward<std::pair<std::string, Args>>(list)...);
        return ret;
    }

    template <class Arg>
    inline InputParameters param_list(std::pair<std::string, Arg> &&p) {
        InputParameters ret;
        param_set(ret, std::forward<std::pair<std::string, Arg>>(p));
        return ret;
    }
}  // namespace utopia

#endif  // UTOPIA_INPUT_PARAMETERS_HPP
