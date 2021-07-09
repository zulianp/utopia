#ifndef UTOPIA_DIRICHLET_BOUNDARY_HPP
#define UTOPIA_DIRICHLET_BOUNDARY_HPP

#include "utopia_Describable.hpp"
#include "utopia_Input.hpp"

#include <string>
#include <vector>

namespace utopia {

    class DirichletBoundary : public Configurable, public Describable {
    public:
        class Condition : public Configurable, public Describable {
        public:
            std::string name;
            double value;
            int component{0};
            int side{-1};

            Condition() = default;
            Condition(std::string name, double value, const int component)
                : name(std::move(name)), value(value), component(component) {}

            void read(Input &in) override {
                in.get("name", name);
                in.get("value", value);
                in.get("var", component);
                in.get("side", side);

                if (name.empty() && side != -1) {
                    name = "surface_" + std::to_string(side);
                }
            }

            void describe(std::ostream &os) const override {
                os << "name:\t" << name << '\n';
                os << "value:\t" << value << '\n';
                os << "var:\t" << component << '\n';
                os << "side:\t" << side << '\n';
            }
        };

        template <class Mapper>
        void convert_user_space_names(const Mapper &mapper) {
            for (auto &c : conditions) {
                if (!c.name.empty() && c.side == mapper.invalid_id()) {
                    auto id = mapper.id_from_user_space_name(c.name);
                    if (id != mapper.invalid_id()) {
                        c.side = id;
                        c.name = mapper.name_from_id(id);
                    } else {
                        assert(false && "FIXME");
                    }
                }
            }
        }

        void read(Input &in) override {
            in.get_all([this](Input &in) {
                conditions.emplace_back();
                auto &c = conditions.back();
                c.read(in);
            });
        }

        void describe(std::ostream &os) const override {
            os << "DirichletBoundary:\n";
            for (auto &c : conditions) {
                c.describe(os);
            }
        }

        std::vector<Condition> conditions;
    };
}  // namespace utopia

#endif  // UTOPIA_DIRICHLET_BOUNDARY_HPP
