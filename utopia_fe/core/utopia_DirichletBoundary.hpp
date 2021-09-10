#ifndef UTOPIA_DIRICHLET_BOUNDARY_HPP
#define UTOPIA_DIRICHLET_BOUNDARY_HPP

#include "utopia_Describable.hpp"
#include "utopia_Input.hpp"
#include "utopia_SimulationTime.hpp"
#include "utopia_SymbolicFunction.hpp"
#include "utopia_make_unique.hpp"

#include <string>
#include <vector>

namespace utopia {

    class DirichletBoundary : public Configurable, public Describable {
    public:
        class Condition : public Configurable, public Describable {
        public:
            virtual ~Condition() = default;
            virtual void update(const SimulationTime<double> &) {}
            virtual double value() const = 0;

            Condition() = default;
            Condition(std::string name, const int component) : name(std::move(name)), component(component) {}

            void read(Input &in) override {
                in.get("name", name);
                in.get("var", component);
                in.get("side", side);

                if (name.empty() && side != -1) {
                    name = "surface_" + std::to_string(side);
                }
            }

            void describe(std::ostream &os) const override {
                os << "name:\t" << name << '\n';
                os << "var:\t" << component << '\n';
                os << "side:\t" << side << '\n';
            }

            std::string name;
            int side{-1};
            int component{0};
        };

        class UniformCondition : public Condition {
        public:
            using Super = Condition;

            double value_;

            inline double value() const override { return value_; }

            UniformCondition() = default;
            UniformCondition(std::string name, double value, const int component)
                : Condition(std::move(name), component), value_(value) {}

            void read(Input &in) override {
                Super::read(in);
                in.get("value", value_);
            }

            void describe(std::ostream &os) const override {
                Super::describe(os);
                os << "value:\t" << value_ << '\n';
            }
        };

        class TimeDependentCondition : public Condition {
        public:
            using Super = Condition;
            static constexpr const char *name() { return "TimeDependentCondition"; }

            std::unique_ptr<utopia::SymbolicFunction> expr_;
            double t_{0};

            inline double value() const override { return expr_->eval(0, 0, 0, t_); }

            TimeDependentCondition() : Condition(), expr_(utopia::make_unique<SymbolicFunction>("0")) {}

            TimeDependentCondition(std::string name, std::string expr, const int component)
                : Condition(std::move(name), component),
                  expr_(utopia::make_unique<SymbolicFunction>(std::move(expr))) {}

            void read(Input &in) override {
                Super::read(in);

                std::string expr;
                in.get("value", expr);
                if (!expr.empty()) {
                    *expr_ = utopia::symbolic(expr);
                }
            }

            void describe(std::ostream &os) const override {
                Super::describe(os);
                os << "value:\t" << expr_->to_string() << '\n';
            }
        };

        template <class Mapper>
        void convert_user_space_names(const Mapper &mapper) {
            for (auto &c_ptr : conditions) {
                assert(c_ptr);

                auto &c = *c_ptr;
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
                std::string type;
                in.get("type", type);

                std::shared_ptr<Condition> c;
                if (type == TimeDependentCondition::name()) {
                    c = std::make_shared<TimeDependentCondition>();
                } else {
                    c = std::make_shared<UniformCondition>();
                }

                c->read(in);
                this->add(c);
            });
        }

        void describe(std::ostream &os) const override {
            os << "DirichletBoundary:\n";
            for (auto &c : conditions) {
                if (c) c->describe(os);
            }
        }

        using ConstIterator = std::vector<std::shared_ptr<Condition>>::const_iterator;

        inline ConstIterator begin() const { return conditions.cbegin(); }
        inline ConstIterator end() const { return conditions.cend(); }

        void add(const UniformCondition &cond) { conditions.push_back(std::make_shared<UniformCondition>(cond)); }
        void add(const std::shared_ptr<Condition> &cond) { conditions.push_back(cond); }

        void update(const SimulationTime<double> &time) {
            for (auto &c_ptr : conditions) {
                c_ptr->update(time);
            }
        }

    private:
        std::vector<std::shared_ptr<Condition>> conditions;
    };
}  // namespace utopia

#endif  // UTOPIA_DIRICHLET_BOUNDARY_HPP
