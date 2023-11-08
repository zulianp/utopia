#ifndef UTOPIA_DIRICHLET_BOUNDARY_HPP
#define UTOPIA_DIRICHLET_BOUNDARY_HPP

#include "utopia_Describable.hpp"
#include "utopia_Input.hpp"
#include "utopia_SimulationTime.hpp"
#include "utopia_SymbolicFunction.hpp"
#include "utopia_make_unique.hpp"

#include "utopia_IOStream.hpp"
#include "utopia_Instance.hpp"
#include "utopia_MPI.hpp"

#include <string>
#include <vector>

namespace utopia {

    template <class Traits>
    class DirichletBoundary : public Configurable, public Describable {
    public:
        using Scalar = typename Traits::Scalar;
        using SizeType = typename Traits::SizeType;
        using Vector = typename Traits::Vector;

        class Condition : public Configurable, public Describable {
        public:
            virtual ~Condition() = default;
            virtual void update(const SimulationTime<double> &) {}
            // virtual double value() const = 0;

            virtual bool is_uniform() const = 0;

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
            static constexpr const char *class_name() { return "UniformCondition"; }

            double value_;

            // inline double value() const override { return value_; }
            virtual double value() const { return value_; }
            bool is_uniform() const override { return true; }

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

#ifdef UTOPIA_WITH_TINY_EXPR
        class VaryingCondition : public Condition {
        public:
            using Super = Condition;
            static constexpr const char *class_name() { return "VaryingCondition"; }

            std::unique_ptr<utopia::SymbolicFunction> expr_;
            double t_{0};
            bool verbose_{false};

            void update(const SimulationTime<double> &t) override {
                t_ = t.get();

                if (verbose_ && mpi_world_rank() == 0) {
                    utopia::out() << "VaryingCondition::update(" << t_ << ")\n";
                }
            }

            inline Scalar eval(const Scalar x, const Scalar y, const Scalar z) const {
                return expr_->eval(x, y, z, t_);
            }

            bool is_uniform() const override { return false; }

            VaryingCondition() : Super(), expr_(utopia::make_unique<SymbolicFunction>("0")) {}
            VaryingCondition(std::string name, const int component)
                : Super(std::move(name), component), expr_(utopia::make_unique<SymbolicFunction>("0")) {}

            void read(Input &in) override {
                Super::read(in);

                std::string expr;
                in.get("value", expr);
                if (!expr.empty()) {
                    *expr_ = utopia::symbolic(expr);

                    if (!expr_->valid()) {
                        Utopia::Abort("VaryingCondition: invalid expression: " + expr);
                    }
                }

                in.get("verbose", verbose_);
            }

            void describe(std::ostream &os) const override {
                Super::describe(os);
                os << "value:\t" << expr_->to_string() << '\n';
            }
        };

        class TimeDependentCondition : public UniformCondition {
        public:
            using Super = UniformCondition;
            static constexpr const char *class_name() { return "TimeDependentCondition"; }

            std::unique_ptr<utopia::SymbolicFunction> expr_;
            double t_{0};

            void update(const SimulationTime<double> &t) override {
                t_ = t.get();
                this->value_ = expr_->eval(0, 0, 0, t_);
            }

            TimeDependentCondition() : Super(), expr_(utopia::make_unique<SymbolicFunction>("0")) {}

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
#endif  // UTOPIA_WITH_TINY_EXPR

        class OverwriteCondition : public Condition {
        public:
            using Super = Condition;

            static constexpr const char *class_name() { return "OverwriteCondition"; }

            void read(Input &in) override { Super::read(in); }

            void describe(std::ostream &os) const override { Super::describe(os); }

            bool is_uniform() const override { return false; }

            OverwriteCondition(const std::shared_ptr<Vector> &vector) : vector(vector) {}
            std::shared_ptr<Vector> vector;
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
#ifdef UTOPIA_WITH_TINY_EXPR
                if (type == TimeDependentCondition::class_name()) {
                    c = std::make_shared<TimeDependentCondition>();
                } else if (type == VaryingCondition::class_name()) {
                    // printf("VaryingCondition\n");
                    c = std::make_shared<VaryingCondition>();
                } else
#endif  // UTOPIA_WITH_TINY_EXPR
                    if (type.empty() || type == UniformCondition::class_name()) {
                        c = std::make_shared<UniformCondition>();
                    } else if (type == OverwriteCondition::class_name()) {
                        c = std::make_shared<OverwriteCondition>(nullptr);
                    } else {
                        return;
                    }

                c->read(in);
                this->add(c);
            });
        }

        void read_with_state(const std::shared_ptr<Vector> &v, Input &in) {
            in.get_all([this, &v](Input &in) {
                std::string type;
                in.get("type", type);

                std::shared_ptr<Condition> c;
#ifdef UTOPIA_WITH_TINY_EXPR
                if (type == TimeDependentCondition::class_name()) {
                    c = std::make_shared<TimeDependentCondition>();
                } else if (type == VaryingCondition::class_name()) {
                    // printf("VaryingCondition\n");
                    c = std::make_shared<VaryingCondition>();
                } else
#endif  // UTOPIA_WITH_TINY_EXPR
                    if (type.empty() || type == UniformCondition::class_name()) {
                        c = std::make_shared<UniformCondition>();
                    } else if (type == OverwriteCondition::class_name()) {
                        c = std::make_shared<OverwriteCondition>(v);
                    } else {
                        return;
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

        using ConstIterator = typename std::vector<std::shared_ptr<Condition>>::const_iterator;

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
