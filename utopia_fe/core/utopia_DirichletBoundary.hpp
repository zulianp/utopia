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

#include <filesystem>
#include <fstream>
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
            virtual bool has_time_derivative() const { return false; }

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
            virtual double time_derivative() const { return 0; }
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

#ifdef UTOPIA_ENABLE_TINY_EXPR
        class VaryingCondition : public Condition {
        public:
            using Super = Condition;
            static constexpr const char *class_name() { return "VaryingCondition"; }

            std::unique_ptr<utopia::SymbolicFunction> expr_;
            std::unique_ptr<utopia::SymbolicFunction> expr_time_derivative_;

            double t_{0};
            bool verbose_{false};

            bool has_time_derivative() const override { return expr_time_derivative_ != nullptr; }

            void update(const SimulationTime<double> &t) override {
                t_ = t.get();

                if (verbose_ && mpi_world_rank() == 0) {
                    utopia::out() << "VaryingCondition::update(" << t_ << ")\n";
                }
            }

            inline Scalar eval(const Scalar x, const Scalar y, const Scalar z) const {
                assert(expr_);
                return expr_->eval(x, y, z, t_);
            }

            inline Scalar eval_time_derivative(const Scalar x, const Scalar y, const Scalar z) const {
                assert(expr_time_derivative_);
                return expr_time_derivative_->eval(x, y, z, t_);
            }

            bool is_uniform() const override { return false; }

            VaryingCondition()
                : Super(),
                  expr_(utopia::make_unique<SymbolicFunction>("0")),
                  expr_time_derivative_(utopia::make_unique<SymbolicFunction>("0")) {}
            VaryingCondition(std::string name, const int component)
                : Super(std::move(name), component), expr_(utopia::make_unique<SymbolicFunction>("0")) {}

            void read(Input &in) override {
                Super::read(in);

                in.get("verbose", verbose_);

                {
                    std::string expr;
                    in.get("value", expr);

                    if (verbose_ && mpi_world_rank() == 0) {
                        utopia::out() << expr << "\n";
                    }

                    if (!expr.empty()) {
                        *expr_ = utopia::symbolic(expr);

                        if (!expr_->valid()) {
                            Utopia::Abort("VaryingCondition: invalid expression: " + expr);
                        }
                    }
                }

                {
                    std::string dexpr;
                    in.get("time_derivative", dexpr);

                    if (verbose_ && mpi_world_rank() == 0) {
                        utopia::out() << dexpr << "\n";
                    }

                    if (!dexpr.empty()) {
                        *expr_time_derivative_ = utopia::symbolic(dexpr);

                        if (!expr_time_derivative_->valid()) {
                            Utopia::Abort("VaryingCondition: invalid expression: " + dexpr);
                        }
                    }
                }
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

        class SignalCondition : public TimeDependentCondition {
        public:
            using Super = TimeDependentCondition;
            static constexpr const char *class_name() { return "SignalCondition"; }

            std::vector<float> time_;
            std::vector<float> signal_;
            double t_{0};
            bool verbose_{false};
            double derivative_{0};

            double interpolate(const double t) const {
                auto iter = std::lower_bound(time_.begin(), time_.end(), t_);
                if (iter == time_.end() || iter - 1 < time_.begin()) {
                    if (!mpi_world_rank()) {
                        utopia::err() << "Requested time: " << t_ << "\n";
                        utopia::err() << "Found index: " << std::distance(time_.begin(), iter) << "\n";
                    }

                    Utopia::Abort("SignalCondition! Invalid time query!");
                }

#ifndef NDEBUG
                const ptrdiff_t n = time_.size();
#endif

                const ptrdiff_t i0 = std::distance(time_.begin(), iter - 1);
                const ptrdiff_t i1 = std::distance(time_.begin(), iter);

                assert(i0 >= 0);
                assert(i0 < n);

                assert(i1 >= 0);
                assert(i1 < n);

                double t0 = time_[i0];
                double t1 = time_[i1];

                assert(t0 != t1);

                double w0 = (t1 - t) / (t1 - t0);
                double w1 = (t - t0) / (t1 - t0);

                double f0 = signal_[i0];
                double f1 = signal_[i1];

                double f = f0 * w0 + f1 * w1;

                if (verbose_ && !mpi_world_rank()) {
                    utopia::out() << "SignalCondition::value(" << t_ << ") = " << f << "\n";
                    utopia::out() << "t0 " << t0 << "\n";
                    utopia::out() << "t1 " << t1 << "\n";
                    utopia::out() << "w0 " << w0 << "\n";
                    utopia::out() << "w1 " << w1 << "\n";
                    utopia::out() << "f0 " << f0 << "\n";
                    utopia::out() << "f1 " << f1 << "\n";
                }

                return f;
            }

            void update(const SimulationTime<double> &t) override {
                t_ = t.get();

                const double f = interpolate(t_);
                const double fp1 = interpolate(t_ + t.delta());
                this->value_ = f;
                this->derivative_ = (fp1 - f) / t.delta();
            }

            double time_derivative() const override { return this->derivative_; }

            bool has_time_derivative() const override { return true; }

            SignalCondition() : Super() {}

            SignalCondition(std::string name, const int component) : Super(std::move(name), "1", component) {}

            static void read_file(const std::string &path, std::vector<float> &data) {
                std::ifstream is(path, std::ios::binary);
                if (!is.good()) {
                    if (!mpi_world_rank()) {
                        utopia::err() << "Unable to read file " << path << "\n";
                    }
                    Utopia::Abort("Aborting..");
                }

                std::streamsize size = static_cast<std::streamsize>(std::filesystem::file_size(path));
                data.resize(size / sizeof(float));
                is.read((char *)&data[0], size);
                is.close();
            }

            void read(Input &in) override {
                Super::read(in);

                std::string time_file;
                in.require("time", time_file);

                std::string signal_file;
                in.require("signal", signal_file);

                read_file(time_file, time_);
                read_file(signal_file, signal_);

                if (time_.size() != signal_.size()) {
                    if (!mpi_world_rank()) {
                        utopia::err() << "Time and signal have different sizes! " << time_.size() << " != "
                                      << "\n";
                    }

                    Utopia::Abort("Aborting..");
                }

                in.get("verbose", verbose_);
                if (verbose_ && !mpi_world_rank()) {
                    size_t n = time_.size();

                    utopia::out() << "time\tsignal\n";

                    for (size_t i = 0; i < n; i++) {
                        utopia::out() << time_[i] << "\t" << signal_[i] << "\n";
                    }
                }
            }

            void describe(std::ostream &os) const override {
                Super::describe(os);
                os << "signal with size " << signal_.size() << "\n";
            }
        };
#endif  // UTOPIA_ENABLE_TINY_EXPR

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
#ifdef UTOPIA_ENABLE_TINY_EXPR
                if (type == TimeDependentCondition::class_name()) {
                    c = std::make_shared<TimeDependentCondition>();
                } else if (type == VaryingCondition::class_name()) {
                    // printf("VaryingCondition\n");
                    c = std::make_shared<VaryingCondition>();
                } else if (type == SignalCondition::class_name()) {
                    c = std::make_shared<SignalCondition>();
                } else
#endif  // UTOPIA_ENABLE_TINY_EXPR
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
#ifdef UTOPIA_ENABLE_TINY_EXPR
                if (type == TimeDependentCondition::class_name()) {
                    c = std::make_shared<TimeDependentCondition>();
                } else if (type == VaryingCondition::class_name()) {
                    c = std::make_shared<VaryingCondition>();
                } else if (type == SignalCondition::class_name()) {
                    c = std::make_shared<SignalCondition>();
                } else
#endif  // UTOPIA_ENABLE_TINY_EXPR
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
            os << "DirichletBoundary(" << conditions.size() << "):\n";
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
