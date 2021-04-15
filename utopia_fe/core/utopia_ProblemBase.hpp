#ifndef UTOPIA_FE_MODEL_BASE_HPP
#define UTOPIA_FE_MODEL_BASE_HPP

#include "utopia_Function.hpp"
#include "utopia_Traits.hpp"

#include "utopia_fe_Environment.hpp"

namespace utopia {

    template <class FunctionSpace>
    class ProblemBase : public Configurable, public Describable {
    public:
        using Matrix = typename Traits<FunctionSpace>::Matrix;
        using Vector = typename Traits<FunctionSpace>::Vector;
        using Scalar = typename Traits<FunctionSpace>::Scalar;
        using Communicator = typename Traits<FunctionSpace>::Communicator;

        virtual bool is_linear() const { return false; }
        virtual bool is_time_dependent() const { return false; }
        virtual bool is_moving_domain() const { return false; }

        virtual bool is_trivial() const { return is_linear() && !is_moving_domain(); }

        virtual bool solve() = 0;
        virtual bool assemble() = 0;

        virtual bool update() { return false; }

        virtual bool integrate() {
            assert(false && "IMPLEMENT ME");
            return false;
        }

        virtual void condensed_system_built() {}

        // Does nothing here
        virtual void increment_time() {}

        virtual bool export_result() const = 0;
        virtual bool complete() const { return true; }

        virtual bool apply_constraints() {
            this->space()->apply_constraints(*this->jacobian(), *this->fun());

            return true;
        }

        void export_tensors() {
            if (export_tensors_) {
                auto ops = operators();

                for (auto& o : ops) {
                    write("load_" + o->name() + ".m", *o);
                }

                write("load_" + fun_->name() + ".m", *fun_);
            }
        }

        virtual bool init() {
            space()->create_matrix(*this->jacobian());
            space()->create_vector(*this->solution());
            space()->create_vector(*this->fun());

            rename(name() + "_jacobian", *this->jacobian());
            rename(name() + "_solution", *this->solution());
            rename(name() + "_fun", *this->fun());
            return true;
        }

        virtual std::vector<std::shared_ptr<Matrix>> operators() {
            std::vector<std::shared_ptr<Matrix>> ret{jacobian()};
            return ret;
        }

        virtual std::vector<std::shared_ptr<Vector>> right_hand_sides() {
            std::vector<std::shared_ptr<Vector>> ret{fun()};
            return ret;
        }

        std::shared_ptr<Matrix> jacobian() const { return jacobian_; }
        std::shared_ptr<Vector> fun() const { return fun_; }
        std::shared_ptr<Vector> solution() const { return solution_; }

        inline void set_space(const std::shared_ptr<FunctionSpace>& space) { space_ = space; }

        std::shared_ptr<FunctionSpace> space() const { return space_; }
        const std::string& name() const { return name_; }

        void read(Input& in) override {
            in.get("name", name_);
            in.get("output_dir", output_dir_);
            in.get("export_tensors", export_tensors_);
        }

        const Path& output_dir() const { return output_dir_; }

        ProblemBase(const std::shared_ptr<FunctionSpace>& space)
            : space_(space),
              jacobian_(std::make_shared<Matrix>()),
              fun_(std::make_shared<Vector>()),
              solution_(std::make_shared<Vector>()) {}

        virtual void set_environment(const std::shared_ptr<Environment<FunctionSpace>>& env) = 0;

        void describe(std::ostream& os) const override {
            os << "name:\t" << name_ << '\n';
            if (space_) {
                os << "space:\t" << space_->name() << '\n';
            }

            os << "output_dir:\t" << output_dir_ << '\n';
            os << "export_tensors:\t" << export_tensors_ << '\n';
        }

    private:
        std::string name_{"no_name"};
        std::shared_ptr<FunctionSpace> space_;
        std::shared_ptr<Matrix> jacobian_;
        std::shared_ptr<Vector> fun_;
        std::shared_ptr<Vector> solution_;

        Path output_dir_{"./"};
        bool export_tensors_{false};
    };

    template <typename T, typename Integer = int>
    class Time : public Configurable {
    public:
        inline Integer step() const { return step_; }
        Integer& step() { return step_; }
        inline T get() const { return time_; }
        T& get() { return time_; }

        inline void update(const T& delta_time) {
            time_ += delta_time;
            step_ += 1;
        }

        inline void read(Input& in) override {
            in.get("step", step_);
            in.get("time", time_);
        }

    private:
        Integer step_{0};
        T time_{0};
    };

    template <class FunctionSpace>
    class TimeDependentProblem : public ProblemBase<FunctionSpace> {
    public:
        using Super = utopia::ProblemBase<FunctionSpace>;
        using Matrix = typename Traits<FunctionSpace>::Matrix;
        using Vector = typename Traits<FunctionSpace>::Vector;
        using Scalar = typename Traits<FunctionSpace>::Scalar;
        using Time = utopia::Time<Scalar>;

        TimeDependentProblem(const std::shared_ptr<FunctionSpace>& space) : Super(space) {}

        void increment_time() override { current_time_.update(delta_time_); }

        bool integrate() override {
            increment_time();
            return this->update();
        }

        inline bool complete() const override {
            return !(current_time_.step() < end_time_.step() && current_time_.get() < end_time_.get());
        }

        void describe(std::ostream& os) const override {
            Super::describe(os);
            os << "delta_time:\t" << delta_time_ << '\n';
            os << "current_time:\t" << current_time_.get() << '\n';
            os << "end_time:\t" << end_time_.get() << '\n';
            os << "integrate_all_before_output:\t" << integrate_all_before_output_ << '\n';
        }

        void read(Input& in) override {
            Super::read(in);

            in.get("time", [this](Input& node) {
                node.get("delta", delta_time_);
                node.get("start", current_time_.get());
                node.get("steps", end_time_.step());

                assert(delta_time_ != 0.0);
                assert(end_time_.step() >= 1);

                end_time_.get() = current_time_.get() + delta_time_ * end_time_.step();
            });

            in.get("integrate_all_before_output", integrate_all_before_output_);
        }

        inline bool is_time_dependent() const override { return true; }
        inline void set_integrate_all_before_output(bool solve_all) { integrate_all_before_output_ = solve_all; }
        inline bool integrate_all_before_output() const { return integrate_all_before_output_; }

        inline Scalar delta_time() const { return delta_time_; }
        inline const Time& current_time() const { return current_time_; }

    private:
        Scalar delta_time_;
        Time current_time_, end_time_;
        bool integrate_all_before_output_{false};
    };
}  // namespace utopia

#endif  // UTOPIA_FE_MODEL_BASE_HPP
