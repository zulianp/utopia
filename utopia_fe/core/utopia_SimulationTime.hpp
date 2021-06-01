#ifndef UTOPIA_TIME_HPP
#define UTOPIA_TIME_HPP

namespace utopia {

    template <typename T, typename Integer = int>
    class SimulationTime : public Configurable {
    public:
        inline Integer step() const { return step_; }
        inline Integer& step() { return step_; }
        inline T get() const { return time_; }
        inline T& get() { return time_; }

        inline T& delta() { return delta_; }
        inline const T& delta() const { return delta_; }

        inline bool finished() const { return time_ >= end_; }

        inline void update(const T& delta_time) {
            delta_ = delta_time;
            time_ += delta_time;
            step_ += 1;
        }

        inline void update() {
            time_ += delta_;
            step_ += 1;
        }

        void read(Input& in) override {
            in.get("step", step_);
            in.get("delta", delta_);
            in.get("start", start_);

            time_ = start_;

            in.get("time", time_);

            Integer n_steps = 1;
            in.get("n_steps", n_steps);
            end_ = start_ + n_steps * delta_;

            in.get("end", end_);
        }

    private:
        Integer step_{0};
        T time_{0};

        T delta_{0.1};
        T start_{0};
        T end_{0};
    };

}  // namespace utopia

#endif  // UTOPIA_TIME_HPP