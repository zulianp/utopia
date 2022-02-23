#include "utopia_Main.hpp"

#include "utopia.hpp"

#include "utopia_SimpleNewton.hpp"

#include "utopia_ImplicitEulerIntegrator.hpp"
#include "utopia_NewmarkIntegrator.hpp"

#include "utopia_stk.hpp"
#include "utopia_stk_intrepid2_OmniAssembler.hpp"

namespace utopia {

    template class NewmarkIntegrator<utopia::stk::FunctionSpace>;
    template class ImplicitEulerIntegrator<utopia::stk::FunctionSpace>;

    template <class FunctionSpace>
    class ResampleTimeSeriesApp : public Configurable {
    public:
        using Matrix_t = typename Traits<FunctionSpace>::Matrix;
        using Vector_t = typename Traits<FunctionSpace>::Vector;
        using Scalar_t = typename Traits<FunctionSpace>::Scalar;

        using IO_t = utopia::IO<FunctionSpace>;

        void read(Input &in) override {
            space_ = std::make_shared<FunctionSpace>();
            input_ = std::make_shared<IO_t>(*space_);
            input_->import_all_field_data(true);
            input_->enable_interpolation_mode();

            in.get("space", [&](Input &node) { valid_ = input_->open_input(node); });

            if (space_->empty()) {
                valid_ = false;
            }

            if (!valid_) {
                return;
            }

            output_ = std::make_shared<IO_t>(*space_);

            in.get("n_time_steps", n_time_steps_);
            in.get("delta_time", delta_time_);
            in.get("output", output_path);

            output_->set_output_path(output_path);
            output_->open_output();
        }

        bool valid() const { return valid_; }

        void run() {
            bool ok;
            for (int t = 0; t < n_time_steps_; ++t) {
                ok = input_->load_time_step(t * delta_time_);
                assert(ok);

                if (!ok) {
                    Utopia::Abort("Failed to load time step");
                }

                ok = output_->write(t, t * delta_time_);
                assert(ok);

                if (!ok) {
                    Utopia::Abort("Failed to write time step");
                }
            }
        }

        std::shared_ptr<FunctionSpace> space_;
        std::shared_ptr<IO_t> input_;
        std::shared_ptr<IO_t> output_;
        int n_time_steps_{2};
        Scalar_t delta_time_{0.1};
        Path output_path{"out.e"};
        bool valid_{false};
    };

}  // namespace utopia

void stk_resample(utopia::Input &in) {
    utopia::ResampleTimeSeriesApp<utopia::stk::FunctionSpace> app;
    app.read(in);

    if (app.valid()) {
        app.run();
    } else {
        utopia::err() << "stk_resample: invalid app setup\n";
    }
}

UTOPIA_REGISTER_APP(stk_resample);
