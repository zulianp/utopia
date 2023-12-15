#include "utopia_fe_config.hpp"

#ifdef UTOPIA_ENABLE_SFEM

#include "utopia_Input.hpp"
#include "utopia_Main.hpp"
#include "utopia_make_unique.hpp"

#include "utopia_SDFResample_impl.hpp"
#include "utopia_sfem_stk_ExtractSurface.hpp"
#include "utopia_sfem_stk_SDFObstacle.hpp"

#include "utopia_stk.hpp"
#include "utopia_stk_FunctionSpace.hpp"
#include "utopia_stk_intrepid2.hpp"

#include <memory>

namespace utopia {
    template <class FunctionSpace>
    class SDFResampleApp : public Configurable {
    public:
        using Vector = typename Traits<FunctionSpace>::Vector;

        void read(Input &in) override {
            space = std::make_shared<FunctionSpace>();
            resample = utopia::make_unique<SDFResample<FunctionSpace>>();

            output_path = "out.e";

            in.get("output_path", output_path);
            in.require("space", *space);
            in.require("resample", *resample);
        }

        bool is_valid() { return space && resample; }

        void run() {
            //
            Vector sdf;
            resample->apply(*space, sdf);
            space->write(output_path, sdf);
        }

        std::shared_ptr<FunctionSpace> space;
        std::unique_ptr<SDFResample<FunctionSpace>> resample;
        Path output_path;
    };
}  // namespace utopia

void stk_sdf_resample(utopia::Input &in) {
    utopia::SDFResampleApp<utopia::stk::FunctionSpace> app;
    app.read(in);
    if (app.is_valid()) {
        app.run();
    } else {
        utopia::err() << "[Error] invalid set-up\n";
    }
}

UTOPIA_REGISTER_APP(stk_sdf_resample);

#endif
