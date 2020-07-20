#include "utopia_Reporter.hpp"

namespace utopia {

    void Reporter::read(Input &in) {
        bool mute = false;
        in.get("mute", mute);

        if (mute) {
            cout_ = utopia::make_unique<NullAppOutputStream>();
            cerr_ = utopia::make_unique<NullAppOutputStream>();
            dev_ = utopia::make_unique<NullAppOutputStream>();
        } else {
            bool dev_mode = true;
            in.get("dev_mode", dev_mode);

            if (!dev_mode) {
                dev_ = utopia::make_unique<NullAppOutputStream>();
            }
        }
    }
}  // namespace utopia