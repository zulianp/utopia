#include "utopia_Reporter.hpp"

namespace utopia {

    void Reporter::read(Input &in) {
        bool mute = false;
        in.get("mute", mute);

        if (mute) {
            cout_ = utopia::make_unique<NullOStream>();
            cerr_ = utopia::make_unique<NullOStream>();
            dev_ = utopia::make_unique<NullOStream>();
        } else {
            bool dev_mode = true;
            in.get("dev_mode", dev_mode);

            if (!dev_mode) {
                dev_ = utopia::make_unique<NullOStream>();
            }
        }
    }

    OStream &out() { return Reporter::instance().cout(); }
    OStream &err() { return Reporter::instance().cerr(); }
    OStream &dev() { return Reporter::instance().dev(); }

}  // namespace utopia
