#include "utopia_mars_Library.hpp"

#include "utopia_fe_base.hpp"
#include "utopia_make_unique.hpp"

// #include "mars.hpp"
#include "mars_env.hpp"

namespace utopia {

    class MarsLibrary::Impl {
    public:
        Impl(int argc, char *argv[]) : env(argc, argv, Traits<UVector>::Communicator::get_default().raw_comm()) {}
        ::mars::Env env;
    };

    MarsLibrary::MarsLibrary() {}
    MarsLibrary::~MarsLibrary() {}
    void MarsLibrary::init(int argc, char *argv[]) { impl_ = utopia::make_unique<Impl>(argc, argv); }

    int MarsLibrary::finalize() {
        impl_ = nullptr;
        return 0;
    }
    std::string MarsLibrary::name() const { return "Mars"; }
    std::string MarsLibrary::version_string() const { return 0; }

}  // namespace utopia
