#include "utopia_ui.hpp"

#include "utopia_Base.hpp"

#include "utopia_Convertible.hpp"
#include "utopia_InputParameters.hpp"
#include "utopia_JSONInput.hpp"
#include "utopia_Path.hpp"
#include "utopia_XMLInput.hpp"

#ifdef UTOPIA_WITH_YAML_CPP
#include "utopia_YAMLInput.hpp"
#endif  // UTOPIA_WITH_YAML_CPP

namespace utopia {
    template class Convertible<double>;
    template class Convertible<float>;
    template class Convertible<long>;
    template class Convertible<int>;
    template class Convertible<bool>;
    template class Convertible<std::string>;

    std::unique_ptr<Input> open_istream(const Path &path) {
        if (path.extension() == "xml") {
            // auto ret = make_unique<XMLInput>();
            // ret->open(path);
            // return ret;

            // portable version when compuling with nvcc
            auto ret = new XMLInput();
            auto ret_ptr = std::unique_ptr<Input>(ret);
            ret->open(path);
            return ret_ptr;

        } else
#ifdef UTOPIA_WITH_JSON
            if (path.extension() == "json") {

            auto ret = new JSONInput();
            auto ret_ptr = std::unique_ptr<Input>(ret);
            ret->open(path);
            return ret_ptr;
        } else
#endif  // UTOPIA_WITH_JSON
#ifdef UTOPIA_WITH_YAML_CPP
            if (path.extension() == "yaml" || path.extension() == "yml") {

            auto ret = new YAMLInput();
            auto ret_ptr = std::unique_ptr<Input>(ret);
            ret->open(path);
            return ret_ptr;
        } else
#endif  // UTOPIA_WITH_YAML_CPP
        {
            std::cerr << "[Error] format not supported" << std::endl;
            return nullptr;
        }
    }
}  // namespace utopia
