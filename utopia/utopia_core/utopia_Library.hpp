#ifndef UTOPIA_LIBRARY_HPP
#define UTOPIA_LIBRARY_HPP

#include <string>

namespace utopia {

    class Library {
    public:
        virtual ~Library() {}
        virtual void init(int argc, char *argv[]) = 0;
        virtual int finalize() = 0;
        virtual std::string name() const = 0;
        virtual std::string version_string() const = 0;
    };

}  // namespace utopia

#endif  // UTOPIA_LIBRARY_HPP
