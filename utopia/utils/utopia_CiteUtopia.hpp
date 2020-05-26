#ifndef UTOPIA_CITE_UTOPIA_HPP
#define UTOPIA_CITE_UTOPIA_HPP

#include "utopia_AuthoredWork.hpp"

namespace utopia {
    // Bibliography

    ///////////////////////////////////////////

    class Utopia2016Git {};

    template <>
    class Cite<Utopia2016Git> {
    public:
        static BibTeX bibtex();
    };

    ///////////////////////////////////////////

    class Kopanicakova2020Recursive {};

    template <>
    class Cite<Kopanicakova2020Recursive> {
    public:
        static BibTeX bibtex();
    };

    ///////////////////////////////////////////

}  // namespace utopia

#endif  // UTOPIA_CITE_UTOPIA_HPP
