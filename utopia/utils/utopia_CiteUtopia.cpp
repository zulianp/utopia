#include "utopia_CiteUtopia.hpp"
#include "utopia_AuthoredWork.hpp"

namespace utopia {

    BibTeX Cite<Utopia2016Git>::bibtex() {
        return BibTeX(
            "@misc{utopia2016git,\n"
            "\tauthor = {Patrick Zulian and Alena Kopani{\\v c}{\\'a}kov{\\'a} and Maria Chiara Giuseppina Nestola "
            "and\n"
            "\t          Andreas Fink and Nur Fadel and Alessandro Rigazzi and\n"
            "\t          Victor Magri and Teseo Schneider and Eric Botter and Jan Mankau and Rolf Krause},\n"
            "\ttitle = {{U}topia: A performance portable {C}++ library for parallel linear and nonlinear algebra. "
            "{G}it "
            "repository},\n"
            "\thowpublished = {https://bitbucket.org/zulianp/utopia},\n"
            "\tyear = {2016},\n"
            "}\n");
    }

    BibTeX Cite<Kopanicakova2020Recursive>::bibtex() {
        return BibTeX(
            "@article{kopanicakova2020recursive,\n"
            "\ttitle={A recursive multilevel trust region method with application to fully monolithic phase-field "
            "models of brittle fracture},\n"
            "\tauthor={Kopani{\\v{c}}{\\'a}kov{\\'a}, Alena and Krause, Rolf},\n"
            "\tjournal={Computer Methods in Applied Mechanics and Engineering},\n"
            "\tvolume={360},\n"
            "\tpages={112720},\n"
            "\tyear={2020},\n"
            "\tpublisher={Elsevier}\n"
            "}\n");
    }

}  // namespace utopia
