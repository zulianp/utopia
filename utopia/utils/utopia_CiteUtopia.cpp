#include "utopia_CiteUtopia.hpp"
#include "utopia_AuthoredWork.hpp"

namespace utopia {
    BibTeX Cite<Utopia2021>::bibtex() 
    {
        return BibTeX(
            "@article{utopia2021,\n"
            "\tauthor = {Zulian, Patrick and Kopani{\\v{c}}{\\'{a}}kov{\\'{a}}, Alena and Nestola, Maria G C and\n" 
            "\t\tFadel, Nur and Fink, Andreas and VandeVondele, Joost and Krause, Rolf},\n"
            "\ttitle = {Large scale simulation of pressure induced phase‐field fracture propagation using {U}topia},\n"
            "\tjournal = {CCF Transactions on High Performance Computing},\n"
            "\tyear = {2021},\n"
            "\tmonth = {06},\n"
            "\tdoi = {10.1007/s42514-021-00069-6},\n"
            "\turl = {https://doi.org/10.1007/s42514-021-00069-6},\n"
            "\teprint = {https://doi.org/10.1007/s42514-021-00069-6},\n"
            "}\n");
    }

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

    BibTeX Cite<Kopanicakova2020Thesis>::bibtex() {
        return BibTeX(
            "@phdthesis{kopanicakova2021Thesis,\n"
            "\ttitle={Multilevel minimization in trust-region framework - Algorithmic and software developments,\n"
            "\tauthor={Kopani{\\v{c}}{\\'a}kov{\\'a}, Alena},\n"
            "\tyear={2021},\n"
            "\tschool={Universit\\`a della Svizzera italiana},\n"
            "}\n");
    }

}  // namespace utopia
