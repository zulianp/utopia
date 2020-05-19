#ifndef UTOPIA_AUTHORED_WORK_HPP
#define UTOPIA_AUTHORED_WORK_HPP

#include "utopia_Describable.hpp"
#include "utopia_MPI.hpp"

#include <string>
#include <vector>

namespace utopia {

    class BibTeX final {
    public:
        inline void set_citation(const char *cite) { cite_ = cite; }
        inline void describe(std::ostream &os) const { os << cite_; }
        explicit constexpr BibTeX(const char *str) : cite_(str) {}

    private:
        const char *cite_;
    };

    template <class Algorithm>
    class Cite {
    public:
        static constexpr BibTeX bibtex() { return BibTeX("undefined"); }
    };

    class CitationsDB : public Describable {
    public:
        inline void describe(std::ostream &os = std::cout) const override {
            if (authorships_.empty()) return;

            if (mpi_world_rank() == 0) {
                os << "\n\n";
                os << "---------------------------------------------------------------------------------------------\n";
                os << "If you produced results with this run for your article, proceeding, presentation, etc... please "
                      "cite the following work:\n\n";

                for (const auto &a : authorships_) {
                    a.describe(os);
                    os << "\n";
                }

                os << "---------------------------------------------------------------------------------------------\n"
                      "\n";
            }
        }

        inline static void print() { instance().describe(); }

        inline static CitationsDB &instance() {
            static CitationsDB instance_;
            return instance_;
        }

        inline void cite(BibTeX bib) { authorships_.push_back(std::move(bib)); }

    private:
        std::vector<BibTeX> authorships_;
    };

    template <class Algorithm>
    class AuthoredWork {
    public:
        static bool cited;

        AuthoredWork() { cite_init(); }

        void cite_init() {
            if (!cited) {
                cited = true;
                CitationsDB::instance().cite(Cite<Algorithm>::bibtex());
            }
        }

        // virtual BibTeX bibtex() const = 0;
    };  // namespace utopia

    template <class Algorithm>
    bool AuthoredWork<Algorithm>::cited = false;

}  // namespace utopia

#endif  // UTOPIA_AUTHORED_WORK_HPP
