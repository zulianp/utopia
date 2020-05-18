#ifndef UTOPIA_AUTHORED_WORK_HPP
#define UTOPIA_AUTHORED_WORK_HPP

#include "utopia_Describable.hpp"
#include "utopia_MPI.hpp"

#include <string>
#include <vector>

namespace utopia {

    class BibTeX final : public Describable {
    public:
        inline void set_citation(std::string cite) { cite_ = std::move(cite); }
        inline void describe(std::ostream &os) const override { os << cite_; }

        BibTeX(const char *str) : cite_(str) {}

    private:
        std::string cite_;
    };

    class CitationsDB : public Describable {
    public:
        inline void describe(std::ostream &os = std::cout) const override {
            if (authorships_.empty()) return;

            if (mpi_world_rank() == 0) {
                os << "---------------------------------------------------\n";
                os << "If you use this runs for your articles, journals, presentations, etc. Please cite the following "
                      "work:\n";

                for (const auto &a : authorships_) {
                    a.describe(os);
                }

                os << "---------------------------------------------------\n";
            }
        }

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

        AuthoredWork() {
            if (!cited) {
                cited = true;
                CitationsDB::instance().cite(std::move(bibtex()));
            }
        }

        virtual BibTeX bibtex() const = 0;
    };

    template <class Algorithm>
    bool AuthoredWork<Algorithm>::cited = false;

}  // namespace utopia

#endif  // UTOPIA_AUTHORED_WORK_HPP
