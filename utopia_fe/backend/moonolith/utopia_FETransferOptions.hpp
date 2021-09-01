#ifndef UTOPIA_FE_TRANSFER_OPTIONS_HPP
#define UTOPIA_FE_TRANSFER_OPTIONS_HPP

#include "utopia_Input.hpp"
#include "utopia_Options.hpp"

#include <cassert>
#include <utility>
#include <vector>

namespace utopia {

    // FIXME move somewhere public
    class FETransferOptions : public Configurable {
    public:
        virtual ~FETransferOptions() = default;

        void read(Input &in) override {
            if (!Options()  //
                     .add_option("verbose", verbose, "Verbose output during computation.")
                     .add_option(
                         "n_var", n_var, "Number of dimensions of vector function. Useful with tensor-product spaces.")
                     .add_option(
                         "use_reference_frame", use_reference_frame, "Use reference element frame for intersections.")
                     .add_option("clear_non_essential_matrices",
                                 clear_non_essential_matrices,
                                 "Keeps only the final transfer matrix in memory and deletes the rest.")
                     .add_option("export_tensors", export_tensors_, "Exports tensors to disk.")
                     .add_option("has_covering", has_covering, "Constructs lagrange multiplier in intersections.")
                     .add_option("chop_tol", chop_tol, "Chop numeric entries close below a certain tolerance")
                     .parse(in)) {
                return;
            }

            in.get("coupling", [this](Input &in) {
                in.get_all([this](Input &in) {
                    int from = -1, to = -1;

                    in.get("from", from);
                    in.get("to", to);

                    assert(from != -1);
                    assert(to != -1);

                    tags.emplace_back(from, to);
                });
            });
        }

#ifndef NDEBUG
        bool verbose{true};
#else
        bool verbose{false};
#endif
        bool has_covering{true};
        int n_var{1};
        std::vector<std::pair<int, int>> tags;
        bool clear_non_essential_matrices{true};
        bool export_tensors_{false};
        bool use_reference_frame{false};
        double chop_tol{0.};
    };

}  // namespace utopia

#endif  // UTOPIA_FE_TRANSFER_OPTIONS_HPP