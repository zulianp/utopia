#ifndef UTOPIA_SEMI_GEOMETRIC_HIERARCHY_HPP
#define UTOPIA_SEMI_GEOMETRIC_HIERARCHY_HPP

#include "utopia_ILU.hpp"
#include "utopia_MonotoneMultigrid.hpp"

#include "utopia_fe_Core.hpp"

#include "utopia_BoundingBoxMultiLevel.hpp"

#include <memory>
#include <vector>

namespace utopia {

    template <class FunctionSpace>
    class SemiGeometricHierarchy {
    public:
        using Vector = typename Traits<FunctionSpace>::Vector;
        using Matrix = typename Traits<FunctionSpace>::Matrix;
        using Scalar = typename Traits<FunctionSpace>::Scalar;
        using SizeType = typename Traits<FunctionSpace>::SizeType;
        using IndexSet = typename Traits<FunctionSpace>::IndexSet;
        using Comm = typename Traits<FunctionSpace>::Communicator;

        using Transfer = utopia::Transfer<Matrix, Vector>;
        using IPTransfer = utopia::IPTransfer<Matrix, Vector>;
        using CoarseSpaceGen = utopia::BoundingBoxMultiLevelFunctionSpaceGenerator<FunctionSpace>;

        virtual ~SemiGeometricHierarchy() = default;

        void generate_coarse_spaces(Input &in, FunctionSpace &finest) {
            int n_auto_levels = space_hierarchy_.empty() ? 1 : 0;
            in.get("n_auto_levels", n_auto_levels);
            if (n_auto_levels == 0) return;

            std::vector<std::shared_ptr<FunctionSpace>> coarse_spaces;
            CoarseSpaceGen gen;
            gen.read(in);

            if (gen.generate_spaces(finest, n_auto_levels, coarse_spaces)) {
                space_hierarchy_.insert(space_hierarchy_.begin(), coarse_spaces.begin(), coarse_spaces.end());
            } else {
                Utopia::Abort("Failed to create coarse spaces!");
            }
        }

        void read_hierarchy(Input &in) {
            if (fine_space_) {
                block_size_ = fine_space_->n_var();
            }

            bool export_coarse_meshes = false;
            in.get("export_coarse_meshes", export_coarse_meshes);
            in.get("block_size", block_size_);

            in.get("coarse_spaces", [this](Input &array_node) {
                array_node.get_all([this](Input &node) {
                    auto space = std::make_shared<FunctionSpace>();
                    space->read(node);

                    assert(!space->empty());
                    space_hierarchy_.push_back(space);
                });
            });

            in.get("clear_spaces_after_init", clear_spaces_after_init_);

            if (space_hierarchy_.empty()) {
                generate_coarse_spaces(in, *fine_space_);
            } else {
                generate_coarse_spaces(in, *space_hierarchy_[0]);
            }

            if (export_coarse_meshes) {
                int num_space = 0;
                for (auto &s : space_hierarchy_) {
                    s->mesh().write("mesh_" + std::to_string(num_space++) + ".e");
                }
            }
        }

        template <class FineLevelTransfer>
        bool generate_transfer_operators(std::vector<std::shared_ptr<Transfer>> &transfers) {
            const int n_levels = space_hierarchy_.size();

            if (n_levels == 0) {
                assert(false);
                return false;
            }

            if (!fine_space_) {
                assert(false);
                Utopia::Abort("SemiGeometricHierarchy fine space must be set!");
            }

            transfers.clear();

            InputParameters params;
            params.set("n_var", fine_space_->n_var());
            params.set("chop_tol", 1e-8);

            FETransfer<FunctionSpace> transfer;
            transfer.read(params);

            for (int l = 0; l < n_levels - 1; ++l) {
                if (!transfer.init(space_hierarchy_[l], space_hierarchy_[l + 1])) {
                    assert(false);
                    Utopia::Abort("SemiGeometricHierarchy failed to set-up transfer operator!");
                }

                transfers.push_back(transfer.template build_transfer<IPTransfer>());
            }

            if (!transfer.init(space_hierarchy_[n_levels - 1], fine_space_)) {
                assert(false);
                Utopia::Abort("SemiGeometricHierarchy failed to set-up transfer operator!");
            }

            auto t = transfer.template build_transfer<FineLevelTransfer>();
            fine_space_->apply_constraints(*t->I_ptr(), 0.0);

            transfers.push_back(t);

            if (clear_spaces_after_init_) {
                space_hierarchy_.clear();
            }

            return true;
        }

        inline std::shared_ptr<FunctionSpace> &fine_space() { return fine_space_; }
        inline void set_fine_space(const std::shared_ptr<FunctionSpace> &fine_space) { fine_space_ = fine_space; }
        inline std::vector<std::shared_ptr<FunctionSpace>> &space_hierarchy() { return space_hierarchy_; }

        inline int block_size() const { return block_size_; };
        inline int n_coarse_spaces() const { return space_hierarchy_.size(); }

    private:
        std::shared_ptr<FunctionSpace> fine_space_;
        std::vector<std::shared_ptr<FunctionSpace>> space_hierarchy_;
        int block_size_{1};
        bool clear_spaces_after_init_{false};
    };

}  // namespace utopia

#endif  // UTOPIA_SEMI_GEOMETRIC_HIERARCHY_HPP
