#ifndef UTOPIA_COMPOSITE_LOCAL_ASSEMBLER_HPP
#define UTOPIA_COMPOSITE_LOCAL_ASSEMBLER_HPP

#include "utopia_LocalAssembler.hpp"

#include <iostream>
#include <memory>
#include "libmesh/dense_matrix.h"

namespace utopia {

    class MixedBidirectionalLocalAssembler final : public LocalAssembler {
    public:
        using Elem = LocalAssembler::Elem;
        using FEType = LocalAssembler::FEType;
        using Matrix = LocalAssembler::Matrix;

        MixedBidirectionalLocalAssembler(const std::shared_ptr<LocalAssembler> &master_to_slave,
                                         const std::shared_ptr<LocalAssembler> &slave_to_master);

        ~MixedBidirectionalLocalAssembler();

        bool assemble(const Elem &master,
                      FEType master_type,
                      const Elem &slave,
                      FEType slave_type,
                      Matrix &mat) override;

        bool assemble(const Elem &master,
                      FEType master_type,
                      const Elem &slave,
                      FEType slave_type,
                      std::vector<Matrix> &mat) override;

        int n_forms() const override;

        Type type(const int index) const override;

        void print_stats(std::ostream &os = std::cout) const override;
        void set_return_false_only_if_both_fail(const bool val) { return_false_only_if_both_fail_ = val; }

    private:
        std::vector<std::shared_ptr<LocalAssembler>> assemblers_;
        bool return_false_only_if_both_fail_;
    };
}  // namespace utopia

#endif  // UTOPIA_COMPOSITE_LOCAL_ASSEMBLER_HPP