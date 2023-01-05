#ifndef UTOPIA_KOKKOS_MASS_NEW_HPP
#define UTOPIA_KOKKOS_MASS_NEW_HPP

#include "utopia_fe_base.hpp"

#include "utopia_kokkos_FE.hpp"
#include "utopia_kokkos_FEAssembler.hpp"
#include "utopia_kokkos_Field.hpp"
#include "utopia_kokkos_SubdomainValue.hpp"

#include "utopia_kokkos_MassOp.hpp"
#include "utopia_kokkos_Material.hpp"

#include "utopia_Traits.hpp"

namespace utopia {

    namespace kokkos {

        template <class FunctionSpace, class FE_, class DiffusionCoefficient = typename FE_::Scalar>
        class MassNew : public utopia::Material<FunctionSpace, FE_> {
        public:
            using FE = FE_;

            using Super = utopia::Material<FunctionSpace, FE_>;
            using Scalar = typename Traits<FunctionSpace>::Scalar;
            using DynRankView = typename Super::Field::DynRankView;

            using Op = utopia::kokkos::kernels::MassOp<Scalar, Scalar, typename FE::Function, typename FE::Measure>;
            using LumpedOp = utopia::kokkos::kernels::LumpedOp<Op>;
            using Lump = utopia::kokkos::kernels::Lump<DynRankView>;
            using SubdomainValue = utopia::kokkos::SubdomainValue<FE>;

            class Params : public Configurable {
            public:
                void read(Input &in) override {
                    in.get("density", density);
                    in.get("n_components", n_components);
                    in.get("lumped", lumped);
                    in.get("verbose", verbose);

                    in.get("density_function", [this](Input &node) {
                        density_function = std::make_shared<SubdomainValue>(1.0);
                        node.get("function", [this](Input &inner_node) { density_function->read(inner_node); });
                    });

                    if (verbose) {
                        utopia::out() << "-----------------------------\n";
                        utopia::out() << "Mass\n";
                        utopia::out() << "lumped:\t" << lumped << '\n';
                        utopia::out() << "density:\t" << density << '\n';

                        if (n_components) {
                            utopia::out() << "n_components:\t" << n_components << '\n';
                        } else {
                            utopia::out() << "n_components:\tauto\n";
                        }

                        utopia::out() << "-----------------------------\n";
                    }
                }

                Params(const Scalar &density = Scalar(1.0)) : density(density) {}
                Params(const Params &) = default;

                Scalar density;
                int n_components{0};
                bool lumped{false};
                std::shared_ptr<SubdomainValue> density_function;

                // Testing an printing
                bool verbose{false};
            };

            void read(Input &in) override {
                Super::read(in);
                params_.read(in);
            }

            void initialize(const std::shared_ptr<FunctionSpace> &space) {
                Super::initialize(space);
                fix_params(*space);
            }

            void fix_params(FunctionSpace &space) {
                if (params_.n_components == 0) {
                    params_.n_components = space.n_var();
                }
            }

            MassNew(Params op = Params()) : Super(), params_(std::move(op)) {}

            inline int n_vars() const override { return params_.n_components; }
            inline std::string name() const override { return "Mass"; }

            inline bool has_hessian() const override { return true; }
            inline bool is_linear() const override { return true; }
            inline bool is_operator() const override { return true; };

            inline bool has_gradient() const override { return false; }
            inline bool has_value() const override { return false; }

            inline Op make_op() {
                auto &&assembler = this->assembler();
                assert(assembler);
                auto fe = assembler->fe();
                assert(params_.n_components == 1);

                return Op(params_.density, fe.fun(), fe.measure(), params_.n_components);
            }

            inline LumpedOp make_lumped_op() { return LumpedOp(make_op()); }

            template <class TOp, int BlockSize>
            class BlockOp : public TestTrialOp {
            public:
                static const int NComponentsTest = BlockSize;
                static const int NComponentsTrial = BlockSize;

                BlockOp(const TOp &op) : op(op) {}

                UTOPIA_INLINE_FUNCTION void operator()(const int cell,
                                                       const int i,
                                                       const int j,
                                                       StaticMatrix<Scalar, BlockSize, BlockSize> &block) const {
                    for (int subi = 0; subi < BlockSize; ++subi) {
                        for (int subj = 0; subj < BlockSize; ++subj) {
                            block(subi, subj) = op(cell, i, j, subi, subj);
                        }
                    }
                }

                TOp op;
            };

            template <int BlockSize, class TOp>
            inline BlockOp<TOp, BlockSize> make_block_op(const TOp &op) {
                return BlockOp<TOp, BlockSize>(op);
            }

            template <class TOp, int BlockSize>
            class DiagBlockOp : public TestTrialOp {
            public:
                static const int NComponentsTest = BlockSize;
                static const int NComponentsTrial = BlockSize;

                DiagBlockOp(const TOp &op) : op(op) {}

                UTOPIA_INLINE_FUNCTION void operator()(const int cell,
                                                       const int i,
                                                       const int j,
                                                       StaticMatrix<Scalar, BlockSize, BlockSize> &block) const {
                    Scalar val = op(cell, i);
                    for (int subi = 0; subi < BlockSize; ++subi) {
                        block(subi, subi) = val;
                    }
                }

                TOp op;
            };

            template <int BlockSize, class TOp>
            inline DiagBlockOp<TOp, BlockSize> make_diag_block_op(const TOp &op) {
                return DiagBlockOp<TOp, BlockSize>(op);
            }

            template <class TOp, class ElementTags, int BlockSize>
            class ScaledBlockOp : public TestTrialOp {
            public:
                static const int NComponentsTest = BlockSize;
                static const int NComponentsTrial = BlockSize;

                ScaledBlockOp(const TOp &op, const ElementTags &tags, SubdomainValue sv) : op(op), tags(tags), sv(sv) {}

                UTOPIA_INLINE_FUNCTION void operator()(const int cell,
                                                       const int i,
                                                       const int j,
                                                       StaticMatrix<Scalar, BlockSize, BlockSize> &block) const {
                    op(cell, i, j, block);

                    auto tag = tags(cell);
                    auto val = sv.value(tag);
                    block *= val;
                }

                TOp op;
                ElementTags tags;
                SubdomainValue sv;
            };

            template <int BlockSize, class TOp, class ElementTags>
            inline ScaledBlockOp<TOp, ElementTags, BlockSize> make_scaled_op(const TOp &op,
                                                                             const ElementTags &tags,
                                                                             const SubdomainValue &sv) {
                return ScaledBlockOp<TOp, ElementTags, BlockSize>(op, tags, sv);
            }

            template <int BlockSize>
            void aux_matrix_assemble(AssemblyMode mode) {
                auto &&assembler = this->assembler();
                assert(assembler);

                if (params_.density_function) {
                    auto &&fe = assembler->fe();
                    auto &&tags = fe.element_tags();

                    if (params_.lumped) {
                        assembler->assemble_matrix_eij_block(
                            "MassNew::hessian",
                            mode,
                            make_scaled_op<BlockSize>(
                                make_diag_block_op<BlockSize>(make_lumped_op()), tags, *params_.density_function));
                    } else {
                        assembler->assemble_matrix_eij_block(
                            "MassNew::hessian",
                            mode,
                            make_scaled_op<BlockSize>(
                                make_block_op<BlockSize>(make_op()), tags, *params_.density_function));
                    }
                } else {
                    if (params_.lumped) {
                        assembler->assemble_matrix_eij_block(
                            "MassNew::hessian", mode, make_diag_block_op<BlockSize>(make_lumped_op()));
                    } else {
                        assembler->assemble_matrix_eij_block(
                            "MassNew::hessian", mode, make_block_op<BlockSize>(make_op()));
                    }
                }
            }

            bool hessian_assemble(AssemblyMode mode) override {
                UTOPIA_TRACE_REGION_BEGIN("MassNew::hessian");

                switch (params_.n_components) {
                    case 1: {
                        aux_matrix_assemble<1>(mode);
                        break;
                    }
                    case 2: {
                        aux_matrix_assemble<2>(mode);
                        break;
                    }
                    case 3: {
                        aux_matrix_assemble<3>(mode);
                        break;
                    }
                    case 4: {
                        aux_matrix_assemble<4>(mode);
                        break;
                    }
                    default: {
                        assert(false && "Add case!");
                        Utopia::Abort();
                    }
                }

                UTOPIA_TRACE_REGION_END("MassNew::hessian");
                return true;
            }

            bool value_assemble(AssemblyMode mode) override {
                assert(false);
                return false;
            }

            template <int BlockSize>
            void aux_apply_assemble(utopia::kokkos::Field<FE> &field, AssemblyMode mode) {
                auto &&assembler = this->assembler();
                assert(assembler);

                if (params_.density_function) {
                    auto &&fe = assembler->fe();
                    auto &&tags = fe.element_tags();

                    if (params_.lumped) {
                        assembler->assemble_apply_ei_block(
                            "MassNew::apply",
                            mode,
                            make_scaled_op<BlockSize>(
                                make_diag_block_op<BlockSize>(make_lumped_op()), tags, *params_.density_function),
                            field);
                    } else {
                        assembler->assemble_apply_ei_block(
                            "MassNew::apply",
                            mode,
                            make_scaled_op<BlockSize>(
                                make_block_op<BlockSize>(make_op()), tags, *params_.density_function),
                            field);
                    }
                } else {
                    if (params_.lumped) {
                        assembler->assemble_apply_ei_block(
                            "MassNew::apply", mode, make_diag_block_op<BlockSize>(make_lumped_op()), field);
                    } else {
                        assembler->assemble_apply_ei_block(
                            "MassNew::apply", mode, make_block_op<BlockSize>(make_op()), field);
                    }
                }
            }

            // Matrix free hessian application
            bool apply_assemble(utopia::kokkos::Field<FE> &field, AssemblyMode mode) override {
                UTOPIA_TRACE_REGION_BEGIN("MassNew::apply");

                switch (params_.n_components) {
                    case 1: {
                        aux_apply_assemble<1>(field, mode);
                        break;
                    }
                    case 2: {
                        aux_apply_assemble<2>(field, mode);
                        break;
                    }
                    case 3: {
                        aux_apply_assemble<3>(field, mode);
                        break;
                    }
                    case 4: {
                        aux_apply_assemble<4>(field, mode);
                        break;
                    }
                    default: {
                        assert(false && "Add case!");
                        Utopia::Abort();
                    }
                }

                UTOPIA_TRACE_REGION_END("MassNew::apply");
                return true;
            }

            // NVCC_PRIVATE :
            Params params_;
        };

    }  // namespace kokkos
}  // namespace utopia

#endif  // UTOPIA_KOKKOS_MASS_NEW_HPP
