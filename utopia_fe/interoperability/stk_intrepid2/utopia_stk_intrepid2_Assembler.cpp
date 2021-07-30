#include "utopia_stk_intrepid2_Assembler.hpp"

#include "utopia_intrepid2_Field.hpp"
#include "utopia_intrepid2_Gradient.hpp"
#include "utopia_intrepid2_Mass.hpp"
#include "utopia_intrepid2_SubdomainFunction.hpp"
#include "utopia_intrepid2_Transport.hpp"
#include "utopia_stk_intrepid2.hpp"

#include "utopia_make_unique.hpp"

namespace utopia {
    namespace stk {
        class StkIntrepid2Assembler::Impl {
        public:
            std::shared_ptr<FunctionSpace> space;
            std::shared_ptr<Environment> environment;
        };

        void StkIntrepid2Assembler::set_environment(const std::shared_ptr<Environment> &env) {
            impl_->environment = env;
        }

        std::shared_ptr<StkIntrepid2Assembler::Environment> StkIntrepid2Assembler::environment() const {
            assert(impl_->environment);
            return impl_->environment;
        }

        void StkIntrepid2Assembler::set_space(const std::shared_ptr<FunctionSpace> &space) { impl_->space = space; }

        std::shared_ptr<FunctionSpace> StkIntrepid2Assembler::space() const { return impl_->space; }

        StkIntrepid2Assembler::~StkIntrepid2Assembler() = default;

        StkIntrepid2Assembler::StkIntrepid2Assembler(const std::shared_ptr<FE> &fe)
            : Super(fe), impl_(utopia::make_unique<Impl>()) {}

        ////////////////////////////////////////////////////

        class StkIntrepid2ProxyAssembler::Impl {
        public:
            std::shared_ptr<Intrepid2Assembler> assembler;
        };

        void StkIntrepid2ProxyAssembler::set_assembler(const std::shared_ptr<Intrepid2Assembler> &assembler) {
            impl_->assembler = assembler;
        }

        const std::shared_ptr<StkIntrepid2ProxyAssembler::Intrepid2Assembler> &StkIntrepid2ProxyAssembler::assembler()
            const {
            assert(impl_);
            return impl_->assembler;
        }

        StkIntrepid2ProxyAssembler::~StkIntrepid2ProxyAssembler() = default;

        StkIntrepid2ProxyAssembler::StkIntrepid2ProxyAssembler(const std::shared_ptr<FE> &fe)
            : Super(fe), impl_(utopia::make_unique<Impl>()) {}

        void StkIntrepid2ProxyAssembler::set_matrix_accumulator(
            const std::shared_ptr<TensorAccumulator> &matrix_accumulator) {
            Super::set_matrix_accumulator(matrix_accumulator);
            assert(impl_->assembler);
            impl_->assembler->set_matrix_accumulator(matrix_accumulator);
        }

        void StkIntrepid2ProxyAssembler::set_vector_accumulator(
            const std::shared_ptr<TensorAccumulator> &vector_accumulator) {
            Super::set_vector_accumulator(vector_accumulator);
            assert(impl_->assembler);
            impl_->assembler->set_vector_accumulator(vector_accumulator);
        }

        void StkIntrepid2ProxyAssembler::set_scalar_accumulator(
            const std::shared_ptr<TensorAccumulator> &scalar_accumulator) {
            Super::set_scalar_accumulator(scalar_accumulator);
            assert(impl_->assembler);
            impl_->assembler->set_scalar_accumulator(scalar_accumulator);
        }

        void StkIntrepid2ProxyAssembler::ensure_matrix_accumulator() {
            assembler()->ensure_matrix_accumulator();
            Super::set_matrix_accumulator(assembler()->matrix_accumulator());
        }

        void StkIntrepid2ProxyAssembler::ensure_vector_accumulator() {
            assembler()->ensure_vector_accumulator();
            Super::set_vector_accumulator(assembler()->vector_accumulator());
        }

        void StkIntrepid2ProxyAssembler::ensure_scalar_accumulator() {
            assembler()->ensure_scalar_accumulator();
            Super::set_scalar_accumulator(assembler()->scalar_accumulator());
        }

        bool StkIntrepid2ProxyAssembler::apply(const DynRankView &x, DynRankView &y) {
            return this->assembler()->apply(x, y);
        }

        bool StkIntrepid2ProxyAssembler::assemble_matrix() { return assembler()->assemble_matrix(); }

        /////////////////////////////////////////////////////////////////

    }  // namespace stk
}  // namespace utopia
