#include "utopia_stk_intrepid2_Discretization.hpp"

namespace utopia {

    class Discretization<stk_FS_t, stk_FE_t>::Impl {
    public:
        std::shared_ptr<stk_FE_t> fe;
    };

    Discretization<stk_FS_t, stk_FE_t>::Discretization(const std::shared_ptr<FunctionSpace> &space,
                                                       const std::shared_ptr<FE> &fe)
        : Super(space), impl_(utopia::make_unique<Impl>()) {
        impl_->fe = fe;
    }

    Discretization<stk_FS_t, stk_FE_t>::~Discretization() {}

    void Discretization<stk_FS_t, stk_FE_t>::read(Input &in) {
        // TODO
    }

    ////////////////////////////////////////////////////////////////////////////////////

    void Discretization<stk_FS_t, stk_FE_t>::create(std::vector<FE_ptr> &fe, int order, const Part &part) {
        // TODO
    }
    void Discretization<stk_FS_t, stk_FE_t>::create_on_boundary(FE &fe, int order, const Part &part) {
        // TODO
    }

    ////////////////////////////////////////////////////////////////////////////////////

    void Discretization<stk_FS_t, stk_FE_t>::convert_field(const Field &in,
                                                           std::vector<std::shared_ptr<FEField>> &out,
                                                           const Part &part) {
        // TODO
    }
    void Discretization<stk_FS_t, stk_FE_t>::convert_field(const std::vector<std::shared_ptr<FEField>> &in,
                                                           Field &out,
                                                           const Part &part) {
        // TODO
    }

    ////////////////////////////////////////////////////////////////////////////////////

    void Discretization<stk_FS_t, stk_FE_t>::global_to_local(const Vector &vector,
                                                             const std::vector<VectorAccumulator> &element_vectors,
                                                             const Part &part,
                                                             const int comp) {
        assert(false);
    }

    void Discretization<stk_FS_t, stk_FE_t>::local_to_global(const std::vector<MatrixAccumulator> &acc,
                                                             AssemblyMode mode,
                                                             Matrix &mat,
                                                             const Part &part) {
        // TODO
    }

    void Discretization<stk_FS_t, stk_FE_t>::local_to_global(const std::vector<VectorAccumulator> &acc,
                                                             AssemblyMode mode,
                                                             Vector &vec,
                                                             const Part &part) {
        // TODO
    }

    void Discretization<stk_FS_t, stk_FE_t>::local_to_global_on_boundary(const std::vector<VectorAccumulator> &acc,
                                                                         AssemblyMode mode,
                                                                         Vector &vec,
                                                                         const Part &part) {
        // TODO
    }

}  // namespace utopia
