#include "utopia_mars_Discretization.hpp"

namespace utopia {

    class Discretization<mars_FS_t, mars_FE_t>::Impl {
    public:
        std::shared_ptr<mars_FE_t> fe;
    };

    Discretization<mars_FS_t, mars_FE_t>::Discretization(const std::shared_ptr<FunctionSpace> &space,
                                                         const std::shared_ptr<FE> &fe)
        : Super(space), impl_(utopia::make_unique<Impl>()) {
        impl_->fe = fe;
    }

    Discretization<mars_FS_t, mars_FE_t>::~Discretization() {}

    void Discretization<mars_FS_t, mars_FE_t>::read(Input & /*in*/) {
        // TODO
    }

    ////////////////////////////////////////////////////////////////////////////////////

    void Discretization<mars_FS_t, mars_FE_t>::create(std::vector<FE_ptr> &fe, int order, const Part &part) {
        this->space()->handler()->create(fe, order, part);
    }

    void Discretization<mars_FS_t, mars_FE_t>::create_on_boundary(std::vector<FE_ptr> &fe,
                                                                  int order,
                                                                  const Part &part) {
        this->space()->handler()->create_on_boundary(fe, order, part);
    }

    ////////////////////////////////////////////////////////////////////////////////////

    void Discretization<mars_FS_t, mars_FE_t>::convert_field(const Field &in,
                                                             std::vector<std::shared_ptr<FEField>> &out,
                                                             const Part &part) {
        this->space()->handler()->convert_field(in, out, part);
    }
    void Discretization<mars_FS_t, mars_FE_t>::convert_field(const std::vector<std::shared_ptr<FEField>> &in,
                                                             Field &out,
                                                             const Part &part) {
        this->space()->handler()->convert_field(in, out, part);
    }

    ////////////////////////////////////////////////////////////////////////////////////

    void Discretization<mars_FS_t, mars_FE_t>::global_to_local(const Vector &vector,
                                                               std::vector<VectorAccumulator> &element_vectors,
                                                               const Part &part,
                                                               const int comp) {
        this->space()->handler()->global_to_local(vector, element_vectors, part, comp);
    }

    void Discretization<mars_FS_t, mars_FE_t>::local_to_global(const std::vector<MatrixAccumulator> &acc,
                                                               AssemblyMode mode,
                                                               Matrix &mat,
                                                               const Part &part) {
        this->space()->handler()->local_to_global(acc, mode, mat, part);
    }

    void Discretization<mars_FS_t, mars_FE_t>::local_to_global(const std::vector<VectorAccumulator> &acc,
                                                               AssemblyMode mode,
                                                               Vector &vec,
                                                               const Part &part) {
        this->space()->handler()->local_to_global(acc, mode, vec, part);
    }

    void Discretization<mars_FS_t, mars_FE_t>::local_to_global_on_boundary(const std::vector<VectorAccumulator> &acc,
                                                                           AssemblyMode mode,
                                                                           Vector &vec,
                                                                           const Part &part) {
        this->space()->handler()->local_to_global_on_boundary(acc, mode, vec, part);
    }

}  // namespace utopia
