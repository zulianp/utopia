#include "utopia_stk_intrepid2_Discretization.hpp"

#include "utopia_stk_intrepid2.hpp"

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
        // FIXME
        FE_ptr fe0;

        if (fe.size() == 1) {
            fe0 = fe[0];
        } else {
            fe0 = std::make_shared<FE>();
        }

        create_fe(*this->space(), *fe0, order);

        fe.clear();
        fe.push_back(fe0);
    }

    void Discretization<stk_FS_t, stk_FE_t>::create_on_boundary(std::vector<FE_ptr> &fe, int order, const Part &part) {
        // FIXME
        FE_ptr fe0;

        if (fe.size() == 1) {
            fe0 = fe[0];
        } else {
            fe0 = std::make_shared<FE>();
        }

        if (part.is_all()) {
            create_fe_on_boundary(*this->space(), *fe0, order);
        } else {
            create_fe_on_boundary(*this->space(), *fe0, part.name, order);
        }

        fe.clear();
        fe.push_back(fe0);
    }

    ////////////////////////////////////////////////////////////////////////////////////

    void Discretization<stk_FS_t, stk_FE_t>::convert_field(const Field &in,
                                                           std::vector<std::shared_ptr<FEField>> &out,
                                                           const Part &part) {
        // FIXME

        if (out.size() == 1) {
            utopia::convert_field(in, *out[0]);
        } else {
            assert(false);
        }
    }
    void Discretization<stk_FS_t, stk_FE_t>::convert_field(const std::vector<std::shared_ptr<FEField>> &in,
                                                           Field &out,
                                                           const Part &part) {
        // FIXME
        assert(in.size() == 1);
        utopia::convert_field(*in[0], out);
    }

    ////////////////////////////////////////////////////////////////////////////////////

    void Discretization<stk_FS_t, stk_FE_t>::global_to_local(const Vector &vector,
                                                             std::vector<VectorAccumulator> &element_vectors,
                                                             const Part &part,
                                                             const int comp) {
        // FIXME
        element_vectors.resize(1);
        utopia::global_to_local(*this->space(), vector, element_vectors[0], comp);
    }

    void Discretization<stk_FS_t, stk_FE_t>::local_to_global(const std::vector<MatrixAccumulator> &acc,
                                                             AssemblyMode mode,
                                                             Matrix &mat,
                                                             const Part &part) {
        // FIXME
        assert(acc.size() == 1);
        utopia::local_to_global(*this->space(), acc[0], mode, mat);
    }

    void Discretization<stk_FS_t, stk_FE_t>::local_to_global(const std::vector<VectorAccumulator> &acc,
                                                             AssemblyMode mode,
                                                             Vector &vec,
                                                             const Part &part) {
        // FIXME
        assert(acc.size() == 1);
        utopia::local_to_global(*this->space(), acc[0], mode, vec);
    }

    void Discretization<stk_FS_t, stk_FE_t>::local_to_global_on_boundary(const std::vector<VectorAccumulator> &acc,
                                                                         AssemblyMode mode,
                                                                         Vector &vec,
                                                                         const Part &part) {
        // TODO
        utopia::side_local_to_global(*this->space(), acc[0], mode, vec, part.name);
    }

}  // namespace utopia
