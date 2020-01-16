#include "utopia_Tpetra_Vector.hpp"
#include "utopia_Logger.hpp"
#include "utopia_Instance.hpp"
#include "utopia_trilinos_Utils.hpp"
#include "utopia_kokkos_Eval_Reduce.hpp"

#include <Tpetra_CrsMatrix_decl.hpp>
#include <MatrixMarket_Tpetra.hpp>

#include <Kokkos_Core.hpp>

#include <cmath>

//FIXME
// - ghosted vector has problematic behaviour e.g.: norm2(ghosted_vec) == norm2(offset_view(ghosted_vec)) which is wrong
// maybe this can help at least for assembly #include <Tpetra_MultiVectorFiller.hpp> or FEMultiVector
namespace utopia {

    template<typename Integer>
    void TpetraVector::generic_select(
        const std::vector<Integer> &index,
        TpetraVector &out
        ) const
    {
        //FIXME does not work in parallel
        auto data = implementation().getData();
        const SizeType offset = implementation().getMap()->getMinGlobalIndex();

        auto map = Teuchos::rcp(
            new map_type(
                Teuchos::OrdinalTraits<Tpetra::global_size_t>::invalid(),
                index.size(),
                0,
                communicator())
        );

        out.init(map);

        if(comm().size() == 1) {
            const SizeType n = index.size();
            auto out_data = out.implementation().getDataNonConst();

            for(SizeType i = 0; i < n; ++i) {
                const SizeType local_index = index[i] - offset;

                if(local_index >= n) {
                    std::cout << local_index << " = " << index[i] << " - " << offset << std::endl;
                }

                assert(local_index < n);
                out_data[i] = data[local_index];
            }
        } else {
            assert(false && "IMPLEMENT ME");
        //     /////////////////////////////////////////////////

        //     std::vector<GO> tpetra_index;
        //     tpetra_index.reserve(index.size());

        //     for(auto i : index) {
        //         tpetra_index.push_back(i);
        //     }

        //     const Teuchos::ArrayView<const GO>
        //        index_view(tpetra_index);

        //      auto import_map = Teuchos::rcp(new map_type(global_size, index_view, 0, comm));

        //     Tpetra::Import<
        //         LO,
        //         GO,
        //         vector_type::node_type> importer(map, import_map);

        //     implementation().doImport(out.implementation(), importer, Tpetra::INSERT);
        }
    }

    void TpetraVector::select(const IndexSet &index, TpetraVector &result) const
    {
        generic_select(index, result);
    }

    void TpetraVector::transform(const Sqrt &op)
    {
        apply(op);
    }

    void TpetraVector::transform(const Pow2 &op)
    {
        apply(op);
    }

    void TpetraVector::transform(const Log &op)
    {
        apply(op);
    }

    void TpetraVector::transform(const Exp &op)
    {
        apply(op);
    }

    void TpetraVector::transform(const Cos &op)
    {
        apply(op);
    }

    void TpetraVector::transform(const Sin &op)
    {
        apply(op);
    }

    void TpetraVector::transform(const Abs &op)
    {
        apply(op);
    }

    void TpetraVector::transform(const Minus &op)
    {
        apply(op);
    }

    void TpetraVector::transform(const Pow &op)
    {
        apply(op);
    }

    void TpetraVector::transform(const Reciprocal<Scalar> &op)
    {
        apply(op);
    }

    void TpetraVector::add_vector(
        const IndexArray &indices,
        const ScalarArray &values)
    {
        const std::size_t n = values.size();
        assert(n == indices.size());

        for(std::size_t i = 0; i < n; ++i) {
            // add(indices[i], values[i]);
            if(!ghosted_vec_.is_null()) {
                ghosted_vec_->sumIntoGlobalValue(indices[i], values[i]);
            } else {
                implementation().sumIntoGlobalValue(indices[i], values[i]);
            }

        }
    }

    void TpetraVector::set_vector(
        const IndexArray &indices,
        const ScalarArray &values)
    {
        const std::size_t n = values.size();
        assert(n == indices.size());

        for(std::size_t i = 0; i < n; ++i) {
            // set(indices[i], values[i]);

            // if(!write_data_.is_null()) {
            //     write_data_[indices[i] - implementation().getMap()->getMinGlobalIndex()] = values[i];
            // } else {
                implementation().replaceGlobalValue(indices[i], values[i]);
            // }
        }
    }

    bool TpetraVector::read(const Teuchos::RCP< const Teuchos::Comm< int > > &comm, const std::string &path)
    {
        typedef Tpetra::Operator<>::scalar_type SC;
        typedef Tpetra::Operator<SC>::local_ordinal_type LO;
        typedef Tpetra::Operator<SC, LO>::global_ordinal_type GO;

        typedef Tpetra::CrsMatrix<SC, LO, GO, NT>         crs_matrix_type;

        try {
            rcp_map_type map;
            vec_ = Tpetra::MatrixMarket::Reader<crs_matrix_type>::readVectorFile(path, comm, map);
        } catch(const std::exception &ex) {
            // is.close();
            std::cout << ex.what() << std::endl;
            return false;
        }

        return !vec_.is_null();
    }

    bool TpetraVector::write(const std::string &path) const {
        typedef Tpetra::CrsMatrix<SC, LO, GO, NT>            crs_matrix_type;
        if(vec_.is_null()) return false;

        try {
            Tpetra::MatrixMarket::Writer<crs_matrix_type>::writeDenseFile(path, vec_, "vec", "");
        } catch(const std::exception &ex) {
            std::cout << ex.what() << std::endl;
            return false;
        }

        return true;
    }

    TpetraVector::Scalar TpetraVector::sum() const
    {
        Scalar ret = KokkosEvalReduce<TpetraVector, Plus>::eval(*this, Plus(), Scalar(0.));
        auto &comm = *communicator();
        Scalar ret_global = 0.;
        Teuchos::reduceAll(comm, Teuchos::REDUCE_SUM, 1, &ret, &ret_global);
        return ret_global;
    }

    TpetraVector::Scalar TpetraVector::min() const
    {
        Scalar ret = KokkosEvalReduce<TpetraVector, Min>::eval(*this, Min(), std::numeric_limits<Scalar>::max());
        auto &comm = *communicator();
        Scalar ret_global = ret;
        Teuchos::reduceAll(comm, Teuchos::REDUCE_MIN, 1, &ret, &ret_global);
        return ret_global;
    }

    TpetraVector::Scalar TpetraVector::max() const
    {
        Scalar ret = KokkosEvalReduce<TpetraVector, Max>::eval(*this, Max(), -std::numeric_limits<Scalar>::max());
        auto &comm = *communicator();
        Scalar ret_global = ret;
        Teuchos::reduceAll(comm, Teuchos::REDUCE_MAX, 1, &ret, &ret_global);
        return ret_global;
    }

    bool TpetraVector::has_nan_or_inf() const
    {
        int ret = KokkosEvalReduce<TpetraVector, IsNaNOrInf>::eval(*this, IsNaNOrInf(), Scalar(0));
        auto &comm = *communicator();
        int ret_global = 0;

        Teuchos::reduceAll(comm, Teuchos::REDUCE_MAX, 1, &ret, &ret_global);
        return ret_global;
    }

    void TpetraVector::ghosted(
        const rcp_comm_type &comm,
        const TpetraVector::GO &local_size,
        const TpetraVector::GO &global_size,
        const IndexArray &ghost_index
    )
    {
        rcp_map_type map(new map_type(global_size, local_size, 0, comm));
        rcp_map_type ghost_map;

        if(comm->getSize() != 0) {
            auto r_begin = local_size == 0 ? 0 : map->getMinGlobalIndex();
            auto r_end   = local_size == 0 ? 0 : (map->getMaxGlobalIndex() + 1);

            Range r = { r_begin, r_end };

            IndexArray filled_with_local;
            filled_with_local.reserve(r.extent() + ghost_index.size());

            for(auto i = r.begin(); i != r.end(); ++i) {
                filled_with_local.push_back(i);
            }

            for(auto g : ghost_index) {
                if(!r.inside(g)) {
                    filled_with_local.push_back(g);
                }
            }

            const Teuchos::ArrayView<const GO>
               local_indices(filled_with_local);

             GO local_entries = local_indices.size();
             GO total_entries = 0;

             Teuchos::reduceAll(*comm, Teuchos::REDUCE_SUM, 1, &local_entries, &total_entries);

             ghost_map = Teuchos::rcp(new map_type(total_entries, local_indices, 0, comm));

        } else {
            ghost_map = map;
        }

        UTOPIA_REPORT_ALLOC("TpetraVector::ghosted");
        ghosted_vec_ = Teuchos::rcp(new vector_type(ghost_map, 1));
        vec_ = ghosted_vec_->offsetViewNonConst(map, 0);

        assert(!vec_.is_null());
    }

    void TpetraVector::update_ghosts()
    {
        if(!has_ghosts()) return;

        auto map = vec_->getMap();
        auto ghost_map = ghosted_vec_->getMap();

        Tpetra::Import<
            LO,
            GO,
            vector_type::node_type> importer(map, ghost_map);

        ghosted_vec_->doImport(*vec_, importer, Tpetra::INSERT);
    }

    void TpetraVector::export_ghosts_add()
    {
        if(!has_ghosts()) return;

        auto map = vec_->getMap();
        auto ghost_map = ghosted_vec_->getMap();

        Tpetra::Export<
            LO,
            GO,
            vector_type::node_type> exporter(ghost_map, map);

        //FIXME can we avoid this? 
        //UTOPIA_REPORT_ALLOC(TpetraVector::export_ghosts_add)
        Teuchos::RCP<vector_type> y(new vector_type(map, 1));

        y->doExport(*ghosted_vec_, exporter, Tpetra::ADD);

        Tpetra::Import<vector_type::local_ordinal_type,
                        vector_type::global_ordinal_type,
                        vector_type::node_type> importer(map, ghost_map);

         ghosted_vec_->doImport(*y, importer, Tpetra::INSERT);
    }

    TpetraVector::TpetraVector(const TpetraVector &other)
    {
        copy(other);
    }

    void TpetraVector::copy(const TpetraVector &other)
    {
        if(&other == this) return;
        if(other.empty()) {
            clear();
            return;
        }

        if(other.has_ghosts()) {
            UTOPIA_REPORT_ALLOC("TpetraVector::copy");
            ghosted_vec_ = Teuchos::rcp(new vector_type(other.ghosted_vec_->getMap(), 1));
            ghosted_vec_->assign(*other.ghosted_vec_);
            vec_ = ghosted_vec_->offsetViewNonConst(other.vec_->getMap(), 0);
        } else {
            if(vec_.is_null() || other.size() != size()) {
                UTOPIA_REPORT_ALLOC("TpetraVector::copy");
                vec_ = (Teuchos::rcp(new vector_type(*other.vec_, Teuchos::Copy)));
            } else {
                vec_->assign(other.implementation());
            }
        }
    }

    void TpetraVector::assign(const TpetraVector &other)
    {
        if(this == &other) return;

        if(other.is_null()) {
            vec_.reset();
            return;
        }

        copy(other);
        return;
    }

    void TpetraVector::assign(TpetraVector &&other)
    {
        if(this == &other) return;
        comm_ = std::move(other.comm_);
        vec_ = std::move(other.vec_);
        ghosted_vec_ = std::move(other.ghosted_vec_);
    }

    TpetraVector &TpetraVector::operator=(const TpetraVector &other)
    {
        if(this == &other) return *this;

        if(other.is_null()) {
            vec_.reset();
            return *this;
        }

        copy(other);
        return *this;
    }

    void TpetraVector::write_unlock(WriteMode mode)
    {
        switch(mode) {
            case utopia::GLOBAL_ADD: {
                export_ghosts_add();
                break;
            }

            case utopia::LOCAL: {
                break;
            }

            default: {
                update_ghosts();
                break;
            }
        }

        // write_data_ = Teuchos::ArrayRCP<Scalar>();
        free_view();
    }

    bool TpetraVector::equals(const TpetraVector &other, const Scalar &tol) const
    {
        TpetraVector diff = *this;
        diff.axpy(-1, other);
        return diff.norm_infty() <= tol;
    }

    void TpetraVector::clear()
    {
        vec_.reset();
        ghosted_vec_.reset();
        view_ptr_.reset();
    }

    void TpetraVector::e_div(const TpetraVector &other)
    {
        KokkosEvalBinary<TpetraVector, Divides>::eval(*this, Divides(), other, *this);
    }

    void TpetraVector::e_min(const TpetraVector &other)
    {
        KokkosEvalBinary<TpetraVector, Min>::eval(*this, Min(), other, *this);
    }

    void TpetraVector::e_max(const TpetraVector &other)
    {
        KokkosEvalBinary<TpetraVector, Max>::eval(*this, Max(), other, *this);
    }

    void TpetraVector::e_div(const Scalar &other)
    {
       KokkosEvalBinary<TpetraVector, Divides>::eval(*this, Divides(), other, *this);
    }

    void TpetraVector::e_min(const Scalar &other)
    {
        KokkosEvalBinary<TpetraVector, Min>::eval(*this, Min(), other, *this);
    }

    void TpetraVector::e_max(const Scalar &other)
    {
        KokkosEvalBinary<TpetraVector, Max>::eval(*this, Max(), other, *this);
    }

    void TpetraVector::values(
        const rcp_comm_type &comm,
        const SizeType &n_local,
        const SizeType &n_global,
        const Scalar &val
    )
    {
        assert(n_local <= n_global);
        assert(n_local >= 0);

        rcp_map_type map;

        UTOPIA_REPORT_ALLOC("TpetraVector::values");
        map = Teuchos::rcp(new map_type(n_global, n_local, 0, comm));
        vec_.reset(new vector_type(map));
        implementation().putScalar(val);

        assert(size() > 0);
    }

    void TpetraVector::values(const SizeType &s, const Scalar &val)
    {
        assert(s > 0);

        //NEW SIZE
        const SizeType local_s = utopia::decompose(comm_, s);

        rcp_map_type map;

        UTOPIA_REPORT_ALLOC("TpetraVector::values");
        map = Teuchos::rcp(new map_type(s, local_s, 0, comm().get()));
        vec_.reset(new vector_type(map));
        implementation().putScalar(val);

        assert(size() > 0);
    }

    void TpetraVector::local_values(const SizeType &s, const Scalar &val)
    {
        assert(s > 0);

        rcp_map_type map;

        UTOPIA_REPORT_ALLOC("TpetraVector::values");
        map = Teuchos::rcp(new map_type(Teuchos::OrdinalTraits<Tpetra::global_size_t>::invalid(), s, 0, comm().get()));
        vec_.reset(new vector_type(map));
        implementation().putScalar(val);

        assert(size() > 0);
    }

    void TpetraVector::c_set(const SizeType &i, const Scalar &value)
    {
        implementation().replaceGlobalValue(i, value);
    }

    void TpetraVector::c_add(const SizeType &i, const Scalar &value)
    {
        if(!ghosted_vec_.is_null()) {
            ghosted_vec_->sumIntoGlobalValue(i, value);
        } else {
            implementation().sumIntoGlobalValue(i, value);
        }
    }

}
