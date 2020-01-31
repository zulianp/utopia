#include "utopia_petsc_Vector_impl.hpp"
#include "utopia_petsc_quirks.hpp"


#include <set>
#include <cstring>
#include <map>

namespace utopia {

 

    void PetscVector::transform(const Sqrt &op)
    {
        op_transform(op);
    }

    void PetscVector::transform(const Pow2 &op)
    {
        op_transform(op);
    }

    void PetscVector::transform(const Log &op)
    {
        op_transform(op);
    }

    void PetscVector::transform(const Exp &op)
    {
        op_transform(op);
    }

    void PetscVector::transform(const Cos &op)
    {
        op_transform(op);
    }

    void PetscVector::transform(const Sin &op)
    {
        op_transform(op);
    }

    void PetscVector::transform(const Abs &op)
    {
        op_transform(op);
    }

    void PetscVector::transform(const Minus &op)
    {
        op_transform(op);
    }


    void PetscVector::transform(const Pow &op)
    {
        op_transform(op);
    }

    void PetscVector::transform(const Reciprocal<Scalar> &op)
    {
        op_transform(op);
    }

    bool PetscVector::has_type(VecType type) const
    {
        PetscBool match = PETSC_FALSE;
        PetscObjectTypeCompare((PetscObject) implementation(), type, &match);
        return match == PETSC_TRUE;
    }

    bool PetscVector::same_type(const PetscVector &other) const
    {
        return has_type(other.type());
    }

    PetscScalar PetscVector::dot(const PetscVector &other) const
    {
        assert(is_consistent());
        assert(other.is_consistent());

        if(!same_type(other)) {
            //other.describe();
            // m_utopia_warning_once("> PetscVector::dot - vectors are not the same");
        }

        PetscScalar result;
        check_error( VecDot(implementation(), other.implementation(), &result) );
        return result;
    }

    void PetscVector::repurpose(MPI_Comm comm,
        VecType type,
        PetscInt n_local,
        PetscInt n_global)
    {

        assert(!immutable_);


        const std::string type_copy = type;

        if(is_null()) {
            UTOPIA_REPORT_ALLOC("PetscVector::repurpose");
            VecCreate(comm, &vec_);
        } else {
            if(comm != PetscObjectComm((PetscObject)vec_) || has_type(type)) {
                destroy();

                UTOPIA_REPORT_ALLOC("PetscVector::repurpose");
                check_error( VecCreate(comm, &vec_) );
            } else {
                PetscInt old_n_global;
                check_error( VecGetSize(vec_, &old_n_global) );

                if(old_n_global == n_global) {
                    PetscInt old_n_local;
                    check_error( VecGetLocalSize(vec_, &old_n_local) );

                    assert(old_n_local == n_local && "We do not handle the local consistency. Explicitly set local sizes in the initialization.");

                    return;
                }
            }
        }

        check_error( VecSetFromOptions(vec_) );

        check_error( VecSetType(vec_, type_copy.c_str()) );

        check_error( VecSetSizes(vec_, n_local, n_global) );

        ghost_values_.clear();

        assert(vec_ != nullptr);
        initialized_ = true;

        if(!is_consistent()) {
            std::cout << "type copy: " << type_copy << " != " << type_override() << "!=" << this->type() << std::endl;
        }

        assert(is_consistent());
    }

    void PetscVector::init(MPI_Comm comm,
        VecType type,
        PetscInt n_local,
        PetscInt n_global)
    {
        assert(vec_ == nullptr);

        UTOPIA_REPORT_ALLOC("PetscVector::repurpose");
        check_error( VecCreate(comm, &vec_) );
        check_error( VecSetFromOptions(vec_) );
        check_error( VecSetType(vec_, type) );

        check_error( VecSetSizes(vec_, n_local, n_global) );

        assert(vec_ != nullptr);
        initialized_ = true;

        assert(is_consistent());
    }

    void PetscVector::nest(
        MPI_Comm comm,
        PetscInt nb,
        IS is[],
        Vec x[],
        const bool use_vec_nest_type)
    {
        if(use_vec_nest_type) {
            destroy();
            UTOPIA_REPORT_ALLOC("PetscVector::nest");
            VecCreateNest(comm, nb, is, x, &vec_);
        } else {

            PetscInt ls = 0, gs = 0;
            PetscInt ls_i = 0, gs_i = 0;


            for(PetscInt i = 0; i < nb; ++i) {
                VecGetLocalSize(x[i], &ls_i);
                ls += ls_i;

                VecGetSize(x[i], &gs_i);
                gs += gs_i;
            }

            repurpose(comm, type_override(), ls, gs);

            // write_lock(LOCAL);

            auto r = range();

            const PetscScalar *a;
            PetscInt local_index = 0;
            for(PetscInt i = 0; i < nb; ++i) {
                VecGetLocalSize(x[i], &ls_i);
                VecGetArrayRead(x[i], &a);

                for(PetscInt k = 0; k < ls_i; ++k) {
                    VecSetValue(implementation(), r.begin() + local_index++, a[k], INSERT_VALUES);
                }

                VecRestoreArrayRead(x[i], &a);
            }

            // write_unlock(LOCAL);

            VecAssemblyBegin(implementation());
            VecAssemblyEnd(implementation());
        }
    }


    void PetscVector::copy_data_from(Vec vec)
    {
        MPI_Comm comm = PetscObjectComm((PetscObject) vec);

        PetscInt n_global;
        check_error( VecGetSize(vec, &n_global) );

        PetscInt n_local;
        check_error( VecGetLocalSize(vec, &n_local) );

        this->repurpose(
            comm,
            type_override(),
            n_local,
            n_global
        );

        VecCopy(vec, raw_type());

        // assert(vec_ != nullptr);
        set_initialized(true);
        // assert(is_consistent());
    }

    void PetscVector::copy_data_to(Vec vec) const
    {
        VecCopy(raw_type(), vec);
    }    

    void PetscVector::ghosted(MPI_Comm comm,
        PetscInt local_size,
        PetscInt global_size,
        const std::vector<PetscInt> &index)
    {

        assert(!immutable_);

        destroy();

        UTOPIA_REPORT_ALLOC("PetscVector::ghosted");
        check_error(
            VecCreateGhost(
                comm,
                local_size,
                global_size,
                static_cast<PetscInt>(index.size()),
                &index[0],
                &vec_)
            );

        //FIXME will this work???
        check_error( VecSetFromOptions(vec_) );

        init_ghost_index(index);
        check_error( VecZeroEntries(vec_) );

        assert(vec_ != nullptr);
        initialized_ = true;
    }


    bool PetscVector::is_mpi() const
    {
        static const std::string seq = "seq";
        const std::string str = type();
        return std::search(begin(str), end(str), begin(seq), end(seq)) == str.end();
    }

    bool PetscVector::has_nan_or_inf() const
    {
        PetscInt m;
        const PetscScalar *x;
        VecGetLocalSize(implementation(), &m);
        VecGetArrayRead(implementation(), &x);

        int has_nan = 0;

        for (PetscInt i = 0; i < m; i++) {
            has_nan = PetscIsInfOrNanScalar(x[i]);
            if(has_nan == 1)
                break;
        }

        VecRestoreArrayRead(implementation(), &x);
        MPI_Comm comm = PetscObjectComm((PetscObject) implementation());

        if(is_mpi()) {
            MPI_Allreduce(MPI_IN_PLACE, &has_nan, 1, MPI_INT, MPI_MAX, comm);
        }

        return has_nan > 0;
    }

    void PetscVector::resize(PetscInt local_size, PetscInt global_size)
    {
        assert(!is_null());

        // if(initialized()) {
        // 	VecSetSizes(implementation(), local_size, global_size);
        // } else {
        MPI_Comm comm = communicator();
        const std::string type = this->type();

        destroy();

        init(comm, type.c_str(), local_size, global_size);
        // }
    }

    void PetscVector::select(
        const PetscIndexSet &index,
        PetscVector &result) const
    {
        MPI_Comm comm = communicator();
        IS is_in;
        VecScatter scatter_context;

        result.repurpose(comm, this->type(), index.size(), PETSC_DETERMINE);

        ISCreateGeneral(comm, index.size(), &index[0], PETSC_USE_POINTER, &is_in);

        UtopiaVecScatterCreate(implementation(), is_in, result.implementation(), nullptr, &scatter_context);

        VecScatterBegin(scatter_context, implementation(), result.implementation(), INSERT_VALUES, SCATTER_FORWARD);
        VecScatterEnd(scatter_context,   implementation(), result.implementation(), INSERT_VALUES, SCATTER_FORWARD);

        ISDestroy(&is_in);
        VecScatterDestroy(&scatter_context);
    }

    void PetscVector::copy_from(Vec vec)
    {
        destroy();

        UTOPIA_REPORT_ALLOC("PetscVector::copy_from");
        VecDuplicate(vec, &implementation());
        VecCopy(vec, implementation());

        assert(vec_ != nullptr);
        set_initialized(true);
        assert(is_consistent());
    }

    void PetscVector::convert_from(const Vec &vec)
    {
        copy_from(vec);
    }

    void PetscVector::convert_to(Vec &vec) const
    {
        check_error( VecCopy(raw_type(), vec) );
    }

    bool PetscVector::read(MPI_Comm comm, const std::string &path)
    {
        destroy();

        PetscViewer fd;

        bool err = check_error( PetscViewerBinaryOpen(comm, path.c_str(), FILE_MODE_READ, &fd) );

        UTOPIA_REPORT_ALLOC("PetscVector::read");
        err = err && check_error( VecCreate(comm, &implementation()) );
        err = err && check_error( VecSetType(implementation(), type_override()) );
        err = err && check_error( VecLoad(implementation(), fd));

        set_initialized(true);
        assert(is_consistent());

        PetscViewerDestroy(&fd);
        return err;
    }

    bool PetscVector::write(const std::string &path) const
    {
        if(is_matlab_file(path)) {
            return write_matlab(path);
        } else {
            return write_binary(path);
        }
    }
    
    bool PetscVector::write_binary(const std::string &path) const
    {
        PetscViewer fd;
        bool err = check_error( PetscViewerBinaryOpen(communicator(), path.c_str(), FILE_MODE_WRITE, &fd) );

        err = err && check_error( VecView(implementation(), fd) );

        PetscViewerDestroy(&fd);
        return err;
    }

    bool PetscVector::write_matlab(const std::string &path) const
    {
        PetscViewer fd;
        bool err = check_error( PetscViewerASCIIOpen(communicator(), path.c_str(), &fd) );

        err = err && check_error( PetscViewerPushFormat(fd, PETSC_VIEWER_ASCII_MATLAB) );
        err = err && check_error( VecView(implementation(), fd) );

        PetscViewerDestroy(&fd);
        return err;
    }

    void PetscVector::select(const Range &global_range, PetscVector &result) const
    {
        assert(!global_range.empty());

        Range rr = range().intersect(global_range);

        result.repurpose(
            communicator(),
            type(),
            PETSC_DECIDE,
            rr.extent()
            );

        result.write_lock(LOCAL);

        for(PetscInt r_this = rr.begin(); r_this < rr.end(); ++r_this) {
            const PetscInt r_selection = r_this - global_range.begin();
            result.set(r_selection, get(r_this));
        }

        result.write_unlock(LOCAL);
    }

    //testing VECSEQCUDA,VECMPICUDA
    bool PetscVector::is_cuda() const
    {
        PetscBool match = PETSC_FALSE;
        PetscObjectTypeCompare((PetscObject) implementation(), VECSEQCUDA, &match);
        if(match == PETSC_TRUE) return true;

        PetscObjectTypeCompare((PetscObject) implementation(), VECMPICUDA, &match);
        return match == PETSC_TRUE;
    }

    void PetscVector::describe() const {

        if(is_root()) {
            std::cout << "is_null    : " << is_null() << "\n";
            std::cout << "initialized: " << initialized() << "\n";
        }

        if(is_null()) return;

        VecView(implementation(), PETSC_VIEWER_STDOUT_(communicator()));
    }

    bool PetscVector::is_root() const
    {
        auto comm = communicator();
        int rank;
        MPI_Comm_rank(comm, &rank);
        return rank == 0;
    }

    bool PetscVector::is_consistent() const
    {
        if(initialized_ && is_null()) {
            return false;
        }

        if(is_null()) {
            return true;
        }

        if(!initialized()) {
            return true;
        }

        //TODO add additional checks
        return true;
    }

    void PetscVector::write_lock(WriteMode mode)
    {
        switch(mode) {
            case GLOBAL_INSERT:
            case GLOBAL_ADD:
            {
                //no-op
                break;
            }
            case LOCAL:
            case AUTO:
            default:
            {
                writeable_ = utopia::make_unique<LocalView>(implementation());
                readable_  = utopia::make_unique<ConstLocalView>(*writeable_);
                break;
            }
        }
    }

    void PetscVector::write_unlock(WriteMode mode)
    {
        switch(mode) {
            case GLOBAL_INSERT:
            case GLOBAL_ADD:
            {
                VecAssemblyBegin(implementation());
                VecAssemblyEnd(implementation());

                set_initialized(true);
                update_ghosts();
                break;
            }
            case LOCAL:
            case AUTO:
            default:
            {
                writeable_ = nullptr;
                readable_  = nullptr;

                if(!initialized_) {
                    VecAssemblyBegin(implementation());
                    VecAssemblyEnd(implementation());

                    set_initialized(true);
                    update_ghosts();
                }

                break;
            }
        }
    }

    void PetscVector::read_and_write_lock(WriteMode mode) { write_lock(mode); }
    void PetscVector::read_and_write_unlock(WriteMode mode) { write_unlock(mode); }

    bool PetscVector::equals(const PetscVector &other, const Scalar &tol) const
    {
        if(this->is_alias(other)) return true;
        
        PetscVector diff = other;
        diff.axpy(-1.0, *this);
        return diff.norm_infty() < tol;
    }

    void PetscVector::e_mul(const PetscVector &other)
    {
         assert(is_consistent());
         assert(other.is_consistent());
         check_error( VecPointwiseMult( implementation(), implementation(), other.implementation()) );
    }

    void PetscVector::e_div(const PetscVector &other)
    {
        assert(is_consistent());
        assert(other.is_consistent());
        check_error( VecPointwiseDivide(raw_type(), raw_type(), other.raw_type() ) );
    }
   
    void PetscVector::e_min(const PetscVector &other)
    {
        check_error( VecPointwiseMin( raw_type(), raw_type(), other.raw_type() ) );
    }

    void PetscVector::e_max(const PetscVector &other)
    {
        check_error( VecPointwiseMax( raw_type(), raw_type(), other.raw_type() ) );
    }

    void PetscVector::e_mul(const Scalar &other)
    {
        element_wise_generic(other, Multiplies(), *this);
    }

    void PetscVector::e_div(const Scalar &other)
    {
        element_wise_generic(other, Divides(), *this);
    }

    void PetscVector::e_min(const Scalar &other)
    {
        element_wise_generic(other, Min(), *this);
    }

    void PetscVector::e_max(const Scalar &other)
    {
        element_wise_generic(other, Max(), *this);
    }

    ///<Scalar>SWAP - swap x and y
    void PetscVector::swap(PetscVector &x) {
         using std::swap;
         swap(comm_, x.comm_);
         swap(vec_, x.vec_);
         swap(initialized_, x.initialized_);
         swap(ghost_values_, x.ghost_values_);
         swap(immutable_, x.immutable_);
    }

    ///<Scalar>SCAL - x = a*x
    void PetscVector::scale(const Scalar &a)
    {
         check_error( VecScale(implementation(), a) );
    }

    ///<Scalar>COPY - copy other into this
     void PetscVector::copy(const PetscVector &other)
    {
         if(this == &other) return;

         assert(!immutable_);

         if(is_compatible(other) && !other.has_ghosts()) {

            // assert(same_type(other) && "TYPE " );
            // assert(this->has_ghosts() && "GHOST" );                

             assert((same_type(other) || this->has_ghosts()) && "Inconsistent vector types. Handle types properly before copying" );
             assert(local_size() == other.local_size() && "Inconsistent local sizes. Handle local sizes properly before copying.");
             PetscErrorHandler::Check(VecCopy(other.vec_, vec_));
             initialized_ = other.initialized_;
             immutable_ = other.immutable_;
             return;
         }

         destroy();

         if(other.vec_) {
             UTOPIA_REPORT_ALLOC("PetscVector::copy");
             PetscErrorHandler::Check(VecDuplicate(other.vec_, &vec_));
             PetscErrorHandler::Check(VecCopy(other.vec_, vec_));
             ghost_values_ = other.ghost_values_;

             initialized_ = other.initialized_;
             immutable_ = other.immutable_;
         } else {
             initialized_ = false;
         }

         return;
    }

    ///<Scalar>AXPY - y = a*x + y
    void PetscVector::axpy(const Scalar &alpha, const PetscVector &x)
    {
         assert(is_consistent());
         assert(x.is_consistent());

         check_error( VecAXPY(implementation(), alpha, x.implementation()) );
    }

    PetscVector::SizeType PetscVector::amax() const
    {
        SizeType idx = 0;
        Scalar val = 0.0;
        VecMax(raw_type(), &idx, &val);
        return idx;
    }
}
