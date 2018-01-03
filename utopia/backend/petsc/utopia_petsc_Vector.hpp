
#ifndef UTOPIA_UTOPIA_PETSCVECTOR_H
#define UTOPIA_UTOPIA_PETSCVECTOR_H

#include "utopia_petsc_Error.hpp"

#include "utopia_Range.hpp"
#include "utopia_Base.hpp"

// #include "petscmat.h"
#include "petscvec.h"

#include <memory>
#include <map>
#include <vector>

namespace utopia {

    class PETScVector {
    private:
        class GhostValues {
        public:
            GhostValues()
            : has_ghosts_(false)
            {}

            void update(Vec &vec)
            {
                if(!has_ghosts()) return;
               
                VecGhostUpdateBegin(vec, INSERT_VALUES, SCATTER_FORWARD);
                VecGhostUpdateEnd(vec,   INSERT_VALUES, SCATTER_FORWARD);
            }

            inline bool has_ghosts() const
            {
                return has_ghosts_;
            }

            inline void init_index(
                const PetscInt n_local,
                const std::vector<PetscInt> &index)
            {
                has_ghosts_ = true;
                const PetscInt n = index.size();
               
                for(PetscInt i = 0; i < n; ++i) {
                    ghost_index_[index[i]] = n_local + i;
                }
            }

            inline PetscInt get_index(const PetscInt &g_index) const
            {
                auto it = ghost_index_.find(g_index);
                
                if(it == ghost_index_.end()) {
                    std::cerr << "[Error] index not present in ghosted vector" << std::endl;
                    assert(false);
                    return -1;
                }

                return it->second;
            }

            inline void clear()
            {
                has_ghosts_ = false;
                ghost_index_.clear();
            }

            std::map<PetscInt, PetscInt> ghost_index_;
            bool has_ghosts_;
        };

    public:


        inline PETScVector()
        : vec_(nullptr), initialized_(false)
        {
#ifndef NDEBUG
            immutable_ = false;
#endif            
        }

        inline ~PETScVector()
        {
            destroy();
        }

        PETScVector(const PETScVector &other)
        {
            if(other.vec_) {
                PETScError::Check(VecDuplicate(other.vec_, &vec_));
                PETScError::Check(VecCopy(other.vec_, vec_));
                initialized_ = other.initialized_;
                ghost_values_ = other.ghost_values_;
            } else {
                vec_ = nullptr;
                initialized_ = false;
            }

#ifndef NDEBUG
            immutable_ = other.immutable_;
#endif 
        }

        inline VecType type() const
        {
            assert(vec_ != nullptr);
            VecType ret;
            VecGetType(vec_, &ret);
            return ret;
        }

        inline PetscInt local_size() const
        {
            PetscInt ret;
            VecGetLocalSize(vec_, &ret);
            return ret;
        }

        inline PetscInt size() const
        {
            PetscInt ret;
            VecGetSize(vec_, &ret);
            return ret;
        }

        inline MPI_Comm communicator() const {
            assert(vec_ != nullptr);
            MPI_Comm comm = PetscObjectComm((PetscObject) implementation());
            assert(comm != MPI_COMM_NULL);
            return comm;
        }

        // assign operator
       inline PETScVector &operator=(const PETScVector &other) {
            if(this == &other) return *this;
            assert(!immutable_);

            initialized_ = other.initialized_;

#ifndef NDEBUG
            immutable_ = other.immutable_;
#endif 

            if(!is_null() && size() == other.size() && !other.has_ghosts()) {
                assert(local_size() == other.local_size() && "Inconsistent local sizes. Handle local sizes properly before copying.");
                PETScError::Check(VecCopy(other.vec_, vec_));
                return *this;
            }

            destroy();

            if(other.vec_) {
                PETScError::Check(VecDuplicate(other.vec_, &vec_));
                PETScError::Check(VecCopy(other.vec_, vec_));
                ghost_values_ = other.ghost_values_;
            }

            return *this;
        }

        inline PETScVector &operator=(PETScVector &&other) {
            if(this == &other) return *this;
            assert(!immutable_);

            destroy();

            initialized_ = other.initialized_;

#ifndef NDEBUG
            immutable_ = other.immutable_;
#endif  
            vec_ = other.vec_;
            other.vec_ = nullptr;
            ghost_values_ = std::move(other.ghost_values_);
            other.initialized_ = false;

#ifndef NDEBUG
            other.immutable_ = false;
#endif 
            return *this;
        }

        inline void destroy() {
            if(vec_) {
                VecDestroy(&vec_);
                vec_ = nullptr;
            }

            initialized_ = false;
            ghost_values_.clear();
        }

        inline Vec &implementation() {
            return vec_;
        }

        inline const Vec &implementation() const {
            assert(vec_ != nullptr);
            return vec_;
        }

        inline void describe() const {
            assert(vec_ != nullptr);
            VecView(vec_, PETSC_VIEWER_STDOUT_(communicator()));
        }

        inline bool is_null() const
        {
            return vec_ == nullptr;
        }

        inline bool initialized() const {
            return initialized_;
        }

        inline void set_initialized(const bool val) 
        {
            initialized_ = val;
        }

        inline PetscInt global_to_local(const PetscInt g_index) const
        {
            PetscInt begin, end;
            VecGetOwnershipRange(vec_, &begin, &end);

            if(!ghost_values_.has_ghosts() || (g_index >= begin && g_index < end)) {
                return g_index - begin;
            }

            return ghost_values_.get_index(g_index);
        }

        inline void init_ghost_index(const std::vector<PetscInt> &index)
        {
            ghost_values_.init_index(local_size(), index);
        }

        inline void get_values(
            const std::vector<PetscInt> &index,
            std::vector<PetscScalar> &values) const
        {
            std::size_t n = index.size();

            values.resize(n);

            if(!ghost_values_.has_ghosts()) {

                VecGetValues(
                    vec_,
                    static_cast<PetscInt>(index.size()),
                    &index[0],
                    &values[0]
                );

            } else {

                const PetscScalar *array;
                Vec local_form;

                VecGhostGetLocalForm(vec_, &local_form);
                VecGetArrayRead(local_form, &array);

                assert(local_form != nullptr);

                for(std::size_t i = 0; i < n; ++i) {
                    auto li = global_to_local(index[i]);
                    values[i] = array[li];
                }

                VecRestoreArrayRead(local_form, &array);
                VecGhostRestoreLocalForm(vec_, &local_form);
            }
        }

        inline bool has_ghosts() const
        {
            return !ghost_values_.has_ghosts();
        }

        inline void update_ghosts()
        {
            ghost_values_.update(vec_);
        }
        
#ifndef NDEBUG
        inline void make_immutable()
        {
            immutable_ = true;
        }
#endif

        //builders
        void repurpose(MPI_Comm comm, VecType type, PetscInt n_local, PetscInt n_global);
        inline void zeros(MPI_Comm comm, VecType type, PetscInt n_local, PetscInt n_global)
        {
            repurpose(comm, type, n_local, n_global);
            VecZeroEntries(vec_);
        }

        void init(MPI_Comm comm, VecType type, PetscInt n_local, PetscInt n_global);
        void ghosted(MPI_Comm comm, PetscInt local_size, PetscInt global_size, const std::vector<PetscInt> &index);

        //ops
        ///this is y
       inline void axpy(const PetscScalar &alpha, const PETScVector &x)
       {
         check_error( VecAXPY(implementation(), alpha, x.implementation()) );
       }

       ///this is y
       inline void axpby(const Scalar alpha, const Vector &x, const Scalar &beta)
       {
         check_error( VecAXPBY(vec_, alpha, beta, x.implementation()) );
       }

    private:
        Vec vec_;
        bool initialized_;

        GhostValues ghost_values_;

        //debug
#ifndef NDEBUG
        bool immutable_;
#endif //NDEBUG     

        inline bool check_error(const PetscInt err) {
            return PETScError::Check(err);
        }   
    };

}

#endif //UTOPIA_UTOPIA_PETSCVECTOR_H
