
#ifndef UTOPIA_UTOPIA_PETSCVECTOR_H
#define UTOPIA_UTOPIA_PETSCVECTOR_H

#include "utopia_petsc_Error.hpp"

#include "utopia_Range.hpp"
#include "utopia_Base.hpp"

#include "petscmat.h"
#include "petscvec.h"

#include <memory>
#include <map>
#include <vector>

namespace utopia {

    class PETScVector {
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
                ghost_index = other.ghost_index;
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

            destroy();

            if(other.vec_) {
                PETScError::Check(VecDuplicate(other.vec_, &vec_));
                PETScError::Check(VecCopy(other.vec_, vec_));

                ghost_index = other.ghost_index;
            }

#ifndef NDEBUG
            immutable_ = other.immutable_;
#endif 

            return *this;
        }

        inline PETScVector &operator=(PETScVector &&other) {
            if(this == &other) return *this;
            assert(!immutable_);

            destroy();

#ifndef NDEBUG
            immutable_ = other.immutable_;
#endif  

            vec_ = other.vec_;
            other.vec_ = nullptr;
            ghost_index = std::move(other.ghost_index);

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
            ghost_index.clear();
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

            if(ghost_index.empty() || (g_index >= begin && g_index < end)) {
                return g_index - begin;
            }

            auto it = ghost_index.find(g_index);
            
            if(it == ghost_index.end()) {
                std::cerr << "[Error] index not present in ghosted vector" << std::endl;
                return begin;
            }

            return it->second;
        }

        inline void init_ghost_index(const std::vector<PetscInt> &index)
        {
            const std::size_t n = index.size();
            const std::size_t n_local = local_size();
            for(std::size_t i = 0; i < n; ++i) {
                ghost_index[index[i]] = n_local+i;
            }
        }

        inline void get_values(
            const std::vector<PetscInt> &index,
            std::vector<PetscScalar> &values) const
        {
            std::size_t n = index.size();

            values.resize(n);

            if(ghost_index.empty()) {

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

        inline void update_ghosts()
        {
            if(ghost_index.empty()) return;
            
            VecGhostUpdateBegin(vec_, INSERT_VALUES, SCATTER_FORWARD);
            VecGhostUpdateEnd(vec_,   INSERT_VALUES, SCATTER_FORWARD);
        }
        
#ifndef NDEBUG
        inline void make_immutable()
        {
            immutable_ = true;
        }
#endif

    private:
        Vec vec_;
        bool initialized_;

        //ghosted
        std::map<PetscInt, PetscInt> ghost_index;


        //debug
#ifndef NDEBUG
        bool immutable_;
#endif //NDEBUG        
    };

}

#endif //UTOPIA_UTOPIA_PETSCVECTOR_H
