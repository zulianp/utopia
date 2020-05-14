#ifndef UTOPIA_PETSC_DM_HPP
#define UTOPIA_PETSC_DM_HPP

#include "utopia_Input.hpp"
#include "utopia_make_unique.hpp"
#include "utopia_Tracer.hpp"
#include "utopia_petsc_Communicator.hpp"
#include "utopia_Path.hpp"
#include "utopia_petsc_Vector.hpp"
#include "utopia_petsc_Matrix.hpp"
#include "utopia_petsc_IO.hpp"

#include <petscdm.h>

namespace utopia {
    class PetscDMBase : public Configurable {
    public:
        ~PetscDMBase() override {}

        template<class SizeType>
        static void dof_ownership_range(const DM &dm, SizeType &begin, SizeType &end)
        {
            Vec v;
            DMGetGlobalVector(dm, &v);
            VecGetOwnershipRange(v, &begin, &end);
            DMRestoreGlobalVector(dm, &v);
        }

        static PetscInt get_dimension(DM dm)
        {
            PetscInt ret;
            DMGetDimension(dm, &ret);
            return ret;
        }

        static void refine(const DM &in, const MPI_Comm &comm, DM &out)
        {
            DMRefine(in, comm, &out);
        }

        DM &raw_type()
        {
            return wrapper_->dm;
        }

        const DM &raw_type() const
        {
            return wrapper_->dm;
        }

        virtual void wrap(DM &dm, const bool delegate_ownership)
        {
            destroy_dm();
            wrapper_->dm = dm;
            wrapper_->owned = delegate_ownership;
        }

        void destroy_dm()
        {
            wrapper_->destroy();
        }

        PetscDMBase(const PetscCommunicator &comm = PetscCommunicator())
        : comm_(comm), wrapper_(utopia::make_unique<Wrapper>())
        {}

        inline PetscCommunicator &comm()
        {
            return comm_;
        }

        inline const PetscCommunicator &comm() const
        {
            return comm_;
        }

        void create_matrix(PetscMatrix &mat) const
        {
            UTOPIA_TRACE_REGION_BEGIN("PetscDMbase::create_matrix(...)");

            mat.destroy();
            DMCreateMatrix(raw_type(), &mat.raw_type());
            UTOPIA_REPORT_ALLOC("PetscDMbase::create_matrix");

            UTOPIA_TRACE_REGION_END("PetscDMbase::create_matrix(...)");
        }

        void create_vector(PetscVector &vec) const
        {
            UTOPIA_TRACE_REGION_BEGIN("PetscDMbase::create_vector(...)");

            vec.destroy();
            DMCreateGlobalVector(raw_type(), &vec.raw_type());
            UTOPIA_REPORT_ALLOC("PetscDMbase::create_vector");

            UTOPIA_TRACE_REGION_END("PetscDMbase::create_vector(...)");
        }

        void create_local_vector(PetscVector &vec) const
        {
            UTOPIA_TRACE_REGION_BEGIN("PetscDMbase::create_local_vector(...)");

            vec.destroy();
            auto err = DMCreateLocalVector(raw_type(), &vec.raw_type()); assert(err == 0);
            UTOPIA_REPORT_ALLOC("PetscDMbase::create_local_vector");

            UTOPIA_TRACE_REGION_END("PetscDMbase::create_local_vector(...)");
        }

        void create_interpolation(const PetscDMBase &target, PetscMatrix &I) const
        {
            UTOPIA_TRACE_REGION_BEGIN("PetscDMbase::create_interpolation(...)");

            I.destroy();
            auto ierr = DMCreateInterpolation(raw_type(), target.raw_type(), &I.raw_type(), nullptr); assert(ierr == 0);
            UTOPIA_REPORT_ALLOC("PetscDMbase::create_local_vector");

            UTOPIA_TRACE_REGION_END("PetscDMbase::create_interpolation(...)");
        }

        void local_to_global(const PetscVector &local,  PetscVector &global) const
        {
            UTOPIA_TRACE_REGION_BEGIN("PetscDMbase::local_to_global(...)");

#if UTOPIA_PETSC_VERSION_GREATER_EQUAL_THAN(3, 11, 0) //DM-INCOMPLETE
            DMLocalToGlobal(raw_type(), local.raw_type(), ADD_VALUES, global.raw_type());
#else
            DMLocalToGlobalBegin(raw_type(), local.raw_type(), ADD_VALUES, global.raw_type());
            DMLocalToGlobalEnd(raw_type(), local.raw_type(), ADD_VALUES, global.raw_type());
#endif

            UTOPIA_TRACE_REGION_END("PetscDMbase::local_to_global(...)");
        }

        void global_to_local(const PetscVector &global, PetscVector &local) const
        {
            UTOPIA_TRACE_REGION_BEGIN("PetscDMbase::global_to_local(...)");

#if UTOPIA_PETSC_VERSION_GREATER_EQUAL_THAN(3, 11, 0) //DM-INCOMPLETE
            DMGlobalToLocal(raw_type(), global.raw_type(), INSERT_VALUES, local.raw_type());
#else
            DMGlobalToLocalBegin(raw_type(), global.raw_type(), INSERT_VALUES, local.raw_type());
            DMGlobalToLocalEnd(raw_type(), global.raw_type(), INSERT_VALUES, local.raw_type());
#endif

            UTOPIA_TRACE_REGION_END("PetscDMbase::global_to_local(...)");
        }

        bool write(const Path &path, const PetscVector &x) const
        {
            // UTOPIA_TRACE_REGION_BEGIN("PetscDMbase::write(...)");
            PetscIO io;
            if(!io.open(comm(), path)) { return false; }
            if(!io.write(*this))       { return false; }
            if(!io.write(x))           { return false; }

            // UTOPIA_TRACE_REGION_END("PetscDMbase::write(...)");
            return true;
        }

        bool write(const Path &path) const
        {
            PetscIO io;
            if(!io.open(comm(), path)) { return false; }
            return io.write(*this);
        }

        virtual void describe() const
        {

        }

        inline PetscInt dm_dim() const
        {
            PetscInt ret;
            DMGetDimension(raw_type(), &ret);
            return ret;
        }

    private:
        PetscCommunicator comm_;

        class Wrapper {
        public:
            Wrapper()
            : dm(nullptr), owned(true)
            {}

            ~Wrapper() {
                destroy();
            }

            void destroy()
            {
                if(owned && dm) {
                    DMDestroy(&dm);
                    dm = nullptr;
                }
            }

            DM dm;
            bool owned;
        };

        std::unique_ptr<Wrapper> wrapper_;
    };
}

#endif //UTOPIA_PETSC_DM_HPP
