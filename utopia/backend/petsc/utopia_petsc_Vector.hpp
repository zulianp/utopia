#ifndef UTOPIA_UTOPIA_PETSCVECTOR_H
#define UTOPIA_UTOPIA_PETSCVECTOR_H

#include "utopia_Base.hpp"
#include "utopia_Range.hpp"
#include "utopia_make_unique.hpp"
#include "utopia_Vector.hpp"
#include "utopia_Tensor.hpp"
#include "utopia_Normed.hpp"
#include "utopia_Reducible.hpp"
#include "utopia_Transformable.hpp"
#include "utopia_ElementWiseOperand.hpp"
#include "utopia_Constructible.hpp"
#include "utopia_Comparable.hpp"
#include "utopia_BLAS_Operands.hpp"
#include "utopia_Allocations.hpp"
#include "utopia_Select.hpp"

#include "utopia_petsc_Base.hpp"
#include "utopia_petsc_ForwardDeclarations.hpp"
#include "utopia_petsc_Error.hpp"
#include "utopia_petsc_Communicator.hpp"
#include "utopia_petsc_IndexSet.hpp"
#include "utopia_petsc_Traits.hpp"

#include <map>
#include <vector>
#include <limits>

#include "petscvec.h"

namespace utopia {

    class PetscVector : 
        public DistributedVector<PetscScalar, PetscInt>,
        public Normed<PetscScalar>,
        public Transformable<PetscScalar>,
        public Reducible<PetscScalar>,
        public Constructible<PetscScalar, PetscInt, 1>,
        public ElementWiseOperand<PetscScalar>,
        public ElementWiseOperand<PetscVector>,
        public Comparable<PetscVector>,
        // public Ranged<PetscVector, 1>,
        public BLAS1Tensor<PetscVector>,
        public Tensor<PetscVector, 1>,
        public Selectable<PetscVector, 1>
    {
    public:
            using Scalar   = PetscScalar;
            using SizeType = PetscInt;
            using Super    = utopia::Tensor<PetscVector, 1>;
            using Constructible = utopia::Constructible<PetscScalar, PetscInt, 1>;

            using Super::Super;
            using Constructible::values;
            using Constructible::local_values;
            using Constructible::local_zeros;
            using Constructible::zeros;

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
                const SizeType n_local,
                const std::vector<SizeType> &index)
            {
                has_ghosts_ = true;
                const SizeType n = index.size();

                for(SizeType i = 0; i < n; ++i) {
                    ghost_index_[index[i]] = n_local + i;
                }
            }

            inline SizeType get_index(const SizeType &g_index) const
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

            std::map<SizeType, SizeType> ghost_index_;
            bool has_ghosts_;
        };

    public:


        ////////////////////////////////////////////////////////////////////
        ///////////////////////// BOILERPLATE CODE FOR EDSL ////////////////
        ////////////////////////////////////////////////////////////////////
        void init_empty()
        {
          vec_ = nullptr; 
          initialized_ = false;
          immutable_ = false;
        }

         template<class Expr>
         PetscVector(const Expression<Expr> &expr)
         {
            init_empty();
             //THIS HAS TO BE HERE IN EVERY UTOPIA TENSOR CLASS
             Super::construct_eval(expr.derived());
         }

         template<class Expr>
         inline PetscVector &operator=(const Expression<Expr> &expr)
         {
             Super::assign_eval(expr.derived());
             return *this;
         }

         void assign(const PetscVector &other) override
         {
             copy(other);
         }

         void assign(PetscVector &&other) override
         {
             comm_ = std::move(other.comm_);
             vec_ = std::move(other.vec_);
             initialized_ = std::move(other.initialized_);
             ghost_values_ = std::move(other.ghost_values_);
             immutable_ = std::move(other.immutable_);
         }

         ////////////////////////////////////////////////////////////////////

        ///////////////////////////////////////////////////////////////////////////
        ////////////// OVERRIDES FOR DistributedObject ////////////////////////////
        ///////////////////////////////////////////////////////////////////////////

        PetscCommunicator &comm() override
        {
            return comm_;
        }

        const PetscCommunicator &comm() const override
        {
            return comm_;
        }

        ///////////////////////////////////////////////////////////////////////////
        ////////////// OVERRIDES FOR VectorBase ////////////////////////////
        ///////////////////////////////////////////////////////////////////////////


        //locks

        inline void read_lock() override {
            readable_  = utopia::make_unique<ConstLocalView>(implementation());
        }

        inline void read_unlock() override 
        {
            readable_ = nullptr;
        }

        void write_lock(WriteMode mode) override ;
        void write_unlock(WriteMode mode) override ;

        void read_and_write_lock(WriteMode mode) override;
        void read_and_write_unlock(WriteMode mode) override;

        //basic mutators
        inline void set(const SizeType &index, const Scalar &value) override
        {
            assert(range().inside(index));
            // check_error(
            // VecSetValues(implementation(), 1, &index, &value, INSERT_VALUES);

            assert((writeable_) && "use Write<Vector> w(vec, LOCAL) before using get. Check if you are using a copy of the vector");
            writeable_->data[index - writeable_->range_begin] = value;
                 // );
        }

        inline void add(const SizeType &index, const Scalar &value) override
        {
            assert(range().inside(index));
            // check_error(
            // VecSetValues(implementation(), 1, &index, &value, ADD_VALUES);
                // );

            assert((writeable_) && "use Write<Vector> w(vec, LOCAL) before using get. Check if you are using a copy of the vector");
            writeable_->data[index - writeable_->range_begin] += value;
        }

        inline Scalar get(const SizeType &index) const override
        {
            // Scalar value;
            // VecGetValues(implementation(), 1, &index, &value);
            // return value;
            assert(range().inside(index));
            assert((readable_) && "use Read<Vector> r(vec) before using get. Check if you are using a copy of the vector");

            return readable_->data[index - readable_->range_begin];
        }


        //print function
        void describe() const override;

        //utility functions
        inline bool empty() const override { return is_null() || size() == 0; }
        inline void clear() override { destroy(); }
        
        inline void set(const Scalar &value) override
        {
            check_error( VecSet(implementation(), value) );
        }

       inline SizeType size() const override
       {
           if(is_null()) {
               return utopia::INVALID_INDEX;
           }

           SizeType ret;
           VecGetSize(implementation(), &ret);
           return ret;
       }


       ///////////////////////////////////////////////////////////////////////////
       ////////////// OVERRIDES FOR DistributedVector ////////////////////////////
       ///////////////////////////////////////////////////////////////////////////


       inline void c_set(const SizeType &i, const Scalar &value) override
       {
            check_error( VecSetValues(implementation(), 1, &i, &value, INSERT_VALUES) );
       }

       inline void c_add(const SizeType &i, const Scalar &value) override
       {
            check_error( VecSetValues(implementation(), 1, &i, &value, ADD_VALUES) );
       }
       
       inline SizeType local_size() const override
       {
           if(is_null()) return 0;

           SizeType ret;
           VecGetLocalSize(implementation(), &ret);
           return ret;
       }

       ///////////////////////////////////////////////////////////////////////////
       ////////////// OVERRIDES FOR Normed ////////////////////////////
       ///////////////////////////////////////////////////////////////////////////

       inline Scalar norm2() const override {
           Scalar val;
           check_error( VecNorm(implementation(), NORM_2, &val) );
           return val;
       }

       inline Scalar norm1() const override {
           Scalar val;
           check_error( VecNorm(implementation(), NORM_1, &val) );
           return val;
       }

       inline Scalar norm_infty() const override {
           Scalar val;
           check_error( VecNorm(implementation(), NORM_INFINITY, &val) );
           return val;
       }

       ///////////////////////////////////////////////////////////////////////////
       ////////////// OVERRIDES FOR ElementWiseOperand ////////////////////////////
       ///////////////////////////////////////////////////////////////////////////

       void e_mul(const PetscVector &other) override;
       void e_div(const PetscVector &other) override;
       void e_min(const PetscVector &other) override;
       void e_max(const PetscVector &other) override;

       //against scalars
       void e_mul(const Scalar &other) override;
       void e_div(const Scalar &other) override;
       void e_min(const Scalar &other) override;
       void e_max(const Scalar &other) override;


       ///////////////////////////////////////////////////////////////////////////
       ////////////// OVERRIDES FOR Transformable ////////////////////////////
       ///////////////////////////////////////////////////////////////////////////

       void transform(const Sqrt &op) override;
       void transform(const Pow2 &op) override;
       void transform(const Log &op)  override;
       void transform(const Exp &op)  override;
       void transform(const Cos &op)  override;
       void transform(const Sin &op)  override;
       void transform(const Abs &op)  override;
       void transform(const Minus &op) override;

       void transform(const Pow &p)  override;
       void transform(const Reciprocal<Scalar> &f) override;


       template<class Operation>
       void aux_transform(const Operation &op);


       ///////////////////////////////////////////////////////////////////////////
       ////////////// OVERRIDES FOR BLAS1Operand ////////////////////////////
       ///////////////////////////////////////////////////////////////////////////

       ///<Scalar>SWAP - swap x and y
       void swap(PetscVector &x) override;
       ///<Scalar>SCAL - x = a*x
       void scale(const Scalar &a) override;
       ///<Scalar>COPY - copy other into this
        void copy(const PetscVector &other) override;
       ///<Scalar>AXPY - y = a*x + y
       void axpy(const Scalar &alpha, const PetscVector &x) override;

        SizeType amax() const; //override;

       ///<Scalar>DOT - dot product
       Scalar dot(const PetscVector &other) const override;
        
       ///////////////////////////////////////////////////////////////////////////
       ////////////// OVERRIDES FOR Reducible ////////////////////////////
       ///////////////////////////////////////////////////////////////////////////

       Scalar reduce(const Plus &) const override
       {
        Scalar result = 0;
        check_error( VecSum(implementation(), &result) );
        return result;
      }

      Scalar reduce(const Min &) const override
      {
        Scalar result = std::numeric_limits<Scalar>::max();
        check_error( VecMin(implementation(), nullptr, &result) );
        return result;
      }

      Scalar reduce(const Max &) const override
      {
        Scalar result = -std::numeric_limits<Scalar>::max();
        check_error( VecMax(implementation(), nullptr, &result) );
        return result;
      }

      ///////////////////////////////////////////////////////////////////////////
      ////////////// OVERRIDES FOR Constructible ////////////////////////////
      ///////////////////////////////////////////////////////////////////////////

      inline void zeros(const SizeType &s) override
      {
        zeros(
            comm().get(),
            type_override(),
            PETSC_DECIDE,
            s
        );
      }

      inline void values(const SizeType &s, const Scalar &val) override
      {
        values(
            comm().get(),
            type_override(),
            PETSC_DECIDE,
            s,
            val
        );
      }

      inline void local_zeros(const SizeType &s) override
      {
        zeros(
            comm().get(),
            type_override(),
            s,
            PETSC_DETERMINE
        );
      }

      inline void local_values(const SizeType &s, const Scalar &val) override
      {
        values(
            comm().get(),
            type_override(),
            s,
            PETSC_DETERMINE,
            val
        );
      }


      ///////////////////////////////////////////////////////////////////////////
      ////////////// OVERRIDES FOR Comparable ////////////////////////////
      ///////////////////////////////////////////////////////////////////////////


      bool equals(const PetscVector &other, const Scalar &tol = 0.0) const override;

      ///////////////////////////////////////////////////////////////////////////

       
        inline PetscVector()
        : vec_(nullptr), initialized_(false)
        {
            immutable_ = false;
        }

        inline ~PetscVector()
        {
            destroy();
        }

        PetscVector(const PetscVector &other)
        {
            if(other.vec_) {
                UTOPIA_REPORT_ALLOC("PetscVector::PetscVector(const PetscVector &)");
                PetscErrorHandler::Check(VecDuplicate(other.vec_, &vec_));
                PetscErrorHandler::Check(VecCopy(other.vec_, vec_));
                initialized_ = other.initialized_;
                ghost_values_ = other.ghost_values_;
            } else {
                vec_ = nullptr;
                initialized_ = false;
            }

            immutable_ = other.immutable_;
        }

        PetscVector(PetscVector &&other)
        : comm_(std::move(other.comm_)),
          vec_(std::move(other.vec_)),
          initialized_(std::move(other.initialized_)),
          ghost_values_(std::move(other.ghost_values_)),
          immutable_(std::move(other.immutable_))
        {
            other.vec_ = nullptr;
        }

        inline std::string name() const
        {
            const char *name;
            PetscObjectGetName((PetscObject)implementation(), &name);
            return name;
        }

        inline void rename(const std::string &name)
        {
            PetscObjectSetName((PetscObject)implementation(), name.c_str());
        }

        inline VecType type() const
        {
            VecType ret;
            VecGetType(implementation(), &ret);
            return ret;
        }

        inline VecType type_override() const
        {
            return VECSTANDARD;
        }

        bool has_type(VecType type) const;

        bool same_type(const PetscVector &other) const;

     

        inline Range range() const override
        {
            SizeType r_begin, r_end;
            VecGetOwnershipRange(implementation(), &r_begin, &r_end);
            return Range(r_begin, r_end);
        }

        inline MPI_Comm communicator() const {
            MPI_Comm comm = PetscObjectComm((PetscObject) implementation());
            assert(comm != MPI_COMM_NULL);
            return comm;
        }

        inline bool is_compatible(const PetscVector &other) const
        {
            return !is_null() && !other.is_null() && size() == other.size();
        }

        // assign operator
        inline PetscVector &operator=(const PetscVector &other) {
            if(this == &other) return *this;
            assert(!immutable_);



            if(is_compatible(other) && !other.has_ghosts()) {
                assert((same_type(other) || this->has_ghosts()) || this->comm().size()==1 && "Inconsistent vector types. Handle types properly before copying" );
                assert(local_size() == other.local_size() && "Inconsistent local sizes. Handle local sizes properly before copying.");
                PetscErrorHandler::Check(VecCopy(other.vec_, vec_));
                initialized_ = other.initialized_;
                immutable_ = other.immutable_;
                return *this;
            }

            destroy();

            if(other.vec_) {
                UTOPIA_REPORT_ALLOC("PetscVector::operator=(const PetscVector &)");
                PetscErrorHandler::Check(VecDuplicate(other.vec_, &vec_));
                PetscErrorHandler::Check(VecCopy(other.vec_, vec_));
                ghost_values_ = other.ghost_values_;

                initialized_ = other.initialized_;
                immutable_ = other.immutable_;
            } else {
                initialized_ = false;
            }

            return *this;
        }

        inline PetscVector &operator=(PetscVector &&other) {
            if(this == &other) return *this;
            assert(!immutable_);

            destroy();

            initialized_ = other.initialized_;
            immutable_ = other.immutable_;
            vec_ = other.vec_;
            other.vec_ = nullptr;
            ghost_values_ = std::move(other.ghost_values_);

            other.initialized_ = false;
            other.immutable_ = false;
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

        inline Vec &raw_type() {
            return vec_;
        }

        inline const Vec &raw_type() const {
            assert(vec_ != nullptr);
            return vec_;
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

        inline SizeType global_to_local(const SizeType g_index) const
        {
            SizeType begin, end;
            VecGetOwnershipRange(implementation(), &begin, &end);

            if(!ghost_values_.has_ghosts() || (g_index >= begin && g_index < end)) {
                assert(g_index < end && "Index out of local range.");
                return g_index - begin;
            }

            return ghost_values_.get_index(g_index);
        }

        inline void init_ghost_index(const std::vector<SizeType> &index)
        {
            ghost_values_.init_index(local_size(), index);
        }

 
        inline const Scalar &operator[](const SizeType index) const
        {
            // Scalar value;
            // VecGetValues(implementation(), 1, &index, &value);
            // return value;
            assert(range().inside(index));
                assert((readable_) && "use Read<Vector> r(vec) before using get. Check if you are using a copy of the vector");

            return readable_->data[index - readable_->range_begin];
        }

        inline void get(
            const std::vector<SizeType> &index,
            std::vector<Scalar> &values) const
        {
            std::size_t n = index.size();

            values.resize(n);

            if(!ghost_values_.has_ghosts()) {

                VecGetValues(
                    implementation(),
                    static_cast<SizeType>(index.size()),
                    &index[0],
                    &values[0]
                    );

            } else {

                const Scalar *array;
                Vec local_form;

                VecGhostGetLocalForm(implementation(), &local_form);
                VecGetArrayRead(local_form, &array);

                assert(local_form != nullptr);

                for(std::size_t i = 0; i < n; ++i) {
                    auto li = global_to_local(index[i]);
                    values[i] = array[li];
                }

                VecRestoreArrayRead(local_form, &array);
                VecGhostRestoreLocalForm(implementation(), &local_form);
            }
        }

        inline void set(
            const std::vector<SizeType> &indices,
            const std::vector<Scalar> &values)
        {
            assert(indices.size() == values.size());
            check_error( VecSetValues(implementation(), indices.size(), &indices[0], &values[0], INSERT_VALUES) );
        }

     
        inline void add_vector(
            const std::vector<SizeType> &indices,
            const std::vector<Scalar> &values)
        {
            assert(indices.size() == values.size());
            check_error( VecSetValues(implementation(), indices.size(), &indices[0], &values[0], ADD_VALUES) );
        }

 

        inline bool has_ghosts() const
        {
            return ghost_values_.has_ghosts();
        }

        inline void update_ghosts()
        {
            ghost_values_.update(vec_);
        }

        inline void make_immutable()
        {
            immutable_ = true;
        }

        //builders
        void repurpose(MPI_Comm comm, VecType type, SizeType n_local, SizeType n_global);
        inline void zeros(MPI_Comm comm, VecType type, SizeType n_local, SizeType n_global)
        {
            repurpose(comm, type, n_local, n_global);
            check_error( VecZeroEntries(implementation()) );

            assert(is_consistent());
        }

        inline void values(MPI_Comm comm, VecType type, SizeType n_local, SizeType n_global, Scalar value)
        {
            repurpose(comm, type, n_local, n_global);
            check_error( VecSet(implementation(), value) );

            assert(is_consistent());
        }

        void init(MPI_Comm comm, VecType type, SizeType n_local, SizeType n_global);
        void ghosted(MPI_Comm comm, SizeType local_size, SizeType global_size, const std::vector<SizeType> &index);
        void nest(MPI_Comm comm, SizeType nb, IS is[], Vec x[], const bool use_vec_nest_type = false);


        //ops
        ///this is y
        // inline void axpy(const Scalar &alpha, const PetscVector &x)
        // {
        //     assert(is_consistent());
        //     assert(x.is_consistent());

        //     check_error( VecAXPY(implementation(), alpha, x.implementation()) );
        // }

        ///this is y
        inline void axpby(const Scalar alpha, const PetscVector &x, const Scalar &beta)
        {
            assert(is_consistent());
            assert(x.is_consistent());

            check_error( VecAXPBY(implementation(), alpha, beta, x.implementation()) );
        }

        inline void zeros()
        {
            assert(initialized_);
            check_error( VecZeroEntries(implementation()) );
        }

   

        inline void e_mul(const PetscVector &other, PetscVector &result) const
        {
            assert(is_consistent());

            if(implementation() != result.vec_ && other.implementation() != result.vec_) {
                //if result is compatibe should not trigger a reallocation
                result.repurpose(communicator(), type(), local_size(), size());
            }

            check_error( VecPointwiseMult(result.implementation(), implementation(), other.implementation()) );

            assert(other.is_consistent());
        }

        // Scalar dot(const PetscVector &other) const;

        inline void e_div(const PetscVector &other, PetscVector &result) const
        {
            if(implementation() != result.vec_ && other.implementation() != result.vec_) {
                //if result is compatibe should not trigger a reallocation
                result.repurpose(communicator(), type(), local_size(), size());
            }

            check_error( VecPointwiseDivide(result.implementation(), implementation(), other.implementation()) );
        }

        inline void abs()
        {
            check_error( VecAbs(implementation()) );
        }

        inline void reciprocal()
        {
            check_error( VecReciprocal(implementation()) );
        }

        // inline void scale(const Scalar factor)
        // {
        //     check_error( VecScale(implementation(), factor) );
        // }

        inline void reciprocal(const Scalar numerator)
        {
            reciprocal();

            if(numerator == 1.) {
                return;
            }

            scale(numerator);
        }

     

        bool has_nan_or_inf() const;
        bool is_mpi() const;

        void resize(SizeType local_size, SizeType global_size);
        void resize(SizeType global_size)
        {
          resize(PETSC_DECIDE, global_size);
        }


        void select(
            const PetscIndexSet &index,
            PetscVector &result) const override;

        void select(const Range &global_range, PetscVector &result) const;


        void copy_from(Vec vec);


        inline bool read(const std::string &path)
        {
            return read(comm().get(), path);
        }

        bool read(MPI_Comm comm, const std::string &path);
 
        bool write(const std::string &path) const;
        bool write_binary(const std::string &path) const;
        bool write_matlab(const std::string &path) const;

        bool is_consistent() const;

        void convert_from(const Vec &mat);
        void convert_to(Vec &mat) const;
        void copy_data_to(Vec vec) const; 
        void copy_data_from(Vec vec); 

        void wrap(Vec &mat);

        inline void ghosted(
          const SizeType &local_size,
          const SizeType &global_size,
          const PetscArray<SizeType> &index
          )
        {
          ghosted(comm().get(), local_size, global_size, index);
        }

        inline std::string get_class() const override {
            return "PetscVector";
        }


        inline bool is_alias(const PetscVector &other) const
        {
            if(is_null()) return false;
            if(other.is_null()) return false;
            return raw_type() == other.raw_type();
        }

    private:
        PetscCommunicator comm_;

        Vec vec_;
        bool initialized_;

        GhostValues ghost_values_;

        //debug
        bool immutable_;

        class LocalView {
        public:
            Vec v;
            SizeType range_begin, range_end;
            Scalar *data;
            PetscErrorCode ierr;

            LocalView(Vec v_in)
            : v(v_in)
            {
                ierr = VecGetArray(v, &data); assert(ierr == 0);
                ierr = VecGetOwnershipRange(v, &range_begin, &range_end); assert(ierr == 0);
            }

            ~LocalView()
            {
                ierr = VecRestoreArray(v, &data); assert(ierr == 0);
            }
        };

        class ConstLocalView {
        public:
            const Vec v;
            const Scalar *data;
            SizeType range_begin, range_end;
            PetscErrorCode ierr;

            ConstLocalView(const LocalView &view)
            : v(nullptr),
            data(view.data),
            range_begin(view.range_begin),
            range_end(view.range_end)
            {}

            ConstLocalView(
                const Scalar * data_in,
                SizeType range_begin_in,
                SizeType range_end_in
                )
            : v(nullptr), data(data_in), range_begin(range_begin_in), range_end(range_end_in)
            {}

            ConstLocalView(const Vec v_in)
            : v(v_in)
            {
                ierr = VecGetArrayRead(v, &data); assert(ierr == 0);
                ierr = VecGetOwnershipRange(v, &range_begin, &range_end); assert(ierr == 0);
            }

            ~ConstLocalView()
            {
                if(v) {
                    ierr = VecGetArrayRead(v, &data); assert(ierr == 0);
                }
            }
        };

        std::unique_ptr<LocalView> writeable_;
        std::unique_ptr<ConstLocalView> readable_;

        inline static bool check_error(const SizeType err) {
            return PetscErrorHandler::Check(err);
        }

        bool is_cuda() const;
        bool is_root() const;

    };

}

#endif //UTOPIA_UTOPIA_PETSCVECTOR_H
