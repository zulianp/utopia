#ifndef UTOPIA_UTOPIA_PETSCVECTOR_H
#define UTOPIA_UTOPIA_PETSCVECTOR_H

#include "utopia_petsc_Error.hpp"

#include "utopia_Range.hpp"
#include "utopia_Base.hpp"
#include "utopia_Range.hpp"

#include "petscvec.h"

#include <map>
#include <vector>
#include <limits>

namespace utopia {
	
	class PetscVector {
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
		
		inline PetscVector()
		: vec_(nullptr), initialized_(false)
		{
#ifndef NDEBUG
			immutable_ = false;
#endif            
		}
		
		inline ~PetscVector()
		{
			destroy();
		}
		
		PetscVector(const PetscVector &other)
		{
			if(other.vec_) {
				PetscErrorHandler::Check(VecDuplicate(other.vec_, &vec_));
				PetscErrorHandler::Check(VecCopy(other.vec_, vec_));
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

		inline std::string name() const
		{
			const char *name;
			PetscObjectGetName((PetscObject)implementation(), &name);
			return name;
		}

		inline void set_name(const std::string &name)
		{
			PetscObjectSetName((PetscObject)implementation(), name.c_str());
		}

		inline VecType type() const
		{
			VecType ret;
			VecGetType(implementation(), &ret);
			return ret;
		}

		bool has_type(VecType type) const;

	 	bool same_type(const PetscVector &other) const;

		inline PetscInt local_size() const
		{
			PetscInt ret;
			VecGetLocalSize(implementation(), &ret);
			return ret;
		}
		
		inline PetscInt size() const
		{
			PetscInt ret;
			VecGetSize(implementation(), &ret);
			return ret;
		}

		inline Range range() const
		{
			PetscInt r_begin, r_end;
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
				assert((same_type(other) || this->has_ghosts()) && "Inconsistent vector types. Handle types properly before copying" );
				assert(local_size() == other.local_size() && "Inconsistent local sizes. Handle local sizes properly before copying.");
				PetscErrorHandler::Check(VecCopy(other.vec_, vec_));
				initialized_ = other.initialized_;
							
#ifndef NDEBUG
				immutable_ = other.immutable_;
#endif 
				return *this;
			}
			
			destroy();
			
			if(other.vec_) {
				PetscErrorHandler::Check(VecDuplicate(other.vec_, &vec_));
				PetscErrorHandler::Check(VecCopy(other.vec_, vec_));
				ghost_values_ = other.ghost_values_;

				initialized_ = other.initialized_;
							
#ifndef NDEBUG
				immutable_ = other.immutable_;
#endif 
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
			VecView(implementation(), PETSC_VIEWER_STDOUT_(communicator()));
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
			VecGetOwnershipRange(implementation(), &begin, &end);
			
			if(!ghost_values_.has_ghosts() || (g_index >= begin && g_index < end)) {
				assert(g_index < end && "Index out of local range.");
				return g_index - begin;
			}
			
			return ghost_values_.get_index(g_index);
		}
		
		inline void init_ghost_index(const std::vector<PetscInt> &index)
		{
			ghost_values_.init_index(local_size(), index);
		}

		inline PetscScalar get(const PetscInt index) const
		{
			PetscScalar value;
			VecGetValues(implementation(), 1, &index, &value);
			return value;
		}
		
		inline void get(const std::vector<PetscInt> &index,
					    std::vector<PetscScalar> &values) const
		{
			std::size_t n = index.size();
			
			values.resize(n);
			
			if(!ghost_values_.has_ghosts()) {
				
				VecGetValues(
							 implementation(),
							 static_cast<PetscInt>(index.size()),
							 &index[0],
							 &values[0]
							 );
				
			} else {
				
				const PetscScalar *array;
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
			const std::vector<PetscInt> &indices,
			const std::vector<PetscScalar> &values) 
		{
			assert(indices.size() == values.size());
			check_error( VecSetValues(implementation(), indices.size(), &indices[0], &values[0], INSERT_VALUES) );
		}

		inline void set(const PetscScalar value)
		{
			check_error( VecSet(implementation(), value) );
		}

		inline void add_vector(
			const std::vector<PetscInt> &indices,
			const std::vector<PetscScalar> &values) 
		{
			assert(indices.size() == values.size());
			check_error( VecSetValues(implementation(), indices.size(), &indices[0], &values[0], ADD_VALUES) );
		}

		inline void set(const PetscInt index, PetscScalar value)
		{
			check_error( VecSetValues(implementation(), 1, &index, &value, INSERT_VALUES) );
		}

		inline void add(const PetscInt index, PetscScalar value)
		{
			check_error( VecSetValues(implementation(), 1, &index, &value, ADD_VALUES) );
		}
		
		inline bool has_ghosts() const
		{
			return ghost_values_.has_ghosts();
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
			check_error( VecZeroEntries(implementation()) );
		}

		inline void values(MPI_Comm comm, VecType type, PetscInt n_local, PetscInt n_global, PetscScalar value)
		{
			repurpose(comm, type, n_local, n_global);
			check_error( VecSet(implementation(), value) );
		}
		
		void init(MPI_Comm comm, VecType type, PetscInt n_local, PetscInt n_global);
		void ghosted(MPI_Comm comm, PetscInt local_size, PetscInt global_size, const std::vector<PetscInt> &index);


		
		//ops
		///this is y
		inline void axpy(const PetscScalar &alpha, const PetscVector &x)
		{
			check_error( VecAXPY(implementation(), alpha, x.implementation()) );
		}

		///this is y
		inline void axpby(const PetscScalar alpha, const PetscVector &x, const PetscScalar &beta)
		{
			check_error( VecAXPBY(implementation(), alpha, beta, x.implementation()) );
		}
		
		inline void zeros()
		{
			assert(initialized_);
			check_error( VecZeroEntries(implementation()) );
		}

		//reductions
		inline PetscReal norm2() const {
			PetscReal val;
			check_error( VecNorm(implementation(), NORM_2, &val) );
			return val;
		}
		
		inline PetscReal norm1() const {
			PetscReal val;
			check_error( VecNorm(implementation(), NORM_1, &val) );
			return val;
		}
		
		inline PetscReal norm_infty() const {
			PetscReal val;
			check_error( VecNorm(implementation(), NORM_INFINITY, &val) );
			return val;
		}

		inline PetscScalar sum() const {
			PetscScalar result = 0;
			check_error( VecSum(implementation(), &result) );
			return result;
		}

		inline PetscScalar min() const {
			PetscScalar result = std::numeric_limits<PetscScalar>::max();
			check_error( VecMin(implementation(), nullptr, &result) );
			return result;
		}

		inline PetscScalar max() const {
			PetscScalar result = -std::numeric_limits<PetscScalar>::max();
			check_error( VecMax(implementation(), nullptr, &result) );
			return result;
		}

		inline void e_mul(const PetscVector &other, PetscVector &result) const
		{
			if(implementation() != result.vec_ && other.implementation() != result.vec_) {
				//if result is compatibe should not trigger a reallocation
				result.repurpose(communicator(), type(), local_size(), size());
			}

			check_error( VecPointwiseMult(result.implementation(), implementation(), other.implementation()) );
		}

		inline PetscScalar dot(const PetscVector &other) const
		{
			PetscScalar result;
			check_error( VecDot(implementation(), other.implementation(), &result) );
			return result;
		}

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

		inline void scale(const PetscScalar factor)
		{
			check_error( VecScale(implementation(), factor) );
		}

		inline void reciprocal(const PetscScalar numerator)
		{
			reciprocal();

			if(numerator == 1.) {	
				return;
			}

			scale(numerator);
		}

		inline void read_lock() {}
		inline void read_unlock() {}

		inline void write_lock() {}
		inline void write_unlock()
		{
			VecAssemblyBegin(implementation());
			VecAssemblyEnd(implementation());
			
			set_initialized(true);
			update_ghosts();
		}

		bool is_nan_or_inf() const;
		bool is_mpi() const;

		void resize(PetscInt local_size, PetscInt global_size);


		void select(
			const std::vector<PetscInt> &index,
			PetscVector &result) const;

		void select(const Range &global_range, PetscVector &result) const;


		void copy_from(Vec vec);

		bool read(MPI_Comm comm, const std::string &path);

		bool write(const std::string &path) const;
		bool write_matlab(const std::string &path) const;

		
	private:
		Vec vec_;
		bool initialized_;
		
		GhostValues ghost_values_;
		
		//debug
#ifndef NDEBUG
		bool immutable_;
#endif //NDEBUG     
		
		inline static bool check_error(const PetscInt err) {
			return PetscErrorHandler::Check(err);
		}
	};
	
}

#endif //UTOPIA_UTOPIA_PETSCVECTOR_H
