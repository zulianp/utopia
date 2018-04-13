#include "utopia_petsc_Vector.hpp"

#include <set>
#include <cstring>

namespace utopia {

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
			other.describe();
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
			VecCreate(comm, &vec_);
		} else {
			if(comm != PetscObjectComm((PetscObject)vec_) || has_type(type)) {
				destroy();

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
		
		check_error( VecCreate(comm, &vec_) );
		check_error( VecSetFromOptions(vec_) );
		check_error( VecSetType(vec_, type) );

		check_error( VecSetSizes(vec_, n_local, n_global) );
		
		assert(vec_ != nullptr);
		initialized_ = true;

		assert(is_consistent());
	}
	
	void PetscVector::ghosted(MPI_Comm comm,
		PetscInt local_size,
		PetscInt global_size,
		const std::vector<PetscInt> &index)
	{

		assert(!immutable_);

		destroy();
		
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

	bool PetscVector::is_nan_or_inf() const
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
		const std::vector<PetscInt> &index,
		PetscVector &result) const
	{
		MPI_Comm comm = communicator();
		IS is_in;
		VecScatter scatter_context;

		result.repurpose(comm, this->type(), index.size(), PETSC_DETERMINE);

		ISCreateGeneral(comm, index.size(), &index[0], PETSC_USE_POINTER, &is_in);
		VecScatterCreate(implementation(), is_in, result.implementation(), nullptr, &scatter_context);

		VecScatterBegin(scatter_context, implementation(), result.implementation(), INSERT_VALUES, SCATTER_FORWARD);
		VecScatterEnd(scatter_context,   implementation(), result.implementation(), INSERT_VALUES, SCATTER_FORWARD);
		
		ISDestroy(&is_in);
		VecScatterDestroy(&scatter_context);
	}

	void PetscVector::copy_from(Vec vec)
	{
		destroy();

		VecDuplicate(vec, &implementation());
		VecCopy(vec, implementation());

		assert(vec_ != nullptr);
		set_initialized(true);
		assert(is_consistent());
	}

	bool PetscVector::read(MPI_Comm comm, const std::string &path)
	{
		destroy();

		PetscViewer fd;
		
		bool err = check_error( PetscViewerBinaryOpen(comm, path.c_str(), FILE_MODE_READ, &fd) );
		
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

		result.write_lock();
		
		for(PetscInt r_this = rr.begin(); r_this < rr.end(); ++r_this) {
			const PetscInt r_selection = r_this - global_range.begin();
			result.set(r_selection, get(r_this));
		}
		
		result.write_unlock();
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
}
