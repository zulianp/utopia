#ifndef UTOPIA_COMMUNICATOR_HPP
#define UTOPIA_COMMUNICATOR_HPP

#include "utopia_Base.hpp"

#ifdef WITH_MPI
#include <mpi.h>
#endif

namespace utopia {

	class Communicator {
	public:
		virtual ~Communicator() {}
		virtual int rank() const = 0;
		virtual int size() const = 0;
		virtual Communicator * clone() const = 0;

// #ifdef WITH_MPI
// 		virtual MPI_Comm get() const = 0;
// #endif //WITH_MPI

	};

	/**
	 * @brief usefull for non-mpi backends as a default object to return/provide
	 */
	class SelfCommunicator final : public Communicator {
	public:
		int rank() const noexcept override { return 0; }
		int size() const noexcept override { return 1; }

		SelfCommunicator * clone() const noexcept override 
		{
			return new SelfCommunicator();
		}

		inline bool conjunction(const bool &val) const {
			return val;
		}

		inline bool disjunction(const bool &val) const {
			return val;
		}

		template<typename T>
		inline T sum(const T &val) const {
			return val;
		}

		template<typename T>
		inline T min(const T &val) const {
			return val;
		}

		template<typename T>
		inline T max(const T &val) const {
			return val;
		}

#ifdef WITH_MPI
		inline MPI_Comm get() const //override 
		{
			return MPI_COMM_SELF;
		}
#endif //WITH_MPI
	};

	template<typename T>
	class MPIType {};


#ifdef WITH_MPI
	template<>
	class MPIType<double> {
	public:
		inline static MPI_Datatype value() { return MPI_DOUBLE; }
	};

	template<>
	class MPIType<float> {
	public:
		inline static MPI_Datatype value() { return MPI_FLOAT; }
	};

	template<>
	class MPIType<long> {
	public:
		inline static MPI_Datatype value() { return MPI_LONG; }
	};

	template<>
	class MPIType<int> {
	public:
		inline static MPI_Datatype value() { return MPI_INT; }
	};

	template<>
	class MPIType<char> {
	public:
		inline static MPI_Datatype value() { return MPI_CHAR; }
	};

	class MPICommunicator : public Communicator {
	public:
		~MPICommunicator() {}

		virtual MPI_Comm get() const = 0;

		/////////////////////////////////////////////

		bool conjunction(const bool &val) const;
		bool disjunction(const bool &val) const;

		inline int rank() const override
		{
			int ret;
			MPI_Comm_rank(get(), &ret);
			return ret;
		}

		inline int size() const override
		{
			int ret;
			MPI_Comm_size(get(), &ret);
			return ret;
		}

		template<typename T>
		inline T sum(const T &val) const {
			T ret = val;
			MPI_Allreduce( MPI_IN_PLACE, &ret, 1, MPIType<T>::value(), MPI_SUM, get() );
			return ret;
		}

		template<typename T>
		inline T min(const T &val) const {
			T ret = val;
			MPI_Allreduce( MPI_IN_PLACE, &ret, 1, MPIType<T>::value(), MPI_MIN, get() );
			return ret;
		}

		template<typename T>
		inline T max(const T &val) const {
			T ret = val;
			MPI_Allreduce( MPI_IN_PLACE, &ret, 1, MPIType<T>::value(), MPI_MAX, get() );
			return ret;
		}
	};

#endif //WITH_MPI

}

#endif //UTOPIA_COMMUNICATOR_HPP
