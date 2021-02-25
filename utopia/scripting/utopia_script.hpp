#ifndef UTOPIA_SCRIPT_HPP
#define UTOPIA_SCRIPT_HPP

//#include "utopia_AbstractVector.hpp"
namespace utopia {
    // Forward declarations


    template <typename Scalar, typename SizeType>
    class AbstractVector;

    template <typename Scalar, typename SizeType>
    class AbstractMatrix;

    template <typename Scalar, typename SizeType>
    class AlgebraFactory;

    template <class Comm, int Order, typename LocalSizeType_, typename SizeType_>
    class Layout;

    // class Clonable
    class Communicator;
    // class SelfCommunicator;


    //template <class Comm, typename LocalSizeType_, typename SizeType_>
    //class Layout;

}  // namespace utopia

namespace scripting {

    using Scalar = double;
    using SizeType = int;
    using LocalSizeType = int;
    using Factory = utopia::AlgebraFactory<Scalar, SizeType>;
    //using Communicator = utopia::Communicator;
    //using Layout = utopia::Layout<Communicator, LocalSizeType, SizeType>;
  

    void init();
    void finalize();
    void print_info();

    class Communicator {
    public:
        using CommunicatorImpl = utopia::Communicator;

        Communicator();
        ~Communicator();
      
    private:
        CommunicatorImpl* impl_;
    };

    class Layout {
        public:
            using LayoutImpl = utopia::Layout<utopia::Communicator, 1, LocalSizeType, SizeType>;
        
        Layout(const Communicator *comm, int Order, LocalSizeType local_size, SizeType global_size);
        ~Layout();

        private:
            const Communicator *comm;
            int Order;
            LocalSizeType local_size;
            SizeType global_size; 
            LayoutImpl * impl_;
    };

    class Vector {
    public:
        using VectorImpl = utopia::AbstractVector<Scalar, SizeType>;

        Vector(); // Layout layout, Scalar value
        ~Vector();
        void print_info();
        void set(const Scalar &val);

    private:
        VectorImpl* impl_;
       // Layout l;
        Scalar value;
    };

    // class Clonable {
    // public:
    //     virtual ~Clonable() = default;

    //     /** @brief This method copies the relevant
    //      * settings but does not have to copy all the state variables
    //      * such as buffers. Maybe we should change its name to somthing
    //      * more suitable ...
    //      */
    //     virtual Clonable* clone() const = 0;
    // };

    // class Communicator : public Clonable {
    // public:
    //     ~Communicator() override = default;
    //     virtual int rank() const = 0;
    //     virtual int size() const = 0;
    //     Communicator *clone() const override = 0;
    //     virtual void barrier() const = 0;

    //     virtual bool same(const Communicator &other) const { return size() == other.size(); }
    //     virtual bool conjunction(const bool &val) const = 0;
    //     virtual bool disjunction(const bool &val) const = 0;
    //     inline bool is_root() const { return rank() == 0; }
    // };

    // class SelfCommunicator : public Communicator {
    // public:
    //     int rank() const noexcept override { return 0; }
    //     int size() const noexcept override { return 1; }

    //     inline static SelfCommunicator world() { return SelfCommunicator(); }

    //     inline static SelfCommunicator self() { return SelfCommunicator(); }

    //     void barrier() const override {}

    //     SelfCommunicator *clone() const noexcept override { return new SelfCommunicator(); }

    //     inline static SelfCommunicator &get_default() {
    //         static SelfCommunicator instance_;
    //         return instance_;
    //     }

    //     inline bool conjunction(const bool &val) const override { return val; }

    //     inline bool disjunction(const bool &val) const override { return val; }
    // };



    // class Layout {
    // public:
   
    //     const Comm &comm() const { return comm_; }

    //     template <typename... Args>
    //     Layout(Comm comm, Args &&... args) : comm_(std::move(comm)) {
    //         init(std::forward<Args>(args)...);
    //     }

    //     Layout() {
    //         for (int i = 0; i < Order; ++i) {
    //             local_size_[i] = 0;
    //             size_[i] = 0;
    //         }
    //     }

    //     template <class OtherComm, typename OtherSizeType, typename OtherLocalSizeType>
    //     Layout(const Layout<OtherComm, Order, OtherLocalSizeType, OtherSizeType> &other) : comm_(other.comm()) {
    //         std::copy(&other.local_size(0), &other.local_size(0) + Order, local_size_);
    //         std::copy(&other.size(0), &other.size(0) + Order, size_);
    //     }

    //     inline void init(const LocalSizeType local_size[Order], SizeType size[Order]) {
    //         std::copy(local_size, local_size + Order, local_size_);
    //         std::copy(size, size + Order, size_);
    //     }

    //     inline void init(const Size &local, const Size &global) {
    //         auto &local_data = local.data();
    //         auto &global_data = global.data();

    //         assert(int(local_data.size()) <= Order);
    //         assert(int(global_data.size()) <= Order);

    //         std::copy(local_data.begin(), local_data.end(), local_size_);
    //         std::copy(global_data.begin(), global_data.end(), size_);
    //     }

    //     inline void init(const LocalSizeType local_size, SizeType size) {
    //         local_size_[0] = local_size;
    //         size_[0] = size;
    //     }

    //     friend void disp(const Layout &layout, std::ostream &os = std::cout) {
    //         os << "comm: " << layout.comm_.rank() << "/" << layout.comm_.size() << "\n";

    //         for (int i = 0; i < Order; ++i) {
    //             os << layout.local_size_[i] << " " << layout.size_[i] << "\n";
    //         }
    //     }

    // private:
    //     SelfCommunicator comm_;
    //     LocalSizeType local_size_[Order];
    //     SizeType size_[Order];
    // };

    
    class SparseMatrix {
    public:
        using MatrixImpl = utopia::AbstractMatrix<Scalar, SizeType>;

        SparseMatrix();
        ~SparseMatrix();
        void print_info();

    private:
        MatrixImpl* impl_;
    };

    

    // inline PetscTraits::Layout layout(Vec v) {
    //     PetscInt n_local, n_global;
    //     VecGetLocalSize(v, &n_local);
    //     VecGetSize(v, &n_global);

    //     MPI_Comm comm = PetscObjectComm((PetscObject)v);
    //     assert(comm != MPI_COMM_NULL);
    //     return layout(PetscTraits::Communicator(comm), n_local, n_global);
    // }


}  // namespace scripting




#endif  // UTOPIA_SCRIPT_HPP
