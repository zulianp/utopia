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

    class Communicator;

}  // namespace utopia

namespace scripting {

    using Scalar = double;
    using SizeType = int;
    using LocalSizeType = int;
    using Factory = utopia::AlgebraFactory<Scalar, SizeType>;
    //using Comm = utopia::Communicator;
    //using Layout = utopia::Layout<Communicator, LocalSizeType, SizeType>;
  

    void init();
    void finalize();
    void print_info();

    class Communicator {
    public:
        using CommunicatorImpl = utopia::Communicator;

        Communicator();
        ~Communicator();

        CommunicatorImpl * get_communicator() const {
            return impl_;
        }
      
    private:
        CommunicatorImpl* impl_;
    };

    class Layout {
        public:
            using LayoutImpl = utopia::Layout<utopia::Communicator, 1, LocalSizeType, SizeType>;
        
        Layout(const Communicator &comm, int Order, LocalSizeType local_size, SizeType global_size);
        ~Layout();

        private:
            LayoutImpl * impl_;
            const Communicator &comm_;
            int Order_;
            LocalSizeType local_size_;
            SizeType global_size_; 
        //     Comm comm_;
        // LocalSizeType local_size_[Order];
        // SizeType size_[Order];
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
        Scalar value_;
    };

    
    class SparseMatrix {
    public:
        using MatrixImpl = utopia::AbstractMatrix<Scalar, SizeType>;

        SparseMatrix();
        ~SparseMatrix();
        void print_info();

    private:
        MatrixImpl* impl_;
    };

    
}  // namespace scripting




#endif  // UTOPIA_SCRIPT_HPP
