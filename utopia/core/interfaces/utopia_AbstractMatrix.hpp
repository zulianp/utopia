#ifndef UTOPIA_ABSTRACT_MATRIX_HPP
#define UTOPIA_ABSTRACT_MATRIX_HPP

#include <iostream>

namespace utopia {

    template<typename Scalar_, typename SizeType_>
    class Traits< AbstractMatrix<Scalar_, SizeType_> > {
    public:
        using Scalar = Scalar_;
        using SizeType = SizeType_;
    };

    //parallel types, collective operations
    template<typename Scalar_, typename SizeType_>
    class AbstractMatrix {
    public:
        //print function
        inline void describe() const
        {
            std::cout << " (matrix)" << std::endl;
            // impl_->describe();
        }

    };
}

#endif //UTOPIA_ABSTRACT_MATRIX_HPP


