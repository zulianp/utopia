#ifndef UTOPIA_SCRIPT_HPP
#define UTOPIA_SCRIPT_HPP

#include <iostream>

// #ifdef EXPORT_PETSC

// #else

// #endif

namespace algebra {
    class SparseMatrix {
    public:
        inline SparseMatrix()
        {
            std::cout << "HI" << std::endl;
        }

        inline ~SparseMatrix()
        {
            std::cout << "BYE" << std::endl;
        }
    };

}

#endif //UTOPIA_SCRIPT_HPP
