#include "utopia_script.hpp"
#include "utopia_Instance.hpp"
#include "utopia_Version.hpp"

#include <iostream>

namespace algebra {

    void init(int argc, char *argv[])
    {
        utopia::Utopia::Init(argc, argv);
    }

    void init()
    {
        //FIXME?
        int argc = 1;
        std::string argv = "utopia_script";
        char * argv_ptr = &argv[0];
        init(argc, &argv_ptr);
    }

    void print_info()
    {
        std::cout << "Utopia\nversion: " << UTOPIA_VERSION << std::endl;
    }

    void finalize()
    {
        utopia::Utopia::Finalize();
    }

    SparseMatrix::SparseMatrix()
    {
        std::cout << "HI" << std::endl;
    }

    SparseMatrix::~SparseMatrix()
    {
        std::cout << "BYE" << std::endl;
    }
}
