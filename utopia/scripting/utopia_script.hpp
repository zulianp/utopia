#ifndef UTOPIA_SCRIPT_HPP
#define UTOPIA_SCRIPT_HPP

namespace algebra {
    // void init(int argc, char *argv[]);

    void init();
    void finalize();
    void print_info();

    class SparseMatrix {
    public:
        SparseMatrix();
        ~SparseMatrix();
        void print_info();
    };
}

#endif //UTOPIA_SCRIPT_HPP
