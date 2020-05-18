// #ifndef UTOPIA_TEXT_OUPUT_HPP
// #define UTOPIA_TEXT_OUPUT_HPP

// #include <fstream>
// #include "utopia_MPI.hpp"
// #include "utopia_Path.hpp"

// namespace utopia {

//     template<class Matrix>
//     bool write_text(const Path &path, const Matrix &mat)
//     {
//         using SizeType = UTOPIA_SIZE_TYPE(Matrix);
//         using Scalar   = UTOPIA_SCALAR(Matrix);

//        int size = utopia::comm_size(mat);
//        int rank = utopia::comm_rank(mat);

//        int nnz = 0;
//        for(SizeType r = 0; r < size; ++r) {
//            if(r == 0) {
//                nnz = 0;
//                each_read(mat, [&nnz](const SizeType, const SizeType, const Scalar) {
//                    ++nnz;
//                });

//                MPI_Allreduce( MPI_IN_PLACE, &nnz, 1, MPI_INT, MPI_SUM, comm );
//            }

//            if(r == rank) {
//                std::ofstream os;

//                if(r == 0) {
//                    os.open(path);
//                    Size s = size(mat);
//                    os << s.get(0) << " " << nnz << "\n";
//                } else {
//                    os.open(path, std::ofstream::out | std::ofstream::app);
//                }

//                if(!os.good()) {
//                    std::cerr << "invalid path: " << path << std::endl;
//                    continue;
//                }

//                each_read(mat, [&os](const SizeType i, const SizeType j, const Scalar value) {
//                    os << i << " " << j << " " << value << "\n";
//                });

//                os.flush();
//                os.close();
//            }

//            MPI_Barrier(comm);
//        }

//        return true;
//    }
// }

// #endif //UTOPIA_TEXT_OUPUT_HPP
