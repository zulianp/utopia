#include "utopia.hpp"
#include <string>

using namespace utopia;

//faster with sizes larger than 10
void mat_mat_mult(const int rows_A, const int cols_A, const int cols_B, const double * A, const double * B, double *work, double *C)
{
	for(int j = 0; j < cols_B; ++j) {
		const double * B_j = &B[j];

		for(int k = 0; k < cols_A; ++k) {
			work[k] = B_j[k * cols_B];
		}

		for(int i = 0; i < rows_A; ++i) {
			const int offset_i = i * cols_A;
			const int offset   = i * cols_B + j;
			
			const double * A_i = &A[offset_i];
			double temp = 0;
			
			for(int k = 0; k < cols_A; ++k) {
				temp += A_i[k] * work[k];
			}

			C[offset] = temp;
		}
	}
}

void mat_mat_mult_naive(const int rows_A, const int cols_A, const int cols_B, const double * A, const double * B, double *C)
{
	for(int j = 0; j < cols_B; ++j) {
		const double * B_j = &B[j]; 

		for(int i = 0; i < rows_A; ++i) {
			
			const int offset_i = i * cols_A;
			const int offset   = i * cols_B + j;

			const double * A_i = &A[offset_i];
			double temp = 0;

			for(int k = 0; k < cols_A; ++k) {
				temp += A_i[k] * B_j[k * cols_B];
			}

			C[offset] = temp;
		}
	}
}

void print_vector(const std::vector<double> &vec)
{
	for(auto v : vec) {
		std::cout << v << " ";
	}

	std::cout << std::endl;
}

//performance of mat*mat
int main(const int argc, char * argv[])
{
	Utopia::Init(argc, argv);

	if(argc < 2) {
		std::cout << "usage: rows_A [cols_A] [cols_B]" << std::endl;
		return 1;
	}
	const int rows_A = std::stoi(argv[1]);
	const int cols_A = std::stoi(argc > 2 ? argv[2] : argv[1]);
	const int cols_B = std::stoi(argc > 3 ? argv[3] : argv[1]);

	std::cout << rows_A << " x " << cols_A << " x " << cols_B << std::endl;

	std::vector<double> A(rows_A * cols_A, 1.0);
	std::vector<double> B(cols_A * cols_B, 1.0);
	std::vector<double> C(rows_A * cols_B, 0.0);
	std::vector<double> work(cols_A, 0.0);

	Chrono c;
	c.start();
	mat_mat_mult(rows_A, cols_A, cols_B, &A[0], &B[0], &work[0], &C[0]);
	c.stop();

	std::cout << "cached cols: " << std::endl;
	c.describe(std::cout);

	// print_vector(C);

	c.start();
	mat_mat_mult_naive(rows_A, cols_A, cols_B, &A[0], &B[0], &C[0]);
	c.stop();

	std::cout << "naive method: " << std::endl;
	c.describe(std::cout);

	// print_vector(C);

	Matrixd m_A = values(rows_A, cols_A, 1.0);
	Matrixd m_B = values(cols_A, cols_A, 1.0);
	Matrixd m_C;

	c.start();
	m_C = m_A * m_B;
	c.stop();

	std::cout << "blas: " << std::endl;
	c.describe(std::cout);

	// disp(m_C);
	return Utopia::Finalize();
}
