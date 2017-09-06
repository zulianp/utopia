#include "utopia_PerformanceTest.hpp"
#include "utopia.hpp"
#ifdef WITH_CUDA
#include  "test_problems/utopia_TestProblems.hpp"
#endif
#ifdef WITH_EIGEN_3
#include <Eigen/Dense>
#endif //WITH_EIGEN_3

namespace utopia {

	template<class M, class V>
	auto make_test_expr_1(const M &m, const V &v1,  const V &v2, const V &v3) -> decltype( abs(m * v1 + v2 - v3) )
	{
		return abs(m * v1 + v2 - v3);
	}

	template<class M, class V>
	auto make_test_expr_2(const M &m, const V &v1,  const V &v2, const V &v3) -> decltype( 10. * v1 + 5. * v2 )
	{
	     std::cout<<"I am in make_test_expr_2"<<std::endl;
             return 10. * v1 + 5. * v2;
	}

	template<class M, class V>
	auto make_test_expr_3(const M &m, const V &v1,  const V &v2, const V &v3) -> decltype( abs(m * pow2(v1) + sqrt(v2) - v3) )
	{
		return abs(m * pow2(v1) + sqrt(v2) - v3);
	}

	template<class M, class V>
	auto make_test_expr_4(const M &m, const V &v1,  const V &v2, const V &v3) -> decltype( abs(transpose(0.1 * (m * m) - m) * (m) * v3) )
	{
		return abs(transpose(0.1 * (m * m) - m) * (m) * v3);
	}

	template<class M, class V>
	auto make_test_expr_5(const M &m, const V &v1,  const V &v2, const V &v3) -> decltype( (0.1 * m - m) * v1 )
	{
		return (0.1 * m - m) * v1;
	}

	template<class M, class V>
	auto make_test_expr_6(const M &m, const V &v1,  const V &v2, const V &v3) -> decltype( (0.1 * pow2(v1) - pow2(v2)) + abs(v3) )
	{
		return (0.1 * pow2(v1) - pow2(v2)) + abs(v3);
	}

	static const int N_SIZES   = 10;
	static const int N[] 	   = { 10, 100, 1000, 5000, 10000, 50000, 100000, 500000, 1000000, 10000000 };
	static const int N_MIXED[] = { 10, 50,  100,  250, 	500,   750,   1000,   2000,   5000, 10000   };
	static const int N_RUNS    = 1;
	static const bool verbose  = true;


#ifdef WITH_CUDA

        template<class Matrix, class Vector>
        void test_poisson(const std::string &experiment_name)
        {
           Matrix A;
           Vector rhs;
           Matrix sol;
           Poisson_1D<Matrix,Vector> test_cuda;
           Chrono c;
       
                       
            for(int k = 0; k < N_SIZES; ++k) {
                    const int n = N[k];
                   // c.start();
                    test_cuda.getOperators(n,A,rhs);
                    std::cout<<"I am trying to multiply"<<std::endl;
                    c.start();
                    //disp((local_size(A)).get(0));
                    //disp((local_size(A)).get(1));
                    //disp(local_size(rhs));
                    sol=A*A;
                    c.stop();
                    if(verbose){
                            std::cout<<"problem size = "<< n << ",\t";
                            c.describe(std::cout);
                    }
             }
          
        }

#else
	template<class Matrix, class Vector>
	void test_program_inlined(const std::string &experiment_name)
	{

		for(int k = 0; k < N_SIZES; ++k) {
			const int n = N[k];

			Matrix m  = values(n, n, 0.001);
			Vector v1 = values(n, 0.1); 
			Vector v2 = values(n, 0.2);
			Vector v3 = values(n, 0.3);

			Chrono c;
			c.start();

			Vector res;
			double checksum = .0;      

			{
				Read<Matrix> r_m(m);
				Read<Vector> r_v1(v1);
				Read<Vector> r_v2(v2);
				Read<Vector> r_v3(v3);

				inline_eval(make_test_expr_2(m, v1, v2, v3), res); checksum += fabs(sum(res));
				inline_eval(make_test_expr_6(m, v1, v2, v3), res); checksum += fabs(sum(res));
			}

			c.stop();
			if(verbose) {
				std::cout<< "run inline " << ", check: " << checksum << ",\t";
				c.describe(std::cout);
			}
		}
	}

	template<class Matrix, class Vector>
	void test_program_mixed(const std::string &experiment_name)
	{
		for(int k = 0; k < N_SIZES; ++k) {
			const int n = N_MIXED[k];

			Matrix m  = values(n, n, 0.001);
			Vector v1 = values(n, 0.1);
			Vector v2 = values(n, 0.2); 
			Vector v3 = values(n, 0.3);

			double avg_time = 0;
			Chrono c;
			for(int i = 0; i < N_RUNS; ++i) {		 

				c.start();
				double checksum = .0;
				Vector res;
				res = make_test_expr_1(m, v1, v2, v3); checksum += fabs(sum(res));
				res = make_test_expr_3(m, v1, v2, v3); checksum += fabs(sum(res));
				res = make_test_expr_4(m, v1, v2, v3); checksum += fabs(sum(res));
				res = make_test_expr_5(m, v1, v2, v3); checksum += fabs(sum(res));
				c.stop();

				if(verbose) {
					std::cout<< "mixed run " << i << ", n: " << n << ", check: " << checksum << ",\t";
					avg_time += c.get_seconds();
					c.describe(std::cout);
				}
			}

			if(verbose) {
				avg_time /= N_RUNS;
				std::cout << "===================\n";
				std::cout << "exp: " << experiment_name << ",mixed," << n << "," << avg_time << "\n";
				std::cout << "===================\n";
			}
		}
	}

	template<class Matrix, class Vector>
	void test_program_vectors(const std::string &experiment_name)
	{
		for(int k = 0; k < N_SIZES; ++k) {
			const int n = N[k];

			Matrix m;
			Vector v1 = values(n, 0.1);
			Vector v2 = values(n, 0.2); 
			Vector v3 = values(n, 0.3);

			double avg_time = 0;
			Chrono c;
			for(int i = 0; i < N_RUNS; ++i) {		 

				c.start();
				double checksum = .0;
				Vector res;
				res = make_test_expr_2(m, v1, v2, v3); checksum += fabs(sum(res));
				res = make_test_expr_6(m, v1, v2, v3); checksum += fabs(sum(res));
				c.stop();

				if(verbose) {
					std::cout<< "vectors run " << i << ", n: " << n << ", check: " << checksum << ",\t";
					avg_time += c.get_seconds();
					c.describe(std::cout);
				}
			}


			if(verbose) {
				avg_time /= N_RUNS;
				std::cout << "===================\n";
				std::cout << "exp: " << experiment_name << ",vectors," << n << "," << avg_time << "\n";
				std::cout << "===================\n";
			}
		}
	}



 template<class Matrix, class Vector>
        void test_program(const std::string &experiment_name)
        {
           test_program_vectors<Matrix,Vector>(experimet_name);
           test_program_mixed<Matrix,Vector>(experiment_name);
        }





#endif

#ifdef WITH_EIGEN_3
	void test_program_eigen_3_vectors(const std::string &experiment_name)
	{
		typedef Eigen::MatrixXd Matrix;
		typedef Eigen::VectorXd Vector;

		for(int k = 0; k < N_SIZES; ++k) {
			const int n = N[k];

			Matrix m;
			Vector v1 = Vector::Constant(n, 0.1);
			Vector v2 = Vector::Constant(n, 0.2); 
			Vector v3 = Vector::Constant(n, 0.3);

			double avg_time = 0;
			Chrono c;
			for(int i = 0; i < N_RUNS; ++i) {		 
				c.start();
				double checksum = .0;
				Vector res;
				res = 10. * v1 + 5. * v2; 												  checksum += fabs(res.sum());
				res = (0.1 * v1.cwiseProduct(v1) - v2.cwiseProduct(v2) ) + v3.cwiseAbs(); checksum += fabs(res.sum());
				c.stop();

				if(verbose) {
					std::cout<< "vectors run " << i << ", n: " << n << ", check: " << checksum << ",\t";
					avg_time += c.get_seconds();
					c.describe(std::cout);
				}
			}

			if(verbose) {
				avg_time /= N_RUNS;
				std::cout << "===================\n";
				std::cout << "exp: " << experiment_name << ",vectors," << n << "," << avg_time << "\n";
				std::cout << "===================\n";
			}
		}
	}

	void test_program_eigen_3_mixed(const std::string &experiment_name)
	{
		typedef Eigen::MatrixXd Matrix;
		typedef Eigen::VectorXd Vector;

		for(int k = 0; k < N_SIZES; ++k) {
			const int n = N_MIXED[k];

			Matrix m  = Matrix::Constant(n, n, 0.001);
			Vector v1 = Vector::Constant(n, 0.1);
			Vector v2 = Vector::Constant(n, 0.2); 
			Vector v3 = Vector::Constant(n, 0.3);

			double avg_time = 0;
			Chrono c;
			for(int i = 0; i < N_RUNS; ++i) {		 
				c.start();
				double checksum = .0;
				Vector res;
				res = (m * v1 + v2 - v3).cwiseAbs(); 							  checksum += fabs(res.sum());
				res = (m * v1.cwiseProduct(v1) + v2.cwiseSqrt() - v3).cwiseAbs(); checksum += fabs(res.sum());
				res = ( (0.1 * (m * m) - m).transpose() * (m) * v3 ).cwiseAbs();  checksum += fabs(res.sum());
				res = (0.1 * m - m) * v1; 										  checksum += fabs(res.sum());
				c.stop();

				if(verbose) {
					std::cout<< "mixed run " << i << ", n: " << n << ", check: " << checksum << ",\t";
					avg_time += c.get_seconds();
					c.describe(std::cout);
				}
			}

			if(verbose) {
				avg_time /= N_RUNS;
				std::cout << "===================\n";
				std::cout << "exp: " << experiment_name << ",mixed," << n << "," << avg_time << "\n";
				std::cout << "===================\n";
			}
		}
	}

	void run_performance_test_eigen_3()
	{
		test_program_eigen_3_vectors("eigen3");
		test_program_eigen_3_mixed("eigen3");
	}
#endif //WITH_EIGEN_3	


	void run_performance_test()
	{
		if(verbose) {
			std::cout << "Begin: performance test" << std::endl; 
		}

#ifdef WITH_UTOPIA_OPENCL
		if(verbose) {
			std::cout << "------------------------------------\n";
			std::cout << "OpenCL: " << std::endl; 
			CLStats::instance().clear();
		}

		// CLContext::instance().set_current_device(2);
		// CLContext::instance().describe_current_setup();
		
		test_program<CLMatrixd, CLVectord>("opencl");

		if(verbose) {
			CLStats::instance().describe(std::cout);
		}
		
#endif //WITH_UTOPIA_OPENCL		

#ifdef WITH_BLAS
		if(verbose) {
			std::cout << "------------------------------------\n";
			std::cout << "Blas: " << std::endl; 
		}
		test_program_inlined<Matrixd, Vectord>("blas_inlined"); 
		test_program<Matrixd, Vectord>("blas");
#endif		

#ifdef WITH_PETSC
#ifdef WITH_CUDA
                if(verbose) {
                        std::cout << "------------------------------------\n";
                        std::cout << "PETSC_WITH_CUDA: " << std::endl;
                }
                test_poisson<DSMatrixd, DVectord>("petsc_with_CUDA");
#else

		if(verbose) {
			std::cout << "------------------------------------\n";
			std::cout << "PETSC: " << std::endl; 
		}
		test_program<DMatrixd, DVectord>("petsc");
#endif //WITH_CUDA		
#endif //WITH_PETSC



#ifdef WITH_EIGEN_3
		if(verbose) {
			std::cout << "------------------------------------\n";
			std::cout << "Eigen 3: " << std::endl; 
		}
		run_performance_test_eigen_3();		
#endif //WITH_EIGEN_3		
		if(verbose) {
			std::cout << "------------------------------------\n";
			std::cout << "End: performance test" << std::endl;
		}
	}

}
