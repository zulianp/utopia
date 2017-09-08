#include "utopia_PerformanceTestCuda.hpp"
#include "utopia_CUDABackendTPL.hpp"
#include "utopia_CUDAMatrix.hpp"
#include "utopia_CUDAVector.hpp"

namespace utopia{

//   int KK=7;
  
   static int SIZES[7] = {16, 32, 64, 128, 256, 512, 1024};
 
   void test_CUDA(const std::string &experiment_name)
  {
          std::cout << "Running " << experiment_name << std::endl;
          int times = 10;
	  CUDAVector<double> v;
	  CUDAMatrix<double> m;
	  CUDAVector<double> resv;
	  CUDAMatrix<double> resm; 

	  for (int sz=0; sz < 7; sz ++) {

            double * v_h = (double *) malloc (SIZES[sz]*sizeof(double));

	    m.values.reserve(SIZES[sz]*SIZES[sz]);
	    m.values.resize(SIZES[sz]*SIZES[sz]);

	    v.values.reserve(SIZES[sz]); 
	    v.values.resize(SIZES[sz]);
 
	    resv.values.reserve(SIZES[sz]); 
	    resv.values.resize(SIZES[sz]);
 
	    resm.values.reserve(SIZES[sz]*SIZES[sz]);
	    resm.values.resize(SIZES[sz]*SIZES[sz]);

            //fill matrices and vectors
            cuda_double::build_values(SIZES[sz], 2, v);
            cuda_double::build_values(SIZES[sz], SIZES[sz], 2, m);
            
            double * v_ptr = thrust::raw_pointer_cast(&v.values[0]);
            

            //thrust::copy(v_ptr, v_ptr+SIZES[sz], v_h);            
            cudaMemcpy(v_h, v_ptr, SIZES[sz]*sizeof(double), cudaMemcpyDeviceToHost);
            
  //          for (int j=0; j<SIZES[sz]; ++j) std::cout << "v_h[0] = " << v_h[j] << std::endl; 

	    utopia::Chrono c;
	    //test matrix-vector multiplication
	    c.start();
	    for (int i = 0; i < times; ++i) {
	      cuda_double::mat_vec_mul( m, v, resv);
	    }
	    c.stop();
            std::cout << "WARNING we are performing the same opeartion for "<< times <<" times\n";

	    std::cout << "CUDA Matrix-Vector Results\n";
            std::cout << "Size Problem\n";
            std::cout << "Size Matrix Array==>"<<SIZES[sz]*SIZES[sz]<<"\n";
            std::cout << "Size Problem\n";
            std::cout << "Size Vector Array==>"<<SIZES[sz]<<"\n";
	    c.describe(std::cout);

	    // test matrix-matrix multiplication
	    c.start();
	    for (int i = 0; i < times; ++i) {
	      cuda_double::mat_mat_mul( m, m, SIZES[sz], resm);
	    }
	    c.stop();
	    std::cout << "CUDA Matrix-Matrix Results\n";
            std::cout << "Size Problem\n";
            std::cout << "Size Matrix Array==>"<<SIZES[sz]*SIZES[sz]<<"\n";
	    c.describe(std::cout);
	  

	  }

  }

   void test_CPU(const std::string &experiment_name)
  {
   int times = 10;
    std::cout << "Running " << experiment_name << std::endl;  
      for (int sz=0; sz < 7; sz ++) {
        int length_v=SIZES[sz];
        int length_m=SIZES[sz]*SIZES[sz]; 
        double v[length_v];
        double m[length_m];
        double res_v[length_v];
        double res_m[length_m];
        for (int k=0; k<length_v; ++k) v[k]=2;
        for (int k=0; k<length_m; ++k) m[k]=2;
        utopia::Chrono c;
        c.start();
        for (int ii = 0; ii < times; ++ii) {
        //     std::cout <<"ciao"<<std::endl;
             for (int l=0; l<length_v; ++l){
                 double res=0;
                  for (int k=0; k<length_v; ++k){
                        double m_x=m[l*length_v+k];
                        double v_x=v[k];
                        res+= m_x*v_x;
                      }
                      res_v[l]=res;
                 }
        }
       c.stop();
       std::cout << "WARING we are performing the same opeartion for "<< times <<" times\n";
       
       std::cout << "CPU Matrix-Vector Results\n";
       std::cout << "Size Problem\n";
       std::cout << "Size Matrix Array==>"<<SIZES[sz]*SIZES[sz]<<"\n";
       std::cout << "Size Problem\n";
       std::cout << "Size Vector Array==>"<<SIZES[sz]<<"\n";
       c.describe(std::cout);
      
//       for (int j=0; j<SIZES[sz]; ++j) std::cout << "v_res[0] = " << res_v[j] << std::endl;       

       c.start();
      for (int ii = 0; ii < times; ++ii){
           for (int l=0; l<length_v; ++l){
                 for (int ll=0; ll<length_v; ++ll){   
                      double res=0;
                      for(int k=0; k<length_v; ++k){
                           double m_x=m[l*length_v+k];
                           double m_y=m[k*length_v+ll];
                           res+= m_x*m_y;
                       }  
                   res_m[l*length_v+ll]=res;
               }
          }
      }
     c.stop();
     std::cout << "CPU Matrix-Matrix Results\n";
     std::cout << "Size Problem\n";
     std::cout << "Size Matrix Array==>"<<SIZES[sz]*SIZES[sz]<<"\n";
     c.describe(std::cout);
              
   }

}
  //int size = 
  //cudaMalloc(&m , WIDTH*WIDTH*sizeof (int) ) ;
  //cudaMalloc(&m , WIDTH*WIDTH*sizeof (int) ) 
  void run_performance_CUDA_test(){          
 
	   test_CUDA("CUDA");

           test_CPU("CPU");
     } 

}
