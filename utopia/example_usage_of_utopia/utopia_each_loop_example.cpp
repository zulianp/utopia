// #define UTOPIA_DISABLE_UNROLLING
#include "utopia.hpp"
#include <cmath>
#include <algorithm>

using namespace utopia;

template<class Vector>
static void run(const int n)
{
	std::cout << "Utopia" << std::endl;
	Vector v = local_values(n, 1.);
	std::cout << backend(v).info().get_name() << std::endl;

	Chrono c;
	c.start();

	each_write(v, [](const SizeType i) -> double {
		return 10.;
	});

	c.stop();
	std::cout << c << std::endl;
}

static void run_raw(const int n)
{
	std::cout << "Raw" << std::endl;
	std::vector<double> v(n, 1.);

	Chrono c;
	c.start();
	for(auto &vi : v) {
		vi = 10.;
	}
	c.stop();
	std::cout << c << std::endl;
}


static void run_stl_each(const int n)
{
	std::cout << "ForEach" << std::endl;
	std::vector<double> v(n, 1.);

	Chrono c;
	c.start();
	
	std::for_each(std::begin(v), std::end(v), [](double &value) {
		value = 10.;
	});

	c.stop();
	std::cout << c << std::endl;
}


static void run_for(const int n)
{
	std::cout << "For" << std::endl;
	std::vector<double> v(n, 1.);

	Chrono c;
	c.start();

	For<64>::apply(
		0,
		n,
		[&v](const std::size_t i) {
			v[i] = 10.;
		}
		);

	c.stop();
	std::cout << c << std::endl;
}

static void run_all(const int n)
{
	run_raw(n);
	run_stl_each(n);
	run_for(n);
	//if it has compiled with blas or petsc WITH_BLAS or WITH_PETSC macros are available (if you want to make it compile no matter the utopia installation)
#ifdef WITH_PETSC    
//run with petsc types 
	run<DVectord>(n);
#endif //WITH_PETSC 

#ifdef WITH_BLAS    
//run with blas types 
	run<Vectord>(n);
#endif //WITH_BLAS    

#ifdef WITH_TRILINOS    
//run with trilinos types 
	run<TVectord>(n);
#endif //WITH_TRILINOS    
}

int main(int argc, char** argv)
{
	using namespace utopia;

	Utopia::Init(argc, argv);

	if(argc < 2) {
		static const int N = 90000000;
		run_all(N);
		return Utopia::Finalize();
	}

	auto n = atol(argv[1]);
	run_all(n);
	return Utopia::Finalize();
}