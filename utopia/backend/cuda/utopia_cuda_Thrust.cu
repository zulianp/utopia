
#include <thrust/host_vector.h>
#include <thrust/device_vector.h>
#include <thrust/transform.h>
#include <thrust/sequence.h>
#include <thrust/copy.h>
#include <thrust/fill.h>
#include <thrust/replace.h>
#include <thrust/functional.h>

void run_stupid()
{
    thrust::device_vector<int> d(40);
    thrust::fill(d.begin(), d.end(), 1);

      //Replace all the ones in Y with tens
    thrust::replace(d.begin(), d.end(), 1, 4);
    // thrust::copy(d.begin(), d.end(), std::ostream_iterator<int>(std::cout, "\n"));
}

void run_stupid_2(thrust::device_vector<double> &d)
{
    d.resize(40);
    thrust::fill(d.begin(), d.end(), 1);

    //Replace all the ones in Y with tens
    thrust::replace(d.begin(), d.end(), 1, 5);
    // thrust::copy(d.begin(), d.end(), std::ostream_iterator<int>(std::cout, "\n"));
}
