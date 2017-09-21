#ifndef UTOPIA_VOLUME_TRANSFER_BENCHMARK_HPP
#define UTOPIA_VOLUME_TRANSFER_BENCHMARK_HPP 

namespace libMesh {
	class LibMeshInit;
}

namespace utopia {
	void run_volume_transfer_benchmark(libMesh::LibMeshInit &init);
}

#endif //UTOPIA_VOLUME_TRANSFER_BENCHMARK_HPP
