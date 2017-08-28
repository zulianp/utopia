// #ifndef LIBMESH_CUTLIBPP_ADAPTERS_HPP
// #define LIBMESH_CUTLIBPP_ADAPTERS_HPP 


// #include "libmesh/mesh_inserter_iterator.h"

// #include "utopia_BoxAdapter.hpp"

// #include "moonolith_tree.hpp"
// #include "moonolith_n_tree_mutator_factory.hpp"
// #include "moonolith_n_tree_with_span_mutator_factory.hpp"
// #include "moonolith_n_tree_with_tags_mutator_factory.hpp"
// #include "moonolith_profiler.hpp"
// #include "moonolith_redistribute.hpp"
// #include "par_moonolith.hpp"

// #include <array>
// #include <vector>

// namespace utopia {
// 		template<typename T>
// 	inline void Print(const std::vector<T> &v, std::ostream &os)
// 	{
// 		for(auto i : v) {
// 			os << i << " ";
// 		}

// 		os << "\n";
// 	}

// 	static std::ostream &logger()
// 	{
// 		return moonolith::logger();
// 	}

// 		// class BoxAdapter : public moonolith::Serializable, public moonolith::Describable, public Box {
// 		// public:
// 		// 	void read(moonolith::InputStream &is) override 
// 		// 	{
// 		// 		auto &min = get_min();
// 		// 		auto &max = get_max();

// 		// 		int n;
// 		// 		is >> n;

// 		// 		reset(n);

// 		// 		for (int i = 0; i < n; ++i) {
// 		// 			min.el(i) << is;
// 		// 			max.el(i) << is;
// 		// 		}
// 		// 	}

// 		// 	void write(moonolith::OutputStream &os) const override
// 		// 	{
// 		// 		const int n = get_dims();
// 		// 		auto &min = get_min();
// 		// 		auto &max = get_max();

// 		// 		os << n;

// 		// 		for (int i = 0; i < n; ++i) {
// 		// 			os << min(i);
// 		// 			os << max(i);
// 		// 		}
// 		// 	}

// 		// 	void describe(std::ostream &os) const override
// 		// 	{
// 		// 		print(os);
// 		// 	}


// 		// 	inline bool isEmpty() const
// 		// 	{
// 		// 		return empty();
// 		// 	}

// 		// 	inline double getMinAt(const int coord) const
// 		// 	{
// 		// 		return get_min(coord);
// 		// 	}

// 		// 	inline double getMaxAt(const int coord) const
// 		// 	{
// 		// 		return get_max(coord);
// 		// 	}

// 		// 	inline void setMinAt(const int coord, const double value)
// 		// 	{
// 		// 		get_min().el(coord) = value;
// 		// 	}

// 		// 	inline void setMaxAt(const int coord, const double value)
// 		// 	{
// 		// 		get_max().el(coord) = value;
// 		// 	}

// 		// 	inline void clear()
// 		// 	{
// 		// 		reset();
// 		// 	}

// 		// 	inline int nDims() const {

// 		// 		return get_dims();
// 		// 	}
// 		// };


// 		// template<int Dimension>
// 		// class BoxBoxAdapter : public moonolith::Describable, public moonolith::Serializable {
// 		// public:
// 		// 	typedef utopia::BoxAdapter StaticBound;

// 		// 	void read(moonolith::InputStream &is)
// 		// 	{
// 		// 		is >> static_;
// 		// 		bool is_empty;
// 		// 		is >> is_empty;
// 		// 		if(!is_empty) { is >> dynamic_; };
// 		// 	}


// 		// 	void write(moonolith::OutputStream &os) const
// 		// 	{
// 		// 		os << static_;
// 		// 		bool is_empty = dynamic_.isEmpty();
// 		// 		os << is_empty;
// 		// 		if(!is_empty) { os << dynamic_; }
// 		// 	}

// 		// 	bool intersects(const BoxBoxAdapter &bound) const
// 		// 	{
// 		// 		return static_.intersects(bound.static_) && dynamic_.intersects(bound.dynamic_);
// 		// 	}

// 		// 	bool intersects(const BoxBoxAdapter &bound, const double tol) const
// 		// 	{
// 		// 		return static_.intersects(bound.static_, tol) && dynamic_.intersects(bound.dynamic_, tol);
// 		// 	}

// 		// 	bool intersects(const BoxAdapter &bound) const
// 		// 	{
// 		// 		return static_.intersects(bound);
// 		// 	}

// 		// 	inline double getMinAt(const int coord) const
// 		// 	{
// 		// 		return static_.getMinAt(coord);
// 		// 	}

// 		// 	inline double getMaxAt(const int coord) const
// 		// 	{
// 		// 		return static_.getMaxAt(coord);
// 		// 	}

// 		// 	inline void setMinAt(const int coord, const double value)
// 		// 	{
// 		// 		static_.setMinAt(coord, value);
// 		// 	}

// 		// 	inline void setMaxAt(const int coord, const double value)
// 		// 	{
// 		// 		static_.setMaxAt(coord, value);
// 		// 	}

// 	 //       //expands to contain the union of this and CompositeBound
// 		// 	BoxBoxAdapter &operator +=(const BoxBoxAdapter &bound)
// 		// 	{
// 		// 		static_ += bound.static_;
// 		// 		if(dynamic_.isEmpty()) {
// 		// 			dynamic_ = bound.dynamic_;
// 		// 		} else if(!bound.dynamic_.isEmpty()) {
// 		// 			dynamic_ += bound.dynamic_;
// 		// 		}
// 		// 		return *this;
// 		// 	}

// 		// 	bool isEmpty() const
// 		// 	{
// 		// 		return static_.isEmpty();
// 		// 	}

// 		// 	void clear()
// 		// 	{
// 		// 		static_.reset(Dimension);
// 		// 		dynamic_.reset(Dimension);
// 		// 	}

// 		// 	BoxBoxAdapter()
// 		// 	{
// 		// 		clear();
// 		// 	}

// 		// 	void describe(std::ostream &os) const
// 		// 	{
// 		// 		os << "Static bound:\n"  << static_  << "\n";
// 		// 		os << "Dynamic bound:\n";
// 		// 		dynamic_.describe(os);
// 		// 		os << "\n";
// 		// 	}

// 		// 	inline BoxAdapter &static_bound() { return static_; }
// 		// 	inline const BoxAdapter &static_bound() const { return static_; }

// 		// 	inline BoxAdapter &dynamic_bound() { return dynamic_; }
// 		// 	inline const BoxAdapter &dynamic_bound() const { return dynamic_; }

// 		// private:
// 		// 	BoxAdapter static_;
// 		// 	BoxAdapter dynamic_;
// 		// };


// 		template<int Dimension>
// 	class ElementAdapter : public moonolith::Serializable {
// 	public:
// 		inline int tag() const
// 		{
// 			return tag_;
// 		}

// 		const BoxBoxAdapter<Dimension> &bound() const
// 		{
// 			return bound_;
// 		}

// 		BoxBoxAdapter<Dimension> &bound()
// 		{
// 			return bound_;
// 		}

// 		void apply_read_write(moonolith::Stream &stream) override
// 		{
// 			stream & bound_;
// 			stream & element_;
// 			stream & element_handle_;
// 		}

// 		ElementAdapter(LibMeshFESpaceBase &fe, const libMesh::dof_id_type &element, const long element_handle, const int tag)
// 		: fe_(&fe), element_(element), element_handle_(element_handle), tag_(tag), dof_map_(nullptr)
// 		{   
// 			assert(element < fe.mesh().n_elem());

// 			libMesh::Elem * e = fe.mesh().elem(element);

// 			for (libMesh::dof_id_type i = 0; i < e->n_nodes(); ++i) {
// 				const libMesh::Point &p = fe.mesh().node(e->node(i));

// 				    for(int d = 0; d < Dimension; ++d) {
//                     p_a[d] = p(d);
//                 }

// 				bound_.static_bound()  += p;
// 				bound_.dynamic_bound() += p;


// 			}

// 		}



// 		ElementAdapter()
// 		:fe_(nullptr) , element_(-1), element_handle_(-1), tag_(-1), dof_map_(nullptr) {}


// 		inline long handle() const
// 		{
// 			return element_handle_;
// 		}


// 		inline long element() const
// 		{
// 			return element_;
// 		}


// 		libMesh::Elem * get()
// 		{
// 			assert(fe_);

// 				//std::cout<<"I AM IN GET"<<std::endl;
			
// 			assert(element_ < fe_->mesh().n_local_elem());

// 			return fe_->mesh().elem(element_);
// 		}

// 		const libMesh::Elem * get() const
// 		{
// 			assert(fe_);
// 			return fe_->mesh().elem(element_);
// 		}


// 		inline const  LibMeshFESpaceBase  &space() const
// 		{
// 			assert(fe_);
// 			return *fe_;
// 		}


// 		void set_dof_map(std::vector<long> * ptr)
// 		{

// 			dof_map_ = ptr;
// 		}

// 		inline const std::vector<long> &dof_map() const
// 		{
// 			assert(dof_map_);
// 			return *dof_map_;
// 		}


// 	private:
// 		LibMeshFESpaceBase * fe_;
// 		libMesh::dof_id_type element_;
// 		long element_handle_;
// 		int tag_;
// 		BoxBoxAdapter<Dimension> bound_;
// 		std::vector<long> * dof_map_;
// 	};


//     template<int Dimension>
// 	class SurfaceElementAdapter : public moonolith::Serializable {
// 	public:
// 		inline int tag() const
// 		{
// 			return tag_;
// 		}

// 		const BoxBoxAdapter<Dimension> &bound() const
// 		{
// 			return bound_;
// 		}

// 		BoxBoxAdapter<Dimension> &bound()
// 		{
// 			return bound_;
// 		}

// 		void apply_read_write(moonolith::Stream &stream) override
// 		{
// 			stream & bound_;
// 			stream & element_;
// 			stream & element_handle_;
// 		}

//         SurfaceElementAdapter(LibMeshFESpaceBase &fe, const libMesh::dof_id_type &element, const long element_handle, const int tag, /*std::vector<long> &map,*/ const libMesh::Real blow_up)
// 		: fe_(&fe), element_(element), element_handle_(element_handle), tag_(tag), dof_map_(nullptr)
// 		{
// 			assert(element < fe.mesh().n_elem());

// 			libMesh::Elem * e = fe.mesh().elem(element);

// 			const int dim = fe.mesh().mesh_dimension();

// 			libMesh::Point o, u, v, n, c, p;

// 			for(uint side = 0; side < e->n_sides(); ++side) {

// 				libMesh::UniquePtr<libMesh::Elem> s = e->side(side);


// 				if(!s->on_boundary()) continue;


// 				auto side_ptr = e->build_side_ptr(side);

// 				compute_side_normal(dim, *side_ptr, n);

// 				if(fix_normal_orientation(*e, side, n)) {
// 					std::cerr << "[Warning] face with wrong orientation detected\n" << std::endl;
// 				}

// 				assert( n.contract(side_ptr->centroid()-e->centroid()) > 0 );

// 				for(uint i = 0; i < side_ptr->n_nodes(); ++i) {
// 					c = side_ptr->point(i);

// 					n *= blow_up;
// 					p = c;
// 					p += n;

// 					std::array<double, Dimension> p_a;

// 					for(int d = 0; d < Dimension; ++d) {
// 						p_a[d] = p(d);
// 					}

// 					bound_.static_bound()  += p_a;
// 					bound_.dynamic_bound() += p_a;
// 					p = c;
// 					n *= 0.01;
// 					p -=n;

// 					for(int d = 0; d < Dimension; ++d) {
// 						p_a[d] = p(d);
// 					}

// 					bound_.static_bound()  += p_a;
// 					bound_.dynamic_bound() += p_a;
// 				}
// 			}

// 		}



// 		SurfaceElementAdapter()
// 		:fe_(nullptr) , element_(-1), element_handle_(-1), tag_(-1), dof_map_(nullptr) {}


// 		inline long handle() const
// 		{
// 			return element_handle_;
// 		}


// 		inline long element() const
// 		{
// 			return element_;
// 		}


// 		libMesh::Elem * get()
// 		{
// 			assert(fe_);
// 			assert(element_ < fe_->mesh().n_local_elem());

// 			return fe_->mesh().elem(element_);
// 		}

// 		const libMesh::Elem * get() const
// 		{
// 			assert(fe_);
// 			return fe_->mesh().elem(element_);
// 		}


// 		inline const  LibMeshFESpaceBase  &space() const
// 		{
// 			assert(fe_);
// 			return *fe_;
// 		}


// 		void set_dof_map(std::vector<long> * ptr)
// 		{

// 			dof_map_ = ptr;
// 		}

// 		inline const std::vector<long> &dof_map() const
// 		{
// 			assert(dof_map_);
// 			return *dof_map_;
// 		}

// 	private:
// 		LibMeshFESpaceBase * fe_;
// 		libMesh::dof_id_type element_;
// 		long element_handle_;
// 		int tag_;
// 		BoxBoxAdapter<Dimension> bound_;
// 		std::vector<long> * dof_map_;
// 	};

// 	template<int _Dimension>
// 	class TreeTraits {
// 	public:
// 		enum {
// 			Dimension = _Dimension
// 		};

// 		typedef utopia::BoxBoxAdapter<Dimension> Bound;
// 		typedef utopia::ElementAdapter<Dimension> DataType;
// 	};

//     template<int _Dimension>
// 	class SurfaceTreeTraits {
// 	public:
// 		enum {
// 			Dimension = _Dimension
// 		};

// 		typedef utopia::BoxBoxAdapter<Dimension> Bound;
// 		typedef utopia::SurfaceElementAdapter<Dimension> DataType;
// 	};

// 	template<int Dimension>
// 	class LibMeshTree : public moonolith::Tree< TreeTraits<Dimension> > {
// 	public:
// 		typedef TreeTraits<Dimension> Traits;

// 		LibMeshTree() {};

// 		static std::shared_ptr<LibMeshTree> New(const std::shared_ptr<moonolith::Predicate> &predicate,
// 			const int maxElementsXNode=moonolith::DEFAULT_REFINE_MAX_ELEMENTS,
// 			const int maxDepth=moonolith::DEFAULT_REFINE_DEPTH
// 			) {
// 			using namespace moonolith;

// 			auto tree = std::make_shared<LibMeshTree>();
// 			auto factory = std::make_shared<NTreeWithTagsMutatorFactory<LibMeshTree>>(predicate);
// 			factory->set_refine_params(maxElementsXNode, maxDepth);
// 			tree->set_mutator_factory(factory);
// 			return tree;
// 		}
// 	};

//     template<int Dimension>
// 	class SurfaceLibMeshTree : public moonolith::Tree< SurfaceTreeTraits<Dimension> > {
// 	public:
// 		typedef SurfaceTreeTraits<Dimension> Traits;

// 		SurfaceLibMeshTree() {};

// 		static std::shared_ptr<SurfaceLibMeshTree> New(const std::shared_ptr<moonolith::Predicate> &predicate,
// 			const int maxElementsXNode=moonolith::DEFAULT_REFINE_MAX_ELEMENTS,
// 			const int maxDepth=moonolith::DEFAULT_REFINE_DEPTH
// 			) {
// 			using namespace moonolith;
// 			auto tree = std::make_shared<SurfaceLibMeshTree>();
// 			auto factory = std::make_shared<NTreeWithTagsMutatorFactory<SurfaceLibMeshTree>>(predicate);
// 			factory->set_refine_params(maxElementsXNode, maxDepth);
// 			tree->set_mutator_factory(factory);
// 			return tree;
// 		}
// 	};

// 	class ElementDofMap : public moonolith::Serializable {
// 	public:
// 		void read(moonolith::InputStream &is) override
// 		{
// 			int n;
// 			is >> n;
// 			global.resize(n);
// 			is.read(&global[0], n);
// 		}

// 		void write(moonolith::OutputStream &os) const override
// 		{
// 			int n = global.size();
// 			os << n;
// 			os.write(&global[0], n);
// 		}


// 		inline bool empty() const
// 		{
// 			return global.empty();
// 		}

// 		std::vector<long> global;


// 	};

//    // add new SurfaceElementAdapter

// }

// #endif //LIBMESH_CUTLIBPP_ADAPTERS_HPP
