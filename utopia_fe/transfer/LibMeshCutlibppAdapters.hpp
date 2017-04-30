#ifndef LIBMESH_CUTLIBPP_ADAPTERS_HPP
#define LIBMESH_CUTLIBPP_ADAPTERS_HPP 

#include "cutlibpp.hpp"
#include "cutlibpp_Base.hpp"
#include "cutlibpp_Tree.hpp"
#include "cutlibpp_NTreeMutatorFactory.hpp"
#include "cutlibpp_NTreeWithSpanMutatorFactory.hpp"
#include "cutlibpp_NTreeWithTagsMutatorFactory.hpp"
#include "cutlibpp_API.hpp"
#include "libmesh/mesh_inserter_iterator.h"
#include "express_Profiler.hpp"
#include "express_Redistribute.hpp"
#include "MapSparseMatrix.hpp"

namespace utopia {
		template<typename T>
		inline void Print(const std::vector<T> &v, std::ostream &os)
		{
			for(auto i : v) {
				os << i << " ";
			}

			os << "\n";
		}

		static std::ostream &logger()
		{
			return express::Express::Instance().logger().os();
		}

		class BoxAdapter : public cutk::Serializable, public cutk::Describable, public Box {
		public:
			void read(cutk::InputStream &is) override 
			{
				auto &min = get_min();
				auto &max = get_max();

				int n;
				is >> n;

				reset(n);

				for (int i = 0; i < n; ++i) {
					min.el(i) << is;
					max.el(i) << is;
				}
			}

			void write(cutk::OutputStream &os) const override
			{
				const int n = get_dims();
				auto &min = get_min();
				auto &max = get_max();

				os << n;

				for (int i = 0; i < n; ++i) {
					os << min(i);
					os << max(i);
				}
			}

			void describe(std::ostream &os) const override
			{
				print(os);
			}


			inline bool isEmpty() const
			{
				return empty();
			}

			inline double getMinAt(const int coord) const
			{
				return get_min(coord);
			}

			inline double getMaxAt(const int coord) const
			{
				return get_max(coord);
			}

			inline void setMinAt(const int coord, const double value)
			{
				get_min().el(coord) = value;
			}

			inline void setMaxAt(const int coord, const double value)
			{
				get_max().el(coord) = value;
			}

			inline void clear()
			{
				reset();
			}

			inline int nDims() const {

				return get_dims();
			}
		};


		template<int Dimension>
		class BoxBoxAdapter : public cutk::Describable, public cutk::Serializable {
		public:
			typedef utopia::BoxAdapter StaticBound;

			void read(cutk::InputStream &is)
			{
				is >> static_;
				bool is_empty;
				is >> is_empty;
				if(!is_empty) { is >> dynamic_; };
			}


			void write(cutk::OutputStream &os) const
			{
				os << static_;
				bool is_empty = dynamic_.isEmpty();
				os << is_empty;
				if(!is_empty) { os << dynamic_; }
			}

			bool intersects(const BoxBoxAdapter &bound) const
			{
				return static_.intersects(bound.static_) && dynamic_.intersects(bound.dynamic_);
			}

			bool intersects(const BoxBoxAdapter &bound, const double tol) const
			{
				return static_.intersects(bound.static_, tol) && dynamic_.intersects(bound.dynamic_, tol);
			}

			bool intersects(const BoxAdapter &bound) const
			{
				return static_.intersects(bound);
			}

			inline double getMinAt(const int coord) const
			{
				return static_.getMinAt(coord);
			}

			inline double getMaxAt(const int coord) const
			{
				return static_.getMaxAt(coord);
			}

			inline void setMinAt(const int coord, const double value)
			{
				static_.setMinAt(coord, value);
			}

			inline void setMaxAt(const int coord, const double value)
			{
				static_.setMaxAt(coord, value);
			}

	       //expands to contain the union of this and CompositeBound
			BoxBoxAdapter &operator +=(const BoxBoxAdapter &bound)
			{
				static_ += bound.static_;
				if(dynamic_.isEmpty()) {
					dynamic_ = bound.dynamic_;
				} else if(!bound.dynamic_.isEmpty()) {
					dynamic_ += bound.dynamic_;
				}
				return *this;
			}

			bool isEmpty() const
			{
				return static_.isEmpty();
			}

			void clear()
			{
				static_.reset(Dimension);
				dynamic_.reset(Dimension);
			}

			BoxBoxAdapter()
			{
				clear();
			}

			void describe(std::ostream &os) const
			{
				os << "Static bound:\n"  << static_  << "\n";
				os << "Dynamic bound:\n";
				dynamic_.describe(os);
				os << "\n";
			}

			inline BoxAdapter &staticBound() { return static_; }
			inline const BoxAdapter &staticBound() const { return static_; }

			inline BoxAdapter &dynamicBound() { return dynamic_; }
			inline const BoxAdapter &dynamicBound() const { return dynamic_; }

		private:
			BoxAdapter static_;
			BoxAdapter dynamic_;
		};


		template<int Dimension>
		class ElementAdapter : public cutk::Serializable {
		public:
			inline int tag() const
			{
				return tag_;
			}

			const BoxBoxAdapter<Dimension> &getBound() const
			{
				return bound_;
			}

			BoxBoxAdapter<Dimension> &getBound()
			{
				return bound_;
			}

			void applyRW(cutk::Stream &stream) 
			{
				stream & bound_;
				stream & element_;
				stream & element_handle_;
			}

			ElementAdapter(LibMeshFESpaceBase &fe, const libMesh::dof_id_type &element, const long element_handle, const int tag)
			: fe_(&fe), element_(element), element_handle_(element_handle), tag_(tag), dof_map_(nullptr)
			{   
                assert(element < fe.mesh().n_elem());
               
				libMesh::Elem * e = fe.mesh().elem(element);
                
				for (libMesh::dof_id_type i = 0; i < e->n_nodes(); ++i) {
					const libMesh::Point &p = fe.mesh().node(e->node(i));
					bound_.staticBound()  += p;
					bound_.dynamicBound() += p;


				}

			}



			ElementAdapter()
			:fe_(nullptr) , element_(-1), element_handle_(-1), tag_(-1), dof_map_(nullptr) {}


			inline long handle() const
			{
				return element_handle_;
			}


			inline long element() const
			{
				return element_;
			}


			libMesh::Elem * get()
			{
				assert(fe_);

				//std::cout<<"I AM IN GET"<<std::endl;
			
                assert(element_ < fe_->mesh().n_local_elem());

				return fe_->mesh().elem(element_);
			}

			const libMesh::Elem * get() const
			{
				assert(fe_);
				return fe_->mesh().elem(element_);
			}


			inline const  LibMeshFESpaceBase  &space() const
			{
				assert(fe_);
				return *fe_;
			}


			void set_dof_map(std::vector<long> * ptr)
			{

				dof_map_ = ptr;
			}

			inline const std::vector<long> &dof_map() const
			{
				assert(dof_map_);
				return *dof_map_;
			}


		private:
			LibMeshFESpaceBase * fe_;
			libMesh::dof_id_type element_;
			long element_handle_;
			int tag_;
			BoxBoxAdapter<Dimension> bound_;
			std::vector<long> * dof_map_;
		};
    
    
    template<int Dimension>
    class SurfaceElementAdapter : public cutk::Serializable {
    public:
        inline int tag() const
        {
            return tag_;
        }
        
        const BoxBoxAdapter<Dimension> &getBound() const
        {
            return bound_;
        }
        
        BoxBoxAdapter<Dimension> &getBound()
        {
            return bound_;
        }
        
        void applyRW(cutk::Stream &stream)
        {
            stream & bound_;
            stream & element_;
            stream & element_handle_;
        }
        
        SurfaceElementAdapter(LibMeshFESpaceBase &fe, const libMesh::dof_id_type &element, const long element_handle, const int tag, /*std::vector<long> &map,*/ const libMesh::Real blow_up)
        : fe_(&fe), element_(element), element_handle_(element_handle), tag_(tag), dof_map_(nullptr)
        {
            assert(element < fe.mesh().n_elem());
            
            libMesh::Elem * e = fe.mesh().elem(element);
            
            const int dim = fe.mesh().mesh_dimension();
            
            libMesh::Point o, u, v, n, c, p;
            
            for(uint side = 0; side < e->n_sides(); ++side) {
                
                libMesh::UniquePtr<libMesh::Elem> s = e->side(side);
                
                
                if(!s->on_boundary()) continue;

                
                auto side_ptr = e->build_side_ptr(side);
                
                compute_side_normal(dim, *side_ptr, n);
                
                if(fix_normal_orientation(*e, side, n)) {
                    std::cerr << "[Warning] face with wrong orientation detected\n" << std::endl;
                }
                
                assert( n.contract(side_ptr->centroid()-e->centroid()) > 0 );
                
                for(uint i = 0; i < side_ptr->n_nodes(); ++i) {
                    c = side_ptr->point(i);
//                    bound_.staticBound()  += c;
//                    bound_.dynamicBound() += c;
//                    //instead of normal blow up
//                    bound_.staticBound().enlarge(blow_up);
//                    bound_.dynamicBound().enlarge(blow_up);
//                     n.print();
//                     std::cout << "\n";
                    
                    n *= blow_up;
                    p = c;
                    p += n;
                    
//                     c.print();
//                     std::cout << "n = "<< n<<"\n";
                    
                    bound_.staticBound()  += p;
                    bound_.dynamicBound()  += p;
                    p = c;
                    n *= 0.01;
                    p -=n;
                    bound_.staticBound()  += p;
                    bound_.dynamicBound()  += p;
                }
            }
            
            
            // element_boxes[index].print(std::cout);
            
            
            
            //
            //            for (libMesh::dof_id_type i = 0; i < e->n_nodes(); ++i) {
            //                const libMesh::Point &p = fe.mesh().node(e->node(i));
            //                bound_.staticBound()  += p;
            //                bound_.dynamicBound() += p;
            //               //da modificare
            //
            //            }
            
        }
        
        
        
        SurfaceElementAdapter()
        :fe_(nullptr) , element_(-1), element_handle_(-1), tag_(-1), dof_map_(nullptr) {}
        
        
        inline long handle() const
        {
            return element_handle_;
        }
        
        
        inline long element() const
        {
            return element_;
        }
        
        
        libMesh::Elem * get()
        {
            assert(fe_);
            
            //std::cout<<"I AM IN GET"<<std::endl;
            
            assert(element_ < fe_->mesh().n_local_elem());
            
            return fe_->mesh().elem(element_);
        }
        
        const libMesh::Elem * get() const
        {
            assert(fe_);
            return fe_->mesh().elem(element_);
        }
        
        
        inline const  LibMeshFESpaceBase  &space() const
        {
            assert(fe_);
            return *fe_;
        }
        
        
        void set_dof_map(std::vector<long> * ptr)
        {
            
            dof_map_ = ptr;
        }
        
        inline const std::vector<long> &dof_map() const
        {
            assert(dof_map_);
            return *dof_map_;
        }
        
        
    private:
        LibMeshFESpaceBase * fe_;
        libMesh::dof_id_type element_;
        long element_handle_;
        int tag_;
        BoxBoxAdapter<Dimension> bound_;
        std::vector<long> * dof_map_;
    };



		template<int _Dimension>
		class TreeTraits {
		public:
			enum {
				Dimension = _Dimension
			};

			typedef utopia::BoxBoxAdapter<Dimension> Bound;
			typedef utopia::ElementAdapter<Dimension> DataType;
            //typedef utopia::SurfaceElementAdapter<Dimension> DataType;

		};
    
    
    
        template<int _Dimension>
        class SurfaceTreeTraits {
        public:
            enum {
               Dimension = _Dimension
          };
        
         typedef utopia::BoxBoxAdapter<Dimension> Bound;
         //typedef utopia::ElementAdapter<Dimension> DataType;
         typedef utopia::SurfaceElementAdapter<Dimension> DataType;
        
    };
    
    
    

		template<int Dimension>
		class LibMeshTree : public cutlibpp::Tree< TreeTraits<Dimension> > {
		public:
			typedef TreeTraits<Dimension> Traits;
			
			LibMeshTree() {};


			static cutk::shared_ptr<LibMeshTree> New(const cutk::shared_ptr<cutlibpp::Predicate> &predicate,
				const int maxElementsXNode=cutlibpp::DEFAULT_REFINE_MAX_ELEMENTS,
				const int maxDepth=cutlibpp::DEFAULT_REFINE_DEPTH
				) {
				using namespace cutlibpp;
				cutk::shared_ptr<LibMeshTree> tree = cutk::make_shared<LibMeshTree>();
				cutk::shared_ptr<NTreeWithTagsMutatorFactory < LibMeshTree> > factory =
				cutk::make_shared<NTreeWithTagsMutatorFactory < LibMeshTree> >(predicate);
				factory->setRefineParams(maxElementsXNode, maxDepth);
				std::cout<<"LibMeshTree::maxElements = "<<maxElementsXNode<<std::endl;
				std::cout<<"LibMeshTree::maxDepth = "<<maxDepth<<std::endl;
				tree->setMutatorFactory(factory);
				return tree;
			}


		};
    
    
       template<int Dimension>
       class SurfaceLibMeshTree : public cutlibpp::Tree< SurfaceTreeTraits<Dimension> > {
       public:
            typedef SurfaceTreeTraits<Dimension> Traits;
        
            SurfaceLibMeshTree() {};
        
        
            static cutk::shared_ptr<SurfaceLibMeshTree> New(const cutk::shared_ptr<cutlibpp::Predicate> &predicate,
                const int maxElementsXNode=cutlibpp::DEFAULT_REFINE_MAX_ELEMENTS,
                const int maxDepth=cutlibpp::DEFAULT_REFINE_DEPTH
                                                 ) {
                using namespace cutlibpp;
                cutk::shared_ptr<SurfaceLibMeshTree> tree = cutk::make_shared<SurfaceLibMeshTree>();
                cutk::shared_ptr<NTreeWithTagsMutatorFactory < SurfaceLibMeshTree> > factory =
                cutk::make_shared<NTreeWithTagsMutatorFactory < SurfaceLibMeshTree> >(predicate);
                factory->setRefineParams(maxElementsXNode, maxDepth);
                std::cout<<"LibMeshTree::maxElements = "<<maxElementsXNode<<std::endl;
                std::cout<<"LibMeshTree::maxDepth = "<<maxDepth<<std::endl;
                tree->setMutatorFactory(factory);
                return tree;
            }
        
        
        };


		class ElementDofMap : public cutk::Serializable {
		public:
			void read(cutk::InputStream &is) override
			{
				int n;
				is >> n;
				global.resize(n);
				is.read(&global[0], n);
			}

			void write(cutk::OutputStream &os) const override
			{
				int n = global.size();
				os << n;
				os.write(&global[0], n);
			}


			inline bool empty() const
			{
				return global.empty();
			}

            std::vector<long> global;


		};
    
   // add new SurfaceElementAdapter
    
 }

#endif //LIBMESH_CUTLIBPP_ADAPTERS_HPP
