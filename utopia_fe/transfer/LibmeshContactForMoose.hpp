#ifndef LibmeshContactForMoose_HPP
#define LibmeshContactForMoose_HPP

#include "cutlibpp.hpp"
#include "cutlibpp_Base.hpp"
#include "cutlibpp_Tree.hpp"
#include "cutlibpp_NTreeMutatorFactory.hpp"
#include "cutlibpp_NTreeWithSpanMutatorFactory.hpp"
#include "cutlibpp_NTreeWithTagsMutatorFactory.hpp"
#include "cutlibpp_API.hpp"

#include "libmesh/mesh_inserter_iterator.h"
#include "libmesh/elem.h"
#include "libmesh/transient_system.h"
#include "libmesh/fe.h"

#include "express_Profiler.hpp"
#include "express_Redistribute.hpp"

#include "Box.hpp"
#include "MapSparseMatrix.hpp"
#include "utopia_fe.hpp"
#include "MortarAssemble.hpp"
#include "libmesh/serial_mesh.h"
#include "MortarAssemble.hpp"
#include "MortarAssembler.hpp"
#include <cmath>
#include <queue>


namespace utopia {
    using namespace libMesh;

    template<typename T>
    static void normalize(std::vector<T> &vec)
    {
        T len = 0.;
        for(uint d = 0; d < vec.size(); ++d) {
            len += vec[d] * vec[d];
        }
        
        len = std::sqrt(len);
        
        assert(len > 0);
        
        for(uint d = 0; d < vec.size(); ++d) {
            vec[d] /= len;
        }
    }
   
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
        
        ElementAdapter(MeshBase &fe, const libMesh::dof_id_type &element, const long element_handle, const int tag)
        : fe_(&fe), element_(element), element_handle_(element_handle), tag_(tag), dof_map_(nullptr)
        {
            assert(element < fe.n_elem());
            
            libMesh::Elem * e = fe.elem(element);
            
            for (libMesh::dof_id_type i = 0; i < e->n_nodes(); ++i) {
                const libMesh::Point &p = fe.node(e->node(i));
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
            
            assert(element_ < fe_->n_local_elem());
            
            return fe_->elem(element_);
        }
        
        const libMesh::Elem * get() const
        {
            assert(fe_);
            return fe_->elem(element_);
        }
        
        
        inline const  MeshBase  &space() const
        {
            assert(fe_);
            return *fe_;
        }
        
        
        void set_dof_map(std::vector<long> * ptr)
        {
            
            dof_map_ = ptr;
        }
        
        void set_face_id(std::vector<long> * ptr2)
        {
            
            face_id_ = ptr2;
        }
        
        
        inline const std::vector<long> &dof_map() const
        {
            assert(dof_map_);
            return *dof_map_;
        }
        
        
    private:
        MeshBase * fe_;
        libMesh::dof_id_type element_;
        long element_handle_;
        int tag_;
        BoxBoxAdapter<Dimension> bound_;
        std::vector<long> * dof_map_;
        std::vector<long> * face_id_;
        
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
        
        SurfaceElementAdapter(MeshBase &fe, const libMesh::dof_id_type &element, const long element_handle, const int tag, /*std::vector<long> &map,*/ const libMesh::Real blow_up)
        : fe_(&fe), element_(element), element_handle_(element_handle), tag_(tag), dof_map_(nullptr), face_id_(nullptr)
        {
            assert(element < fe.n_elem());
            
            libMesh::Elem * e = fe.elem(element);
            
            const int dim = fe.mesh_dimension();
            
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
        :fe_(nullptr) , element_(-1), element_handle_(-1), tag_(-1), dof_map_(nullptr) , face_id_(nullptr) {}
        
        
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
            
            assert(element_ < fe_->n_local_elem());
            
            return fe_->elem(element_);
        }
        
        const libMesh::Elem * get() const
        {
            assert(fe_);
            return fe_->elem(element_);
        }
        
        
        inline const  MeshBase  &space() const
        {
            assert(fe_);
            return *fe_;
        }
        
        
        void set_dof_map(std::vector<long> * ptr)
        {
            
            dof_map_ = ptr;
        }
        
        void set_face_id(std::vector<long> * ptr2)
        {
            
            face_id_ = ptr2;
        }
        
        inline const std::vector<long> &dof_map() const
        {
            assert(dof_map_);
            return *dof_map_;
        }
        
        inline const std::vector<long> &dof_map_face() const
        {
            assert(face_id_);
            return *face_id_;
        }
        
        
    private:
        MeshBase * fe_;
        libMesh::dof_id_type element_;
        long element_handle_;
        int tag_;
        BoxBoxAdapter<Dimension> bound_;
        std::vector<long> * dof_map_;
        std::vector<long> * face_id_;
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
    
    
    class UtopiaMesh {
    public:
        
        
        explicit UtopiaMesh(const express::Communicator &comm) : comm(comm)
        {
            must_destroy_attached[0] = false;
        }
        
        UtopiaMesh(const std::shared_ptr<MeshBase> &master_slave,
                   const std::shared_ptr<libMesh::DofMap>  &original_dofmap,
                   const std::shared_ptr<const unsigned int> &_var_num)
        {
            
            utopiamesh_.reserve(1);
            
            utopiamesh_.push_back(master_slave);
            
            must_destroy_attached[0] = false;
            
            express::Communicator comm = master_slave->comm().get();
            
            const int n_elements = master_slave->n_elem();
            
            copy_global_dofs(*master_slave, original_dofmap, _var_num,
                             dof_maps_[0], var_type_[0], n_elements, var_number_[0],
                             subdomain_id_[0], side_set_id_[0], side_set_id_tag_[0],
                             face_set_id_global_[0],ownershipRangesFaceID_[0], 101, 102);
            
//            copy_var_number(*original_dofmap, var_number_[0]);
            
            copy_var_order(*original_dofmap, var_order_[0]);
            
        }
        
        
        
        inline std::vector< std::shared_ptr<MeshBase> > &utopiamesh()
        {
            return utopiamesh_;
            
        }
        
        inline const std::vector< std::shared_ptr<MeshBase> > &utopiamesh() const
        {
            return utopiamesh_;
            
        }
        
        
        inline long n_elements() const
        {
            
            long ret=0;
            ret += utopiamesh_[0]->n_elem();
            return ret;
            
        }
        
        
        inline std::vector<ElementDofMap> &dof_map()
        {
            
            return dof_maps_[0];
        }
        
        inline const std::vector<ElementDofMap> &dof_map() const
        {
            
            return dof_maps_[0];
        }
        
        inline void set_must_destroy_attached(const int index, const bool value)
        {
            
            must_destroy_attached[index] = value;
        }
        
        
        inline  std::vector<ElementDofMap> &variable_number()
        {
            
            return var_number_[0];
        }
        
        inline const std::vector<ElementDofMap> &variable_number() const
        {
            
            return var_number_[0];
        }
        
        
        
        inline std::vector<ElementDofMap> &variable_order()
        {
            
            return var_order_[0];
        }
        
        
        
        inline const std::vector<ElementDofMap> &variable_order() const
        {
            
            return var_order_[0];
        }
        
        
        inline std::vector<ElementDofMap> &variable_type()
        {
            
            return var_type_[0];
        }
        
        
        
        inline const std::vector<ElementDofMap> &variable_type() const
        {
            return var_type_[0];
        }
        
        inline std::vector<ElementDofMap> &subdomain_id()
        {
            
            return subdomain_id_[0];
        }
        
        
        
        inline const std::vector<ElementDofMap> &subdomain_id() const
        {
            return subdomain_id_[0];
        }
        
        
        inline std::vector<ElementDofMap> &side_set_id()
        {
            
            return side_set_id_[0];
        }
        
        
        
        inline const std::vector<ElementDofMap> & side_set_id() const
        {
            return side_set_id_[0];
        }
        
        
        inline std::vector<ElementDofMap> &face_set_id_global()
        {
            
            return face_set_id_global_[0];
        }
        
        
        
        inline const std::vector<ElementDofMap> & face_set_id_global() const
        {
            return face_set_id_global_[0];
        }
        
        inline std::vector<ElementDofMap> & n_face_nodes()
        {
            
            return n_face_nodes_[0];
        }
        
        
        
        inline const std::vector<ElementDofMap> & n_face_nodes() const
        {
            return n_face_nodes_[0];
        }
        
        
        inline std::vector<ElementDofMap> &side_set_id_tag()
        {
            
            return side_set_id_tag_[0];
        }
        
        
        
        inline const std::vector<ElementDofMap> & side_set_id_tag() const
        {
            return side_set_id_tag_[0];
        }
        
        inline express::Array<express::SizeType> & ownershipRangesFaceID()
        {
            
            return ownershipRangesFaceID_[0];
        }
        
        
        
        inline const express::Array<express::SizeType>   & ownershipRangesFaceID() const
        {
            return ownershipRangesFaceID_[0];
        }
        
        
    private:
        express::Communicator comm;
        std::vector<std::shared_ptr< MeshBase>>  utopiamesh_;
        std::vector<ElementDofMap> dof_maps_[1];
        std::vector<ElementDofMap> var_number_[1];
        std::vector<ElementDofMap> var_order_[1];
        std::vector<ElementDofMap> var_type_[1];
        std::vector<ElementDofMap> subdomain_id_[1];
        std::vector<ElementDofMap> side_set_id_[1];
        std::vector<ElementDofMap> face_set_id_global_[1];
        std::vector<ElementDofMap> side_set_id_tag_[1];
        std::vector<ElementDofMap> n_face_nodes_[1];
        express::Array<express::SizeType> ownershipRangesFaceID_[1];
        bool must_destroy_attached[1];
        
        
        
        
        inline static void copy_global_dofs(MeshBase &space, const std::shared_ptr<libMesh::DofMap>  &original_dofmap,
                                            const std::shared_ptr<const unsigned int>  &var_num,
                                            std::vector<ElementDofMap> &dof_map, std::vector<ElementDofMap> &variable_type,
                                            const int n_elements, std::vector<ElementDofMap> &variable_number,
                                            std::vector<ElementDofMap> &subdomain_id, std::vector<ElementDofMap> &side_set_id,
                                            std::vector<ElementDofMap> &side_set_id_tag, std::vector<ElementDofMap> &face_set_id_global,
                                            express::Array<express::SizeType>  & ownershipRangesFaceID, int tag_1, int tag_2)
        {
            
            auto &mesh = space;
//          auto &original_dof_map = space.dof_map();
            std::vector<dof_id_type> temp;
            
            
            
            
            express::Communicator comm = mesh.comm().get();
            // express::Array<express::SizeType>
            ownershipRangesFaceID.resize(comm.size()+1);
            ownershipRangesFaceID.allSet(0);
            std::vector<ElementDofMap> face_set_id;
            dof_map.resize(n_elements);
            subdomain_id.resize(n_elements);
            side_set_id.resize(n_elements);
            side_set_id_tag.resize(n_elements);
            face_set_id.resize(n_elements);
            face_set_id_global.resize(n_elements);
            variable_number.resize(1);
            //n_face_nodes.resize(n_elements);
            
            
            
            variable_type.resize(1);
            
            bool first=true;
            
            std::vector<const Node *> elem_nodes;
            int jj_side_id_one = 0;
            int jj_side_id_one_tag = 0;
            int jj_side_id_one_check = 0;
            int offset=0;
            int f_id=0;
            int n_f=0;
            MeshBase::const_element_iterator e_it_s = mesh.active_local_elements_begin();
            const MeshBase::const_element_iterator e_end_s = mesh.active_local_elements_end();
            
            for (; e_it_s != e_end_s; ++e_it_s)
            {
                Elem * elem = *e_it_s;
                
                
                bool  check_side_id_one=true;
                bool  check_side_id_one_tag=true;
                bool  check_side_id_one_check=true;
                bool  check_face_id=true;
                
                for (int side_elem=0; side_elem<elem->n_sides(); side_elem++){
                    if (check_side_id_one==true){
                        side_set_id[elem->id()].global.insert(side_set_id[elem->id()].global.end(),-1);
                        check_side_id_one=false;
                        jj_side_id_one++;
                    }
                }
                
                if (elem->on_boundary()){
                    for (int side_elem=0; side_elem<elem->n_sides(); side_elem++){
                        if ((mesh.get_boundary_info().has_boundary_id(elem,side_elem,tag_1) || mesh.get_boundary_info().has_boundary_id(elem,side_elem,tag_2)) && check_side_id_one_tag==true){
                            side_set_id[elem->id()].global.insert(side_set_id[elem->id()].global.end()-1,mesh.get_boundary_info().boundary_id(elem,side_elem));
                            check_side_id_one_tag=false;
                            jj_side_id_one_tag++;
                        }
                    }
                }
                
                
                
                for (int side_elem=0; side_elem<elem->n_sides(); side_elem++){
                    if (check_side_id_one_check==true){
                        //                        std::cout<<"side_set_id[ "<< elem->id() <<" ] = "<<side_set_id[elem->id()].global.at(0)<<std::endl;
                        check_side_id_one_check=false;
                        jj_side_id_one_check++;
                    }
                }
                
                if (elem->on_boundary()){
                    
                    for(uint side_elem = 0; side_elem < elem->n_sides(); ++side_elem){
                        //                        n_face_nodes[elem->id()].global.insert(n_face_nodes[elem->id()].global.end(), n_f);
                        //                        n_f++;
                        if ((mesh.get_boundary_info().has_boundary_id(elem,side_elem,tag_1) || mesh.get_boundary_info().has_boundary_id(elem,side_elem,tag_2))){
                            face_set_id[elem->id()].global.insert(face_set_id[elem->id()].global.end(),f_id);
                            //std::cout<<"f_id"<<f_id<<std::endl;
                            f_id++;
                            offset++;
                            
                        }
                    }
                }
                else
                {
                    
                    face_set_id[elem->id()].global.insert(face_set_id[elem->id()].global.end(),-1);
                }
                
                
                subdomain_id[elem->id()].global.insert(subdomain_id[elem->id()].global.end(),elem->subdomain_id());
                original_dofmap->dof_indices(elem, temp, *var_num);
                //std::cout<<" TEMP SIZE "<< temp.size() <<std::endl;
                dof_map[elem->id()].global.insert(dof_map[elem->id()].global.end(), temp.begin(), temp.end());
                
                if (first)
                {
                    variable_type[0].global.push_back(elem->type());
                    variable_number[0].global.push_back(*var_num);
                    first=false;


                    
                }
                
            }
            
            
            std::cout<<"size for all"<<face_set_id.size()<<std::endl;
            
            
            ownershipRangesFaceID[comm.rank()+1]+= static_cast<unsigned int>(offset);
            
            //std::cout<<"comm.rank"<<comm.rank()<<std::endl;
            
            comm.allReduce(&ownershipRangesFaceID[0],  ownershipRangesFaceID.size(),  express::MPISum());
            
            std::partial_sum(ownershipRangesFaceID.begin(), ownershipRangesFaceID.end(),
                             ownershipRangesFaceID.begin());
            
            if(comm.isRoot()) {
                std::cout << "ownershipRangesFaceID = "<< ownershipRangesFaceID << std::endl;
            }
            
            
            MeshBase::const_element_iterator e_it_new = mesh.active_local_elements_begin();
            const MeshBase::const_element_iterator e_end_new = mesh.active_local_elements_end();
            
            
            for (; e_it_new != e_end_new; ++e_it_new)
            {
                Elem * elem_new = *e_it_new;
                
                if (elem_new->on_boundary()){
                    //                 std::cout<<"jj<face_set_id[elem_new->id()].global.size()"<<face_set_id[elem_new->id()].global.size()<<std::endl;
                    //                 std::cout<<"face_id"<<face_id.size()<<std::endl;
                    for (int jj=0; jj<face_set_id[elem_new->id()].global.size(); jj++){
                        int i = face_set_id[elem_new->id()].global.at(jj);
                        //std::cout<<"size for el"<<n_face_nodes[elem_new->id()].global.size()<<std::endl;
                        
                        if(i!=-1){
                            //std::cout<<"comm.rank"<<ownershipRangesFaceID[comm.rank()]<<std::endl;
                            int global_id = i + ownershipRangesFaceID[comm.rank()] ;
                            //std::cout<<"global_id"<<i + ownershipRangesFaceID[comm.rank()] <<std::endl;
                            face_set_id_global[elem_new->id()].global.insert(face_set_id_global[elem_new->id()].global.end(),global_id);
                        }
                        else
                            face_set_id_global[elem_new->id()].global.insert(face_set_id_global[elem_new->id()].global.end(),-1);
                    }
                }
            }
        }
        
    
        
        inline static void copy_var_order(DofMap &dofmap, std::vector<ElementDofMap> &variable_order)
        {
            variable_order.resize(1);
            FEType fe_type =  dofmap.variable(0).type();
            variable_order[0].global.push_back(fe_type.order);
        }
        
        
    };

    
    
    template<class Iterator>
    static void write_space(const Iterator &begin, const Iterator &end,MeshBase &space, const std::vector<ElementDofMap> &dof_map, const std::vector<ElementDofMap> &variable_number, const std::vector<ElementDofMap> &variable_order, const std::vector<ElementDofMap> &subdomain_id, const std::vector<ElementDofMap> &side_set_id, const std::vector<ElementDofMap> &face_set_id_global,cutk::OutputStream &os, const int tag_1, const int tag_2)
    {
        const int dim 		  = space.mesh_dimension();
        const long n_elements = std::distance(begin, end);
        
        std::set<long> nodeIds;
        std::map<long, long> mapping;
        std::vector<dof_id_type> dof_array;
        
        for(Iterator it = begin; it != end; ++it) {
            
            const Elem *elem = space.elem(*it);
            
            for(dof_id_type j = 0; j != elem->n_nodes(); ++j) {
                
                nodeIds.insert(elem->node(j));
                
                
            }
        }
        
        long n_nodes = nodeIds.size();
        
        // Estimate for allocation
        os.requestSpace( (n_elements * 8 + n_nodes * dim) * (sizeof(double) + sizeof(long)) );
        
        //WRITE 1
        os << dim;
        
        //std::cout<<"dim write= "<<dim<<std::endl;
        
        int index = 0;
        for (auto nodeId : nodeIds) {
            mapping[nodeId] = index++;
        }
        
        //WRITE 2
        os << n_nodes;
        
        //WRITE 6
        os << n_elements;
        
        //        std::cout<<"elem write= "<<n_elements<<std::endl;
        
        
        
        for(auto node_id : nodeIds){
            
            const Point &p = space.node(node_id);
            
            for(int i = 0; i < dim; ++i) {
                
                //WRITE 3
                os << p(i);
                
            }
            
            //std::cout<<"write_point"<<p<<std::endl;
            
        }
        
        std::vector<dof_id_type> indices_vector;
        
        
        
        for(Iterator it = begin; it != end; ++it) {
            
            const int k = *it;
            
            const Elem *elem = space.elem(*it);
            
            const int e_n_nodes = elem->n_nodes();
            
            const int type = elem->type();
            
            //WRITE 7
            os << type << e_n_nodes;
            
            
            
            for (int i = 0; i != e_n_nodes; ++i) {
                
                auto it = mapping.find(elem->node(i));
                
                assert(it != mapping.end());
                
                int index = it->second;
                
                //WRITE 8
                os << index;
                
            }
            
            
            
            //WRITE 9
            assert(!dof_map.at(elem->id()).empty());
            
            os << dof_map.at(elem->id());
            
            bool  size=true;
            
            int volume_tag;
            
            volume_tag=subdomain_id[elem->id()].global.at(0);
            
            os << volume_tag;
            
            int side_set_tag;
            
            int face_id;
            
            bool check_side_id_one=true;
            
            //            for (int side_elem=0; side_elem<elem->n_sides(); side_elem++){
            //
            //                if(check_side_id_one==true){
            
            side_set_tag=side_set_id[elem->id()].global.at(0);
            
            
            //                    std::cout<<" write surface role outside= "<< side_set_tag <<std::endl;
            
            os << side_set_tag;
            
            
            // std::cout <<"write value"<< face_set_id_global[elem->id()].global.at(0)<<std::endl;
            
            os << face_set_id_global.at(elem->id());
            
            //            os << n_face_nodes.at(elem->id())
            //
            //                    check_side_id_one=false;
            //                }
            //            }
        }
        //
        //
        
        //WRITE 11
        os << variable_number.at(0);
        
        //WRITE 12
        os << variable_order.at(0);
        
    }

    template<class Iterator>
    static void write_element_selection(const Iterator &begin, const Iterator &end, const UtopiaMesh &utopiamesh, cutk::OutputStream &os)
    {
        
        
        auto m = utopiamesh.utopiamesh()[0];
        
        write_space(begin, end, *m, utopiamesh.dof_map(), utopiamesh.variable_number(), utopiamesh.variable_number(), utopiamesh.subdomain_id(), utopiamesh.side_set_id(), utopiamesh.face_set_id_global(), os, 101, 102);
        
        
        
        
        
    }
    
    
    static void read_space(cutk::InputStream &is, cutk::shared_ptr<MeshBase> & space,
                           std::vector<ElementDofMap> &dof_map,
                           std::vector<ElementDofMap> &variable_number,
                           std::vector<ElementDofMap> &variable_order,
                           std::vector<ElementDofMap> &subdomain_id,
                           std::vector<ElementDofMap> &side_set_id,
                           std::vector<ElementDofMap> &face_set_id_global,
                           const libMesh::Parallel::Communicator &comm,
                           int tag_1, int tag_2)
    {
        using namespace std;
        
        
        //READ 1
        int dim;
        is >> dim;
        std::cout<<"I am reading the mesh "<<std::endl;
        
        
        //READ 2
        long n_nodes;
        is >> n_nodes;
        
        //READ 6
        long n_elements;
        is >> n_elements;
        
        std::cout<<"elem read= "<<n_elements<<std::endl;
        auto mesh_ptr = std::make_shared<SerialMesh>(comm,dim);
        
        mesh_ptr->reserve_nodes(n_nodes);
        
        for (long iii = 0; iii != n_nodes; ++iii) {
            
            Point p;
            
            for(int j = 0; j < dim; ++j) {
                //READ 3
                is >> p(j);
            }
            
            mesh_ptr->add_point(p);
            //std::cout<<"read_point"<<p<<std::endl;
            
        }
        
        
        
        dof_map.resize(n_elements);
        
        subdomain_id.resize(n_elements);
        
        side_set_id.resize(n_elements);
        
        face_set_id_global.resize(n_elements);
        
        face_set_id_global.resize(n_elements);
        
        
        for(long i = 0; i !=n_elements; ++i) {
            
            //READ 7
            
            int type, e_n_nodes;
            
            is >> type >> e_n_nodes;
            
            //std::cout<<"e_n_nodes_read = "<<e_n_nodes<<std::endl;
            
            auto elem =  Elem::build(ElemType(type)).release();
            
            //std::cout<<"n_side_read ="<< elem->n_sides()<<std::endl;
            
            
            int index;
            
            for (int ii = 0; ii != e_n_nodes; ++ii) {
                
                //READ 8
                is >> index;
                //std::cout<<"index = "<<index<<std::endl;
                elem->set_node(ii) = & mesh_ptr->node(index);
                
            }
            
            
            //READ 9
            is >> dof_map.at(i);
            //std::cout<< "dof_map_read = "<<dof_map[i].global.at(0)<<std::endl;
            
            int volume_tag, side_set_tag, face_id;
            
            bool on_boundary=false;
            //std::cout<<"read n_elements = "<<n_elements<<std::endl;
            
            
            is >> volume_tag;
            
            //std::cout<<" read volume role = "<< volume_tag <<std::endl;
            
            subdomain_id[i].global.insert(subdomain_id[i].global.end(),volume_tag);
            
            is >> side_set_tag;
            
            is >> face_set_id_global.at(i);
            
            //std::cout <<"read value"<< face_set_id_global[i].global.at(0)<<std::endl;
            
            side_set_id[i].global.insert(side_set_id[i].global.end(),side_set_tag);
            
            mesh_ptr->add_elem(elem);
            
            libmesh_assert(elem);
            
        }
        
        //READ 11
        variable_number.resize(1);
        is >> variable_number.at(0);
        
        //READ 12
        variable_order.resize(1);
        is >> variable_order.at(0);
        
        
        
        
        
        //EquationSystems es(*mesh_ptr);
        // const System & new_sys = es.get_system(0);
        
        //!!!! dummy parameters
        space = mesh_ptr;
        
    }
    
    
    static void read_spaces(cutk::InputStream &is, UtopiaMesh &utopiamesh, const libMesh::Parallel::Communicator &comm_mesh)
    {
        
        bool has_master, has_slave;
        // is >> has_master >> has_slave;
        
        utopiamesh.utopiamesh().resize(1);
        
        read_space(is, utopiamesh.utopiamesh()[0], utopiamesh.dof_map(), utopiamesh.variable_number(), utopiamesh.variable_order(), utopiamesh.subdomain_id(), utopiamesh.side_set_id(), utopiamesh.face_set_id_global(), comm_mesh, 101, 102);
        
        utopiamesh.set_must_destroy_attached(0,true);
        
    }
    
    
    template<int Dimensions, class Fun>
    static bool SurfaceAssemble(express::Communicator &comm,
                                const std::shared_ptr<MeshBase> &master_slave,
                                const std::shared_ptr<DofMap> &dof_map,
                                const std::shared_ptr<const unsigned int> &_var_num,
                                Fun process_fun,
                                const cutk::Settings &settings,const libMesh::Real search_radius, const int tag_1, const int tag_2)
    {
        
        std::shared_ptr<UtopiaMesh> local_fun_spaces_new = cutk::make_shared<UtopiaMesh>(master_slave, dof_map, _var_num);
        auto predicate = std::make_shared<cutlibpp::MasterAndSlave>();
        predicate->add(tag_1,tag_2);
        
        using namespace cutlibpp;
        using namespace express;
        using namespace cutk;
        
        typedef SurfaceLibMeshTree<Dimensions> NTreeT;
        typedef typename NTreeT::DataContainer DataContainer;
        typedef typename NTreeT::DataType SurfaceAdapter;
        
        long maxNElements = 40;
        long maxDepth = 5;
        
        
        if (!settings.get("max_depth").isNull()) {
            maxDepth = settings.get("max_depth").toInt();
            //std::cout<<"max_depth  = "<< maxDepth  <<std::endl;
        }
        
        const auto &mesh = master_slave;
        
        const int n_elements = mesh->n_elem();
        
        std::cout << "mesh_elem_inside " << n_elements << std::endl;
        
        const Parallel::Communicator &libmesh_comm_mesh = master_slave->comm();
        
        
        
        const int dim_src = master_slave->mesh_dimension();
        const int dim_sla = master_slave->mesh_dimension();
        
        MeshBase::const_element_iterator e_it = mesh->active_elements_begin();
        const MeshBase::const_element_iterator e_end = mesh->active_elements_end();
        std::vector<int> block_id;
        std::vector<int> block_id_def;
        
        int i=0;
        for (; e_it != e_end; ++e_it)
        {
            Elem * elem = *e_it;
            if (i==0){
                block_id_def.push_back(elem->subdomain_id());}
            
            block_id.push_back(elem->subdomain_id());
            if (i>0 && block_id.at(i)!=block_id.at(i-1)){
                block_id_def.push_back(block_id.at(i));
            }
            
            i++;
            
            
        }
        
        
        EXPRESS_EVENT_BEGIN("create_adapters");
        ////////////////////////////////////////////////////////////////////////////////////////////////////
        cutk::shared_ptr<NTreeT> tree = NTreeT::New(predicate, maxNElements, maxDepth);
        
        tree->reserve(n_elements);
        
        std::cout << "nElements = tree->memory().nData()_inside " << n_elements << std::endl;
        
        std::shared_ptr<UtopiaMesh> local_spaces = make_shared<UtopiaMesh>(master_slave, dof_map, _var_num);
        
        int jj=0;
        
        for (auto it = master_slave->active_local_elements_begin(); it!=master_slave->active_local_elements_end(); ++it) {
            
            auto elem=*it;
            
            if(!elem->on_boundary()) {
                continue;
            }
            
            bool check_size=false;
            
            for(uint side_elem = 0; side_elem < elem->n_sides(); ++side_elem){
                if ((predicate->select(master_slave->get_boundary_info().boundary_id(elem, side_elem))) && check_size==false){
                    SurfaceAdapter a(*master_slave, elem->id(), elem->id(), master_slave->get_boundary_info().boundary_id(elem, side_elem), search_radius);
                    //                    std::cout<<"side_set_id[ "<< elem->id() <<" ] = "<<local_spaces->side_set_id()[elem->id()].global.at(0)<<std::endl;
                    assert(!local_spaces->dof_map()[elem->id()].empty());
                    a.set_dof_map(&local_spaces->dof_map()[elem->id()].global);
                    a.set_face_id(&local_spaces->face_set_id_global()[elem->id()].global);
                    tree->insert(a);
                    check_size=true;
                    //                    jj++;
                }
            }
        }
        
        
        //        std::cout<<"jj_tree_elem = "<< jj <<std::endl;
        
        
        tree->getRoot()->getBound().staticBound().enlarge(1e-8);
        
        
        
        
        
        
        ////////////////////////////////////////////////////////////////////////////////////////////////////
        EXPRESS_EVENT_END("create_adapters");
        
        //Just to have an indexed-storage
        std::map<long, cutk::shared_ptr<UtopiaMesh> > utopiamesh;
        std::map<long, std::vector<cutk::shared_ptr<UtopiaMesh> > > migrated_meshes;
        
        
        auto read = [&utopiamesh, &migrated_meshes, block_id, comm, &libmesh_comm_mesh, search_radius]
        (
         const long ownerrank,
         const long senderrank,
         bool is_forwarding, DataContainer &data,
         InputStream &in
         ) {
            
            CHECK_STREAM_READ_BEGIN("vol_proj", in);
            
            cutk::shared_ptr<UtopiaMesh> proc_space = cutk::make_shared<UtopiaMesh>(comm);
            
            
            //std::cout<<"I am in read space"<<std::endl;
            
            read_spaces(in, *proc_space, libmesh_comm_mesh);
            
            if (!is_forwarding) {
                assert(!utopiamesh[ownerrank]);
                utopiamesh[ownerrank] = proc_space;
            } else {
                migrated_meshes[ownerrank].push_back(proc_space);
            }
            
            data.reserve(data.size() + proc_space->n_elements());
            
            for(auto s : proc_space->utopiamesh()){
                int i=0;
                for (int i = 0; i<s->n_elem(); ++i) {
                    auto elem=s->elem(i);
                    int tag =proc_space->side_set_id()[i].global.at(0);
                    data.push_back(SurfaceAdapter(*s, i, i,tag,search_radius));
                    assert(!proc_space->dof_map()[i].empty());
                    assert(!proc_space->side_set_id()[i].empty());
                    data.back().set_dof_map(&proc_space->dof_map()[i].global);
                    data.back().set_face_id(&proc_space->face_set_id_global()[i].global);
                    //std::cout<<"&proc_space->dof_map()[i].global"<<proc_space->face_set_id_global()[i].global.at(0)<<std::endl;
                    
                    
                }
            }
            
            
            
            
            CHECK_STREAM_READ_END("vol_proj", in);
            
            
            
        };
        
        
        
        auto write = [&local_spaces, &utopiamesh, &comm]
        (
         const long ownerrank, const long recvrank,
         const std::vector<long>::const_iterator &begin,
         const std::vector<long>::const_iterator &end,
         const DataContainer &data,
         OutputStream &out) {
            
            CHECK_STREAM_WRITE_BEGIN("vol_proj", out);
            
            
            
            if (ownerrank == comm.rank()) {
                write_element_selection(begin, end, *local_spaces, out);
                
                
            } else {
                auto it = utopiamesh.find(ownerrank);
                assert(it != utopiamesh.end());
                cutk::shared_ptr<UtopiaMesh> spaceptr = it->second;
                assert(std::distance(begin, end) > 0);
                write_element_selection(begin, end, *spaceptr, out);
                
            }
            
            
            
            CHECK_STREAM_WRITE_END("vol_proj", out);
            
        };
        
        
        
        long n_false_positives = 0, n_projections = 0;
        
        
        
        auto fun = [&n_false_positives, &n_projections, &process_fun](
                                                                      SurfaceAdapter &master, SurfaceAdapter &slave) -> bool {
            //std::cout<<"n_intersections"<<n_intersections<<std::endl;
            bool ok = process_fun(master, slave);
            
            if(ok) {
                n_projections++;
                //std::cout<<"n_intersections"<<n_intersections<<std::endl;
                return true;
            } else {
                n_false_positives++;
                return false;
            }
            return true;
            
        };
        

        
        cutk::Settings custom_settings = settings;
        custom_settings.set("disable_redistribution", cutk::Boolean(true));
        custom_settings.set("verbosity_level", cutk::Integer(2));
        
        cutlibpp::search_and_compute(comm, tree, predicate, read, write, fun, custom_settings);
        
        long n_total_candidates = n_projections + n_false_positives;
        
        //std::cout<<"n_total_candidates"<<n_intersections<<std::endl;
        
        long n_collection[3] = {n_projections, n_total_candidates, n_false_positives};
        comm.allReduce(n_collection, 3, express::MPISum());
        
        if (comm.isRoot()) {
            std::cout << "n_intersections: " << n_collection[0]
            << ", n_total_candidates: " 	 << n_collection[1]
            << ", n_false_positives: " 	     << n_collection[2] << std::endl;
        }
        
        return true;
    }
    
    
    template<int Dimensions>
    bool SurfaceAssemble(
                         express::Communicator &comm,
                         const std::shared_ptr<MeshBase> &master_slave,
                         const std::shared_ptr<DofMap> &dof_map,
                         const std::shared_ptr<const unsigned int> &_var_num,
                         DSMatrixd &B,
                         DSMatrixd &orthogonal_trafos,
                         DVectord &gap,
                         DSMatrixd &normals,
                         DVectord &is_contact_node,
                         const cutk::Settings &settings,const libMesh::Real search_radius,
                         const int tag_1, const int tag_2)
    {
        std::shared_ptr<UtopiaMesh> local_fun_spaces_new = cutk::make_shared<UtopiaMesh>(master_slave, dof_map, _var_num);
        
        libMesh::DenseMatrix<libMesh::Real> src_pts;
        libMesh::DenseMatrix<libMesh::Real> dest_pts;
        libMesh::DenseMatrix<libMesh::Real> intersection2;
        Polyhedron src_poly, dest_poly;
        Polyhedron  intersection3,temp_poly;
        Intersector isector;
        
        std::shared_ptr<MeshBase> master_slave_space = master_slave;
        
        auto predicate = std::make_shared<cutlibpp::MasterAndSlave>();
        predicate->add(tag_1,tag_2);
        std::cout<<"tag value 1 = "<<tag_1<<std::endl;
        std::cout<<"tag value 2 = "<<tag_2<<std::endl;
        
        static const double tol = 1e-8;
        
        
        std::vector<libMesh::dof_id_type> master_dofs, slave_dofs;
        libMesh::DenseMatrix<libMesh::Real> elemmat;
        libMesh::DenseMatrix<libMesh::Real> cumulative_elemmat;
        DenseMatrix<Real> side_polygon_1, side_polygon_2;
        DenseMatrix<Real> isect_polygon_1, isect_polygon_2;
        
        std::shared_ptr<Transform> src_trans;
        std::shared_ptr<Transform> dest_trans;
        
        Point n1, n2;
        
        
        
        int skip_zeros = 1;
        
        int slave_dof_n=0;
        
        int master_dof_n=0;
        
        const int dim = master_slave->mesh_dimension();
        
        
        libMesh::Real total_intersection_volume = 0.0;
        libMesh::Real local_element_matrices_sum = 0.0;
        
        MeshBase::const_element_iterator e_it = master_slave->active_local_elements_begin();
        
        const MeshBase::const_element_iterator e_end = master_slave->active_local_elements_end();
        
        
        for (; e_it != e_end; ++e_it)
        {
            Elem * elem = *e_it;
            
            // std::cout <<"subdomain_id = "<< elem->subdomain_id() << "elem_id = " << elem->id() <<std::endl;
            
            bool size_new=true;
            
            int jj = 0;
            
            if (elem->on_boundary()){
                
                for (int side_elem=0; side_elem<elem->n_sides(); side_elem++){
                    
                    if (master_slave->get_boundary_info().has_boundary_id(elem,side_elem, tag_1)){
                        
                        master_dof_n++;
                    }
                    
                    if (master_slave->get_boundary_info().has_boundary_id(elem,side_elem, tag_2)){
                        
                        slave_dof_n++;
                    }
                }
            }
        }
        //
        //        int mat_buffer_row = master_slave->mesh().n_nodes() * master_slave->mesh();
        //
        //        int mat_buffer_col = master_slave->dof_map().n_dofs() * master_slave->mesh().n_elem();
        //
        express::MapSparseMatrix<double> mat_buffer;
        
        express::MapSparseMatrix<double> p_buffer;
        
        express::MapSparseMatrix<double> q_buffer;
        
        express::MapSparseMatrix<double> rel_area_buff;
        
        express::MapSparseMatrix<double> gap_buff;
        
        express::MapSparseMatrix<double> normal_buff;
        
        
        
        std::cout<<"*********** master_slave->dof_map().n_dofs() = "<<  dof_map->n_dofs() <<std::endl;
        
        bool intersected = false;
        
        auto fun = [&](const SurfaceElementAdapter<Dimensions> &master,
                       const SurfaceElementAdapter<Dimensions> &slave) -> bool {
            
            long n_intersections = 0;
            
            using namespace cutlibpp;
            using namespace express;
            using namespace cutk;
            
            
            
            
            //std::cout<<"ciao sn in fun"<<std::endl;
            
            
            const auto &src  = master.space();
            const auto &dest = slave.space();
            
            
            
            
            const auto &src_mesh  = src;
            const auto &dest_mesh = dest;
            
            //dest_mesh.print_info();
            
            libMesh::DenseMatrix<libMesh::Real> elemmat;
            
            const int src_index  = master.element();
            const int dest_index = slave.element();
            
            auto &src_el  = *src_mesh.elem(src_index);
            auto &dest_el = *dest_mesh.elem(dest_index);
            
            const int dim_src = src_mesh.mesh_dimension();
            const int dim_sla = dest_mesh.mesh_dimension();
            
            Box box_1(dim_src), box_2(dim_sla);
            
            QMortar src_ir_ref(dim_src);
            QMortar src_ir(dim_src);
            QMortar dest_ir(dim_sla);
            QMortar dest_ir_ref(dim_sla);
            
            std::vector<long> src_order = local_fun_spaces_new->variable_order()[0].global;
            const int approx_order=src_order[0];
            
            std::shared_ptr<Contact> surface_assemble;
            
            const auto &master_face_id = master.dof_map_face();
            const auto &slave_face_id  = slave.dof_map_face();
            
            
            std::unique_ptr<libMesh::FEBase> master_fe, slave_fe;
            //
            master_fe = libMesh::FEBase::build(src_mesh.mesh_dimension(), FIRST);
            slave_fe  = libMesh::FEBase::build(dest_mesh.mesh_dimension(), FIRST);
            //
            
            typedef Intersector::Scalar Scalar;
            
            
            if(dim_src == 2)  {
                make_polygon(src_el,   src_pts);
                make_polygon(dest_el,  dest_pts);
                src_trans  = std::make_shared<Transform2>(src_el);
                dest_trans = std::make_shared<Transform2>(dest_el);
                
            }
            
            else if(dim_src == 3) {
                make_polyhedron(src_el,  src_poly);
                make_polyhedron(dest_el, dest_poly);
                src_trans  = std::make_shared<Transform3>(src_el);
                dest_trans = std::make_shared<Transform3>(dest_el);
                
            }
            
            bool intersected = false;
            
            for(uint side_1 = 0; side_1 < src_el.n_sides(); ++side_1) {
                
                
                
                //                libMesh::UniquePtr<libMesh::Elem> s_1 = src_el.side(side_1);
                
                //std::cout << std::to_string(master_face_id[0]) << " -> " <<side_ptr_1->on_boundary() << std::endl;
                //if(src_el.neighbor_ptr(side_1)!=nullptr) continue;
                
                if(src_el.neighbor_ptr(side_1) != nullptr) continue;
                
                auto side_ptr_1 = src_el.build_side_ptr(side_1);
                compute_side_normal(dim_src, *side_ptr_1, n1);
                
                box_1.reset();
                enlarge_box_from_side(dim_src, *side_ptr_1, box_1, search_radius);
                
                if(dim_src == 2) {
                    make_polygon(*side_ptr_1, side_polygon_1);
                } else if(dim_src == 3) {
                    make_polygon_3(*side_ptr_1, side_polygon_1);
                } else {
                    assert(false);
                }
                
                
                
                
                for(uint side_2 = 0; side_2 < dest_el.n_sides(); ++side_2) {
                    
                    if(dest_el.neighbor_ptr(side_2) != nullptr) continue;
                    if (!predicate->tagsAreRelated(tag_1, tag_2)) continue;
                    
                    //                    libMesh::UniquePtr<libMesh::Elem> s_2 = dest_el.side(side_2);
                    auto side_ptr_2 = dest_el.build_side_ptr(side_2);
                    //[if(!side_ptr_2->on_boundary()) continue;
                    
                    compute_side_normal(dim_sla, *side_ptr_2, n2);
                    
                    const Real cos_angle = n1.contract(n2);
                    
                    //if the angle is more than 60 degrees ( cos(60/180*pi) == 0.5 ) or has same orientation skip
                    if(cos_angle >= -0.5) {
                        continue;
                    }
                    
                    
                    
                    
                    box_2.reset();
                    enlarge_box_from_side(dim_sla, *side_ptr_2, box_2, search_radius);
                    
                    if(!box_1.intersects(box_2, tol)) {
                        continue;
                    }
                    
                    
                    bool pair_intersected = false;
                    if(dim_sla==2){
                        make_polygon(*side_ptr_2, side_polygon_2);
                        
                        //plot_lines(2, 2, &side_polygon_1.get_values()[0], "in_master/" + std::to_string(master_face_id[0]) + "_" + std::to_string(cos_angle));
                        //plot_lines(2, 2, &side_polygon_2.get_values()[0], "in_slave/" + std::to_string(slave_face_id[0]) + "_" + std::to_string(cos_angle));
                        
                        if(!project_2D(side_polygon_1, side_polygon_2, isect_polygon_1, isect_polygon_2)){
                            continue;
                        }
                        const Scalar dx = side_polygon_2(0, 0) - side_polygon_2(1, 0);
                        const Scalar dy = side_polygon_2(0, 1) - side_polygon_2(1, 1);
                        
                        const Scalar isect_dx = isect_polygon_2(0, 0) - isect_polygon_2(1, 0);
                        const Scalar isect_dy = isect_polygon_2(0, 1) - isect_polygon_2(1, 1);
                        
                        const Scalar area   = std::sqrt(isect_dx*isect_dx + isect_dy*isect_dy);
                        const Scalar weight = area/std::sqrt(dx*dx + dy*dy);
                        
                        const int order = order_for_l2_integral(dim_src, src_el, approx_order, dest_el, approx_order);
                        
                        make_composite_quadrature_on_surf_2D(isect_polygon_1, weight, order, src_ir);
                        
                        make_composite_quadrature_on_surf_2D(isect_polygon_2, weight, order, dest_ir);
                        
                        pair_intersected = true;
                        
                        surface_assemble = std::make_shared<Contact>();
                        surface_assemble->isect_area	   = area;
                        surface_assemble->relative_area    = weight;
                        
                        
                        
                    } else if(dim_src == 3) {
                        make_polygon_3(*side_ptr_2, side_polygon_2);
                        
                        if(!project_3D(
                                       side_polygon_1,
                                       side_polygon_2,
                                       isect_polygon_1,
                                       isect_polygon_2))
                        {
                            continue;
                        }
                        
                        const Scalar area_slave = isector.polygon_area_3(side_polygon_2.m(),  &side_polygon_2.get_values()[0]);
                        const Scalar area   	= isector.polygon_area_3(isect_polygon_2.m(), &isect_polygon_2.get_values()[0]);
                        const Scalar weight 	= area/area_slave;
                        
                        const int order = order_for_l2_integral(dim_src, src_el, approx_order, dest_el, approx_order);
                        
                        make_composite_quadrature_on_surf_3D(isect_polygon_1, weight, order, src_ir);
                        make_composite_quadrature_on_surf_3D(isect_polygon_2, weight, order, dest_ir);
                        
                        pair_intersected = true;
                        
                        surface_assemble = std::make_shared<Contact>();
                        surface_assemble->isect_area	= area;
                        surface_assemble->relative_area = weight;
                        
                    } else {
                        assert(false);
                        return false;
                    }
                    
                    
                    if(pair_intersected) {
                        
                        //////////////////////////////////ASSEMBLY ////////////////////////////////////////
                        //////////////////////////////////////////////////////////////////////////////////////
                        transform_to_reference_surf(*src_trans,  src_el.type(),  src_ir, src_ir_ref);
                        transform_to_reference_surf(*dest_trans, dest_el.type(), dest_ir, dest_ir_ref);
                        
                        master_fe->attach_quadrature_rule(&src_ir_ref);
                        
                        master_fe->reinit(&src_el);
                        
                        slave_fe->attach_quadrature_rule(&dest_ir_ref);
                        
                        slave_fe->get_xyz();
                        
                        slave_fe->reinit(&dest_el);
                        
                        
                        surface_assemble->parent_element_master  = src_index;
                        
                        surface_assemble->id_master 			 = src_el.id();
                        
                        surface_assemble->parent_element_slave   = dest_index;
                        
                        surface_assemble->id_slave 			     = dest_el.id();
                        
                        surface_assemble->coupling.zero();
                        
                        elemmat.zero();
                        
                        mortar_assemble(*master_fe, *slave_fe, elemmat);
                        
                        
                        const libMesh::Point pp = side_ptr_1->point(0);
                        const Real plane_offset = n1.contract(pp);
                        
                        mortar_normal_and_gap_assemble(
                                                       dim,
                                                       *slave_fe,
                                                       n2,
                                                       n1,
                                                       plane_offset,
                                                       surface_assemble->normals,
                                                       surface_assemble->gap);
                        
                        
                        //////////////////////////////////////////////////////////////////////////////////////
                        
                        // use boundary info instead
                        rel_area_buff.add(slave_face_id[0], 0, surface_assemble->relative_area);
                        
                        
                        const auto &master_dofs = master.dof_map();
                        const auto &slave_dofs  = slave.dof_map();
                        
                        
                        
                        //
                        //                std::cout << "master_face_id.size()" << master_face_id.size() <<std::endl;
                        
                        
                        
                        MeshBase::const_element_iterator e_it = master_slave->active_elements_begin();
                        const MeshBase::const_element_iterator e_end = master_slave->active_elements_end();
                        
                        std::vector<int> block_id;
                        std::vector<int> block_id_def;
                        
                        
                        //
                        //                std::cout<<"************************************************** "<<std::endl;
                        //
                        //
                        //                std::cout<<"************* dest_el.id() = "<< dest_el.id() <<std::endl;
                        //
                        //
                        //                std::cout<<"************* src_el.id() = "<< src_el.id() <<std::endl;
                        //
                        //
                        //                std::cout<<"************************************************** "<<std::endl;
                        //
                        
                        
                        int i=0;
                        for (; e_it != e_end; ++e_it)
                        {
                            Elem * elem = *e_it;
                            if (i==0){
                                block_id_def.push_back(elem->subdomain_id());}
                            
                            block_id.push_back(elem->subdomain_id());
                            if (i>0 && block_id.at(i)!=block_id.at(i-1)){
                                block_id_def.push_back(block_id.at(i));
                            }
                            
                            i++;
                            
                            
                        }
                        
                        auto side_dest = dest_el.build_side_ptr(0);
                        
                        
                        int n_nodes_face_dest = side_dest->n_nodes();
                        
                        auto side_src = src_el.build_side_ptr(0);
                        
                        int n_nodes_face_src = side_src->n_nodes();
                        
                        //std::cout<<"addendum"<<addendum<<std::endl;
                        
                        std::vector<bool> node_is_boundary_dest(elemmat.m(), false);
                        std::vector<bool> node_is_boundary_src(elemmat.n(), false);
                        
                        
                        std::vector<double> sum_rows(elemmat.m(), 0.);
                        std::vector<double> sum_cols(elemmat.n(), 0.);
                        
                        
                        for(int i = 0; i < sum_rows.size(); ++i) {
                            for(int j = 0; j < sum_cols.size(); ++j) {
                                sum_rows[i] += std::abs(elemmat(i, j));
                                sum_cols[j] += std::abs(elemmat(i, j));
                            }
                        }
                        
                        
                        //Hacky but it is the only way with libmesh that we know of!!!!
                        for(int i = 0; i < sum_rows.size(); ++i) {
                            if(std::abs(sum_rows[i]) > 1e-16) {
                                node_is_boundary_dest[i] = true;
                            }
                        }
                        //                        int kk=0;
                        //                        for(int j = 0; j < dest_el.n_sides(); ++j){
                        //                            //auto side_dest=build_side_ptr(j);
                        //                            for(int i = 0; i < side_dest->n_nodes(); ++i) {
                        //                                if(dest_mesh.get_boundary_info().has_boundary_id(&dest_el,j,tag_2)) {
                        //                                    node_is_boundary_dest[kk] = true;
                        //                                    kk++;
                        //                                }
                        //                            }
                        //                        }
                        //
                        //                        std::cout<<"query nodes"<<kk<<std::endl;
                        //                        int kkk=0;
                        //                        for(int j = 0; j < src_el.n_sides(); ++j){
                        //                            //auto side_dest=build_side_ptr(j);
                        //                            for(int i = 0; i < side_src->n_nodes(); ++i) {
                        //                                if(src_mesh.get_boundary_info().has_boundary_id(&src_el,j,tag_1)) {
                        //                                    node_is_boundary_src[kkk] = true;
                        //                                    kkk++;
                        //                                }
                        //                            }
                        //                        }
                        
                        
                        //Hacky but it is the only way with libmesh that we know of!!!!
                        for(int i = 0; i < sum_cols.size(); ++i) {
                            if(std::abs(sum_cols[i]) > 1e-16) {
                                node_is_boundary_src[i] = true;
                            }
                        }
                        
                        
                        //                        build_boundary_query(dest_el, *side_dest, 1, node_is_boundary_dest);
                        //
                        
                        //                        build_boundary_query(src_el, *side_src, 1, node_is_boundary_src);
                        
                        //plot_lines(dim_src,  side_src->n_nodes(), &side_polygon_1.get_values()[0] , "master/" + std::to_string(master_face_id[0]));
                        //plot_lines(dim_sla, side_dest->n_nodes(), &side_polygon_2.get_values()[0] , "slave/" + std::to_string(slave_face_id[0]));
                        
                        
                        std::vector<dof_id_type> dof_indices_slave_vec(slave_dofs.size(),   -1);
                        std::vector<dof_id_type> dof_indices_master_vec(master_dofs.size(), -1);
                        
                        int offset = 0;
                        for(uint i = 0; i < node_is_boundary_dest.size(); ++i){
                            if (node_is_boundary_dest[i]) {
                                dof_indices_slave_vec[i] =  slave_face_id[0] * n_nodes_face_dest + offset++;
                            }
                        }
                        
                        offset = 0;
                        for(uint i = 0; i <  node_is_boundary_src.size(); ++i) {
                            if (node_is_boundary_src[i]) {
                                dof_indices_master_vec[i] = master_face_id[0] * n_nodes_face_src + offset++;
                            }
                        }
                        
                        
                        for(int i = 0; i <  slave_dofs.size(); ++i) {
                            const long dof_I = slave_dofs[i];
                            const long dof_J = dof_indices_slave_vec[i];
                            
                            if(node_is_boundary_dest[i]) {
                                p_buffer.setAt(dof_I, dof_J, 1.);
                            }
                        }
                        
                        for(int i = 0; i <  dof_indices_master_vec.size(); ++i) {
                            const long dof_I = master_dofs[i];
                            const long dof_J = dof_indices_master_vec[i];
                            
                            
                            
                            if(node_is_boundary_src[i]) {
                                q_buffer.setAt(dof_I, dof_J, 1.);
                            }
                        }
                        
                        
                        auto partial_sum = std::accumulate(elemmat.get_values().begin(), elemmat.get_values().end(), libMesh::Real(0.0));
                        
                        // const Scalar local_mat_sum = std::accumulate(surface_assemble->coupling.get_values().begin(), surface_assemble->coupling.get_values().end(), libMesh::Real(0.0));
                        
                        local_element_matrices_sum += partial_sum;
                        
                        assert(slave_dofs.size() == elemmat.m());
                        assert(master_dofs.size() == elemmat.n());
                        
                        
                        for(int i = 0; i <  dof_indices_slave_vec.size(); ++i) {
                            
                            if(!node_is_boundary_dest[i]) continue;
                            
                            const long dof_I = dof_indices_slave_vec[i];
                            
                            gap_buff.add(dof_I, 0, surface_assemble->gap(i));
                            
                            for (int k=0; k<dim_sla; ++k)
                            {
                                normal_buff.add(dof_I,k, surface_assemble->normals(i,k));
                            }
                            
                            //                            std::cout<< "************ dof_I_index = "<< dof_I <<std::endl;
                            
                            for(int j = 0; j <  dof_indices_master_vec.size(); ++j) {
                                if(!node_is_boundary_src[j]) continue;
                                
                                //                    std::cout<<"************* J = "<< j <<std::endl;
                                
                                const long dof_J = dof_indices_master_vec[j];
                                
                                //std::cout<< "************ dof_J_index = "<< dof_J <<std::endl;
                                
                                mat_buffer.add(dof_I, dof_J, elemmat(i, j));
                            }
                        }
                        
                        int dim = dim_src;
                        
                        //                std::cout<< "surface_assemble->isect_area = " << surface_assemble->isect_area <<std::endl;
                        //
                        //                std::cout<<" pow(surface_assemble->isect_area, dim/(dim-1.)) * dim = " << pow(surface_assemble->isect_area, dim/(dim-1.)) * dim  <<std::endl;
                        //
                        //                std::cout<<" partial sum = " << partial_sum <<std::endl;
                        //
                        //                std::cout<<" local_element_matrices_sum = " << local_element_matrices_sum <<std::endl;
                        //
                        //
                        
                        // assert(fabs(partial_sum - pow(surface_assemble->isect_area, dim/(dim-1.))) < 1e-8 || (!is_quad(dest_el.type()) && !is_hex(dest_el.type())));
                        intersected = true;
                        
                        //Assemble p
                        //Iteri su elementi
                        //calcoli i=dof_map.index(elemen), j=el_id * n_dofs_x_el + local_dof_offset (local_dof_id e.g., triangolo = {0, 1, 2}, quad = {0, 1, 2, 3})
                        
                    }
                }
            }
            
            return intersected;
            
        };
        
        
        
        if(!SurfaceAssemble<Dimensions>(comm, master_slave, dof_map, _var_num, fun, settings, search_radius, tag_1, tag_2)) {
            return false;
        }
        
        
        // std::cout << mat_buffer << std::endl;
        
        double volumes[1] = { local_element_matrices_sum };
        
        comm.allReduce(volumes, 1, express::MPISum());
        
        const processor_id_type master_proc_id  = master_slave->processor_id();
        
        const dof_id_type n_dofs_on_proc_master = dof_map->n_dofs_on_processor(master_proc_id);
        
        const processor_id_type slave_proc_id   = master_slave->processor_id();
        
        const dof_id_type n_dofs_on_proc_slave  = dof_map->n_dofs_on_processor(slave_proc_id);
        
        if(comm.isRoot()) {
            std::cout << "sum(B): " << volumes[0] <<std::endl;
        }
        
        
        
        express::Array<express::SizeType>  ownershipRangesMaster(comm.size()+1);
        ownershipRangesMaster.allSet(0);
        
        
        express::Array<express::SizeType>  ownershipRangesSlave(comm.size()+1);
        ownershipRangesSlave.allSet(0);
        
        
        const int n_nodes_x_face = master_slave->elem(0)->build_side_ptr(0)->n_nodes();
        express::Array<express::SizeType>  ownershipRangesBTilde = local_fun_spaces_new->ownershipRangesFaceID();
        for(SizeType i = 0; i < ownershipRangesBTilde.size(); ++i) {
            ownershipRangesBTilde[i] *= n_nodes_x_face;
        }
        
        
        ownershipRangesMaster[comm.rank() + 1] += static_cast<unsigned int>(n_dofs_on_proc_master);
        ownershipRangesSlave[comm.rank() +  1] += static_cast<unsigned int>(n_dofs_on_proc_slave);
        
        comm.allReduce(&ownershipRangesMaster[0], ownershipRangesMaster.size(), express::MPISum());
        comm.allReduce(&ownershipRangesSlave[0],  ownershipRangesSlave.size(),  express::MPISum());
        
        std::partial_sum(ownershipRangesMaster.begin(), ownershipRangesMaster.end(),
                         ownershipRangesMaster.begin());
        
        std::partial_sum(ownershipRangesSlave.begin(), ownershipRangesSlave.end(),
                         ownershipRangesSlave.begin());
        
        
        const SizeType local_range_slave_range  = ownershipRangesSlave [comm.rank() + 1] - ownershipRangesSlave [comm.rank()];
        const SizeType local_range_master_range = ownershipRangesMaster[comm.rank() + 1] - ownershipRangesMaster[comm.rank()];
        const SizeType local_range_b_tilde  = ownershipRangesBTilde[comm.rank() + 1] - ownershipRangesBTilde[comm.rank()];
        
        
        
        if(comm.isRoot()) {
            std::cout << "ownershipRangesMaster = "<< ownershipRangesMaster << std::endl;
        }
        
        mat_buffer.finalizeStructure();
        p_buffer.finalizeStructure();
        q_buffer.finalizeStructure();
        rel_area_buff.finalizeStructure();
        gap_buff.finalizeStructure();
        normal_buff.finalizeStructure();
        
        SizeType sizes[6] = {
            mat_buffer.rows(), mat_buffer.columns(),
            p_buffer.rows(), p_buffer.columns(),
            q_buffer.rows(), q_buffer.columns()
        };
        
        comm.allReduce(sizes, 6, express::MPIMax());
        
        SizeType fIndCols = express::Math<SizeType>::Max(sizes[0], sizes[1], sizes[3], sizes[5]);
        
        std::cout << sizes[0] << " " << sizes[1] << " " <<  sizes[3] << " " <<  sizes[5] << std::endl;
        std::cout << "NFaceDOFS=" << fIndCols << std::endl;
        
        mat_buffer.setSize(fIndCols, fIndCols);
        p_buffer.setSize(sizes[2], fIndCols);
        q_buffer.setSize(sizes[4], fIndCols);
        rel_area_buff.setSize(fIndCols, 1);
        gap_buff.setSize(fIndCols, 1);
        normal_buff.setSize(fIndCols, dim);
        
        express::Redistribute< express::MapSparseMatrix<double> > redist(comm.getMPIComm());
        redist.apply(ownershipRangesBTilde, mat_buffer, express::AddAssign<double>());
        redist.apply(ownershipRangesBTilde, gap_buff, express::AddAssign<double>());
        redist.apply(ownershipRangesBTilde, normal_buff, express::AddAssign<double>());
        
        redist.apply(local_fun_spaces_new->ownershipRangesFaceID(), rel_area_buff, express::AddAssign<double>());
        
        redist.apply(ownershipRangesSlave, p_buffer, express::AddAssign<double>());
        redist.apply(ownershipRangesMaster, q_buffer, express::AddAssign<double>());
        
        
        
        SizeType  mMaxRowEntries = mat_buffer.maxEntriesXCol();
        comm.allReduce(&mMaxRowEntries, 1, express::MPIMax());
        
        
        std::cout<<"local_range_b_tilde"<<local_range_b_tilde<<std::endl;
        // const SizeType local_range_master_range2 = ownershipRangesMaster[comm.rank()+1] - ownershipRangesMaster[comm.rank()];
        
        
        //        DVectord relAreaVec = zeros(local_fun_spaces_new->ownershipRangesFaceID()[comm.rank()+ 1] -
        //                                    local_fun_spaces_new->ownershipRangesFaceID()[comm.rank()]);
        
        express::RootDescribe("petsc rel_area_buff assembly begin", comm, std::cout);
        
        express::Array<bool> removeRow(local_range_b_tilde);
        
        if(!removeRow.isNull()) {
            removeRow.allSet(false);
            
            {
                //            utopia::Write<utopia::DVectord> write(relAreaVec);
                for (auto it = rel_area_buff.iter(); it; ++it) {
                    //                relAreaVec.add(it.row(), *it);
                    if(*it < 1 - 1e-8) {
                        const SizeType faceId = it.row();
                        for(int k = 0; k < n_nodes_x_face; ++k) {
                            const SizeType nodeId = faceId * n_nodes_x_face + k;
                            const SizeType index = nodeId - ownershipRangesBTilde[comm.rank()];
                            assert(index < removeRow.size());
                            removeRow[index] = true;
                        }
                    }
                }
            }
        }
        
        express::RootDescribe("petsc mat_buffer assembly begin", comm, std::cout);
        
        DVectord is_contact_node_tilde = local_zeros(local_range_b_tilde);
        DVectord gap_tilde = local_zeros(local_range_b_tilde);
        {
            utopia::Write<utopia::DVectord> write(gap_tilde);
            for (auto it = gap_buff.iter(); it; ++it) {
                
                const SizeType index = it.row() - ownershipRangesBTilde[comm.rank()];
                assert(index < removeRow.size());
                
                //std::cout<< comm.rank() << " row  "<< it.row() << " index " << index << " " << *it << std::endl;
                
                if(!removeRow[index]) {
                    gap_tilde.set(it.row(), *it);
                    is_contact_node_tilde.set(it.row(), 1.0);
                }
            }
        }
        
        
        DSMatrixd normal_tilde = utopia::local_sparse(local_range_b_tilde, dim, dim);
        {
            utopia::Write<utopia::DSMatrixd> write(normal_tilde);
            for (auto it = normal_buff.iter(); it; ++it) {
                
                const SizeType index = it.row() - ownershipRangesBTilde[comm.rank()];
                assert(index < removeRow.size());
                
                //std::cout<< comm.rank() << " row  "<< it.row() << " index " << index << " " << *it << std::endl;
                
                if(!removeRow[index]) {
                    normal_tilde.set(it.row(), it.col(), *it);
                }
            }
        }
        
        
        
        DSMatrixd B_tilde = utopia::local_sparse(local_range_b_tilde, local_range_b_tilde, mMaxRowEntries);
        {
            utopia::Write<utopia::DSMatrixd> write(B_tilde);
            for (auto it = mat_buffer.iter(); it; ++it) {
                
                const SizeType index = it.row() - ownershipRangesBTilde[comm.rank()];
                assert(index < removeRow.size());
                
                //std::cout<< comm.rank() << " row  "<< it.row() << " index " << index << " " << *it << std::endl;
                
                if(!removeRow[index]) {
                    B_tilde.set(it.row(), it.col(), *it);
                }
            }
        }
        
        SizeType  mMaxRowEntries_q = q_buffer.maxEntriesXCol();
        SizeType  mMaxRowEntries_p = p_buffer.maxEntriesXCol();
        
        comm.allReduce(&mMaxRowEntries_p, 1, express::MPIMax());
        comm.allReduce(&mMaxRowEntries_q, 1, express::MPIMax());
        
        comm.barrier();
        express::RootDescribe("petsc p_buffer assembly begin", comm, std::cout);
        
        
        DSMatrixd P = utopia::local_sparse(local_range_slave_range, local_range_b_tilde, mMaxRowEntries_p);
        {
            utopia::Write<utopia::DSMatrixd> write(P);
            for (auto it = p_buffer.iter(); it; ++it) {
                P.set(it.row(), it.col(), *it);
            }
        }
        
        comm.barrier();
        express::RootDescribe("petsc q_buffer assembly begin", comm, std::cout);
        
        std::cout << local_range_master_range << " == " << 81 << std::endl;
        DSMatrixd Q = utopia::local_sparse(local_range_master_range, local_range_b_tilde, mMaxRowEntries_q);
        
       // disp(size(Q));
        {
            utopia::Write<utopia::DSMatrixd> write(Q);
            for (auto it = q_buffer.iter(); it; ++it) {
                std::cout << it.row() << ", " << it.col() << std::endl;
                Q.set(it.row(), it.col(), *it);
                
            }
        }
        
        
        DSMatrixd Q_t = transpose(Q);
        
        
        DSMatrixd B_x = P * B_tilde * Q_t;
        
        
        normals = P * normal_tilde;
        
        DVectord gap_x = P * gap_tilde;
        
        DVectord normals_vec = local_zeros(local_range_slave_range);
        
        is_contact_node = P * is_contact_node_tilde;
        
        typedef Intersector::Scalar Scalar;
        
        
        {
            Write<DVectord> w(normals_vec);
            
            each_read(normals, [&](const SizeType i, const SizeType j, const double value){
                normals_vec.set(i + j, value);
            });
        }
        
        
        bool has_contact = false;
        
        // utopia::Range r = utopia::range(gap);
        each_read(is_contact_node, [&](const SizeType i , const double value){
            if (value > 0)
            {
                is_contact_node.set(i, 1);
                has_contact = true;
                std::cout << "expetected_contact_node: " << i << std::endl;
            }
        });
        
        
        
        orthogonal_trafos = local_sparse(local_range_slave_range , local_range_slave_range , dim);
        
        {
            std::vector<Scalar> normal(dim, 0);
            std::vector<Scalar> H(dim * dim, 0);
            
            Read<DVectord>  r_n(normals_vec);
            Write<DSMatrixd> w_o(orthogonal_trafos);
            Read<DVectord> r_icn(is_contact_node);


            // disp(is_contact_node);

            // MPI_Barrier(PETSC_COMM_WORLD);

            // std::cout << "normals_vec: " << std::endl;
            // disp(local_size(normals_vec));

            // std::cout << "is_contact_node: " << std::endl;
            // disp(local_size(is_contact_node));
            // std::cout << std::flush;

            // MPI_Barrier(PETSC_COMM_WORLD);

            
            bool check_has_contact = false;
            
            utopia::Range r = utopia::range(normals_vec);
            for(uint i = r.begin(); i < r.end(); i += dim) {
                bool use_identity = true;
                
                
                // bool is_cn_i = is_contact_node.get(i/dim) > 0;
                bool is_cn_i = is_contact_node.get(i) > 0;
                
                if(is_cn_i) {
                    std::cout << "actual_contact_node: " << i << std::endl;
                    
                    check_has_contact = true;
                    
                    for(uint d = 0; d < dim; ++d) {
                        normal[d] = normals_vec.get(i + d);
                    }
                    
                    // if(std::abs(normal[0] - 1.) > 1e-8) {
                        use_identity = false;
                        
                        normalize(normal);
                        //-e1 basis vector
                        normal[0] -= 1;
                        normalize(normal);

                        
                        if(dim == 2) {
                            isector.householder_reflection_2(&normal[0], &H[0]);
                        } else {
                            isector.householder_reflection_3(&normal[0], &H[0]);
                        }
                        
                        
                        for(uint di = 0; di < dim; ++di) {
                            for(uint dj = 0; dj < dim; ++dj) {
                                orthogonal_trafos.set(i + di, i + dj, H[di * dim + dj]);
                            }
                        }
                    // }
                }
                
                if(use_identity)
                {
                    for(uint di = 0; di < dim; ++di) {
                        orthogonal_trafos.set(i + di, i + di, 1.);
                    }
                }
            }
            
            if(!check_has_contact == has_contact) {
                std::cerr << "inconsistent contact determination" << std::endl;
            }
        }
        
        

        auto s_gap = local_size(gap_x);
        gap = local_zeros(s_gap);

        static const double LARGE_VALUE = 10000;
        {
            Write<DVectord> w_g(gap);
            Read<DVectord> r_icn(is_contact_node);

            each_read(gap_x, [&](const SizeType i, const double value) {
                const SizeType offset = i;

                if(is_contact_node.get(i) > 0) {
                    gap.set(offset, value);
                } else {
                    gap.set(offset, LARGE_VALUE);
                }
            });
        }
        
        auto s_B_x = local_size(B_x);
        B = local_sparse(s_B_x.get(0), s_B_x.get(1), mMaxRowEntries_q * dim);

        {
            Write<DSMatrixd> w_B(B);
            each_read(B_x, [&](const SizeType i, const SizeType j, const double value) {
                for(SizeType d = 0; d < dim; ++d) {
                    B.set(i + d, j + d, value);    
                } 
            }); 
        }

        // disp("B_tilde:");
        // disp(B.size());
        
        // disp("Q:");
        // disp(Q_t.size());
        
        // disp("P:");
        // disp(P.size());
        
        // disp("B:");
        // disp(B.size());
        
        
        // write("B_tilde.m", B_tilde);
        // write("Q.m", Q);
        // write("P.m", P);

        // write("B.m", B);
        // write("O.m", orthogonal_trafos);
        // write("g.m", gap);
        // write("c.m", is_contact_node);
        
        express::RootDescribe("petsc assembly end", comm, std::cout);
        
        return true;
    }
    
    
    
  inline bool MooseSurfaceAssemble(express::Communicator &comm,
                                  const std::shared_ptr<MeshBase> &mesh,
                                  const std::shared_ptr<DofMap> &dof_map,
                                  const std::shared_ptr<const unsigned int> &_var_num,
                                  DSMatrixd &B,DSMatrixd &orthogonal_trafos,DVectord &gap, DSMatrixd &normals,
                                  DVectord &is_contact_node,
                                  const libMesh::Real search_radius, const int tag_1, const int tag_2)
    {
        
        cutk::Settings settings;
        
//        auto predicate = std::make_shared<cutlibpp::MasterAndSlave>();
//        
//        predicate->add(tag_1,tag_2);
        
        
//        express::Communicator comm = libmesh_comm_.get();
        
        if(mesh->mesh_dimension() == 2) {
            //std::cout<<"Assemble_matrix::I am in assemble"<<std::endl;
            return utopia::SurfaceAssemble<2>(comm, mesh, dof_map, _var_num, B,  orthogonal_trafos, gap, normals, is_contact_node, settings, search_radius, tag_1, tag_2);
        }
        
        
        if(mesh->mesh_dimension() == 3) {
            return utopia::SurfaceAssemble<3>(comm, mesh, dof_map, _var_num, B,  orthogonal_trafos, gap, normals,  is_contact_node, settings, search_radius, tag_1, tag_2);
        }
        
        assert(false && "Dimension not supported!");
        return false;
    }

//
//    template<int Dimensions, class Fun>
//    static bool Assemble(express::Communicator &comm,
//                         const std::shared_ptr<MeshBase> &master,
//                         const std::shared_ptr<MeshBase> &slave,
//                         const std::shared_ptr<DofMap> &dof_master,
//                         const std::shared_ptr<DofMap> &dof_slave,
//                         const std::shared_ptr<const unsigned int> &_from_var_num,
//                         const std::shared_ptr<const unsigned int> &_to_var_num,
//                         Fun process_fun,
//                         const cutk::Settings &settings, bool use_biorth_)
//    {
//        
//        
//        using namespace cutlibpp;
//        using namespace express;
//        using namespace cutk;
//        
//        typedef LibMeshTree<Dimensions> NTreeT;
//        typedef typename NTreeT::DataContainer DataContainer;
//        typedef typename NTreeT::DataType Adapter;
//        
//        long maxNElements = 100;
//        long maxDepth = 5;
//        
//        
//        if (!settings.get("max_depth").isNull()) {
//            maxDepth = settings.get("max_depth").toInt();
//            std::cout<<"max_depth  = "<< maxDepth  <<std::endl;
//        }
//        
//        const auto &master_mesh = master;
//        const auto &slave_mesh  = slave;
//        const int n_elements_master = master_mesh->n_elem();
//        const int n_elements_slave  = slave_mesh->n_elem();
//        const int n_elements 		= n_elements_master + n_elements_slave;
//        
//        
//        const Parallel::Communicator &libmesh_comm_master = master_mesh->comm();
//        const Parallel::Communicator &libmesh_comm_slave = slave_mesh->comm();
//        
//        
//        auto predicate = make_shared<MasterAndSlave>();
//        predicate->add(0, 1);
//        
//        EXPRESS_EVENT_BEGIN("create_adapters");
//        ////////////////////////////////////////////////////////////////////////////////////////////////////
//        cutk::shared_ptr<NTreeT> tree = NTreeT::New(predicate, maxNElements, maxDepth);
//        std::cout<<"n_elements"<<n_elements<<std::endl;
//        tree->reserve(n_elements);
//        
//        
//        std::shared_ptr<Spaces> local_spaces = make_shared<Spaces>(master, slave, dof_master, dof_slave, _from_var_num, _to_var_num);
//        
//        int offset = 0;
//        int space_num = 0;
//        
//        for(auto s : local_spaces->spaces()) {
//            if(s) {
//                
//                bool first = true;
//                for (auto it = s->active_local_elements_begin(); it != s->active_local_elements_end(); ++it) {
//                    auto elem=*it;
//                    Adapter a(*s, elem->id(), offset+elem->id(), space_num);
//                    assert(!local_spaces->dof_map(space_num)[elem->id()].empty());
//                    a.set_dof_map(&local_spaces->dof_map(space_num)[elem->id()].global);
//                    tree->insert(a);
//                }
//                
//                offset += s->n_elem(); //s->mesh().n_active_local_elem();//(*s->mesh().active_local_elements_end())->id();
//                
//            }
//            
//            
//            ++space_num;
//            
//            
//        }
//        
//        //std::cout<<" --------------------------- tree->memory().nData()=" << tree->memory().nData()<<std::endl; ;
//        
//        tree->getRoot()->getBound().staticBound().enlarge(1e-8);
//        
//        ////////////////////////////////////////////////////////////////////////////////////////////////////
//        EXPRESS_EVENT_END("create_adapters");
//        
//        // std::cout<<"-----------------------------------ADAPTERS-------------------------------------------------"<<std::endl;
//        //Just to have an indexed-storage
//        std::map<long, cutk::shared_ptr<Spaces> > spaces;
//        std::map<long, std::vector<cutk::shared_ptr<Spaces> > > migrated_spaces;
//        
//        
//        auto read = [&spaces, &migrated_spaces, comm, &libmesh_comm_master, &libmesh_comm_slave ]
//        (
//         const long ownerrank,
//         const long senderrank,
//         bool is_forwarding, DataContainer &data,
//         InputStream &in
//         ) {
//            
//            
//            //   std::cout<<"------------------------------------AUTO-READ-IN---------------------------------------------"<<std::endl;
//            
//            
//            CHECK_STREAM_READ_BEGIN("vol_proj", in);
//            
//            cutk::shared_ptr<Spaces> proc_space = cutk::make_shared<Spaces>(comm);
//            
//            read_spaces(in, *proc_space, libmesh_comm_master, libmesh_comm_slave);
//            
//            if (!is_forwarding) {
//                assert(!spaces[ownerrank]);
//                spaces[ownerrank] = proc_space;
//            } else {
//                migrated_spaces[ownerrank].push_back(proc_space);
//            }
//            
//            data.reserve(data.size() + 3000);
//            //std::cout<<proc_space->dof_map(0)[0].global<<std::endl;
//            
//            int space_num = 0;
//            long offset = 0;
//            for(auto s : proc_space->spaces()) {
//                if(s) {
//                    for (int i=0; i<s->n_elem(); i++) {
//                        data.push_back(Adapter(*s, i, offset + i, space_num) );
//                        assert(!proc_space->dof_map(space_num)[i].empty());
//                        data.back().set_dof_map(&proc_space->dof_map(space_num)[i].global);
//                    }
//                    
//                    offset += s->n_elem();
//                    
//                }
//                
//                ++space_num;
//                
//            }
//            
//            
//            //    std::cout<<"------------------------------------AUTO-READ-OUT---------------------------------------------"<<std::endl;
//            
//            CHECK_STREAM_READ_END("vol_proj", in);
//            
//            
//            
//        };
//
//        
//        auto write = [&local_spaces, &spaces, &comm]
//        (
//         const long ownerrank, const long recvrank,
//         const std::vector<long>::const_iterator &begin,
//         const std::vector<long>::const_iterator &end,
//         const DataContainer &data,
//         OutputStream &out) {
//            
//            CHECK_STREAM_WRITE_BEGIN("vol_proj", out);
//            
//            //       std::cout<<"------------------------------------AUTO-WRITE-IN-------------------------------------------"<<std::endl;
//            
//            if (ownerrank == comm.rank()) {
//                
//                write_element_selection(begin, end, *local_spaces, out);
//                
//                
//            } else {
//                
//                auto it = spaces.find(ownerrank);
//                assert(it != spaces.end());
//                cutk::shared_ptr<Spaces> spaceptr = it->second;
//                assert(std::distance(begin, end) > 0);
//                write_element_selection(begin, end, *spaceptr, out);
//                
//            }
//            
//            //      std::cout<<"------------------------------------AUTO-WRITE-OUT-------------------------------------------"<<std::endl;
//            //      comm.barrier();
//            
//            CHECK_STREAM_WRITE_END("vol_proj", out);
//            
//        };
//
//        
//        long n_false_positives = 0, n_intersections = 0;
//        
//        //       std::cout<<"------------------------------------SEARCH-COMPUTE-------------------------------------------"<<std::endl;
//        
//        auto fun = [&n_false_positives, &n_intersections, &process_fun](
//                                                                        
//                                                                        Adapter &master, Adapter &slave) -> bool {
//            
//            bool ok = process_fun(master, slave);
//            
//            if(ok) {
//                n_intersections++;
//                
//                return true;
//            } else {
//                
//                n_false_positives++;
//                return false;
//            }
//            return true;
//            
//        };
//
//        
//
//        cutk::Settings custom_settings = settings;
//        custom_settings.set("disable_redistribution", cutk::Boolean(true));
//        custom_settings.set("verbosity_level", cutk::Integer(2));
//        
//        
//        cutlibpp::search_and_compute(comm, tree, predicate, read, write, fun, custom_settings);
//        
//        
//        
//        long n_total_candidates = n_intersections + n_false_positives;
//        
//        long n_collection[3] = {n_intersections, n_total_candidates, n_false_positives};
//        
//        
//        comm.allReduce(n_collection, 3, express::MPISum());
//        
//        
//        if (comm.isRoot()) {
//            std::cout << "n_intersections: " << n_collection[0]
//            << ", n_total_candidates: " 	 << n_collection[1]
//            << ", n_false_positives: " 	     << n_collection[2] << std::endl;
//        }
//        
//        return true;
//    }
//
//
//    template<int Dimensions>
//    bool Assemble(
//                  express::Communicator &comm,
//                  const std::shared_ptr<MeshBase> &master,
//                  const std::shared_ptr<MeshBase> &slave,
//                  const std::shared_ptr<DofMap> &dof_master,
//                  const std::shared_ptr<DofMap> &dof_slave,
//                  const std::shared_ptr<const unsigned int> &_from_var_num,
//                  const std::shared_ptr<const unsigned int> &_to_var_num,
//                  DSMatrixd &B,
//                  const cutk::Settings &settings,bool  use_biorth_)
//    {
//        std::shared_ptr<Spaces> local_fun_spaces = cutk::make_shared<Spaces>(master, slave, dof_master, dof_slave,_from_var_num,_to_var_num);
//        
//        libMesh::DenseMatrix<libMesh::Real> src_pts;
//        libMesh::DenseMatrix<libMesh::Real> dest_pts;
//        libMesh::DenseMatrix<libMesh::Real> intersection2;
//        Polyhedron src_poly, dest_poly;
//        Polyhedron  intersection3,temp_poly;
//        Intersector isector;
//        
//        std::shared_ptr<MeshBase> master_space = master;
//        std::shared_ptr<MeshBase> slave_space  = slave;
//        
//        
//       // std::vector<libMesh::dof_id_type> master_dofs, slave_dofs;
//        libMesh::DenseMatrix<libMesh::Real> elemmat;
//        libMesh::DenseMatrix<libMesh::Real> cumulative_elemmat;
//        
//        std::shared_ptr<Transform> src_trans;
//        std::shared_ptr<Transform> dest_trans;
//        
//        
//        int skip_zeros = 1;
//        
//        
//        libMesh::Real total_intersection_volume = 0.0;
//        libMesh::Real local_element_matrices_sum = 0.0;
//        
//        
//        
//        express::MapSparseMatrix<double> mat_buffer(dof_slave->n_dofs(), dof_master->n_dofs());
//        
////        std::cout<<"dof_slave->n_dofs()"<<dof_slave->n_dofs()<<std::endl;
////        std::cout<<"dof_master->n_dofs()"<<dof_master->n_dofs()<<std::endl;
//        
//        bool intersected = false;
//        
//        double element_setup_time = 0.0;
//        double intersection_time = 0.0;
//        double assembly_time     = 0.0;
//        
//        utopia::Chrono c;
//        
//        auto fun = [&](const ElementAdapter<Dimensions> &master,
//                       const ElementAdapter<Dimensions> &slave) -> bool {
//            
//            c.start();
//            
//            long n_intersections = 0;
//            
//            bool pair_intersected = false;
//            
//            const auto &src  = master.space();
//            
//            const auto &dest = slave.space();
//            
//            const auto &src_mesh  = src;
//            
//            const auto &dest_mesh = dest;
//            
//            const int src_index  = master.element();
//            
//            const int dest_index = slave.element();
//            
//            auto &src_el  = *src_mesh.elem(src_index);
//            
//            auto &dest_el = *dest_mesh.elem(dest_index);
//            
//            const int dim = src_mesh.mesh_dimension();
//            
//            
//            std::unique_ptr<libMesh::FEBase> master_fe, slave_fe;
//            
//            master_fe = libMesh::FEBase::build(src_mesh.mesh_dimension(),  dof_master->variable_type(0));
//            slave_fe  = libMesh::FEBase::build(dest_mesh.mesh_dimension(), dof_slave->variable_type(0));
//            
//            QMortar composite_ir(dim);
//            QMortar src_ir(dim);
//            QMortar dest_ir(dim);
//            
//            
//            const int order = order_for_l2_integral(dim, src_el, dof_master->variable(0).type().order , dest_el,dof_slave->variable(0).type().order);
//            
//            c.stop();
//            element_setup_time += c.get_seconds();
//            c.start();
//            
//            if(dim == 2)  {
//                make_polygon(src_el,   src_pts);
//                make_polygon(dest_el, dest_pts);
//                
//                if(intersect_2D(src_pts, dest_pts, intersection2)) {
//                    total_intersection_volume += fabs(isector.polygon_area_2(intersection2.m(), &intersection2.get_values()[0]));
//                    
//                    const libMesh::Real weight=isector.polygon_area_2(dest_pts.m(), &dest_pts.get_values()[0]);
//                    
//                    make_composite_quadrature_2D(intersection2, weight, order, composite_ir);
//                    pair_intersected = true;
//                    
//                    src_trans  = std::make_shared<Transform2>(src_el);
//                    dest_trans = std::make_shared<Transform2>(dest_el);
//                    pair_intersected = true;
//                }
//            }
//            else if(dim == 3) {
//                make_polyhedron(src_el,  src_poly);
//                make_polyhedron(dest_el, dest_poly);
//                
//                
//                if(intersect_3D(src_poly, dest_poly, intersection3)) {
//                    
//                    total_intersection_volume += isector.p_mesh_volume_3(intersection3);
//                    
//                    const libMesh::Real weight = isector.p_mesh_volume_3(dest_poly);
//                    
//                    make_composite_quadrature_3D(intersection3, weight, order, composite_ir);
//                    src_trans  = std::make_shared<Transform3>(src_el);
//                    dest_trans = std::make_shared<Transform3>(dest_el);
//                    pair_intersected = true;
//                }
//                
//            } else {
//                assert(false);
//                return false;
//            }
//            
//            c.stop();
//            intersection_time += c.get_seconds();
//            c.start();
//            
//            const auto &master_dofs = master.dof_map();
//            const auto &slave_dofs  = slave.dof_map();
//            
//            if(pair_intersected) {
//                
//                
//                transform_to_reference(*src_trans,  src_el.type(),  composite_ir,  src_ir);
//                transform_to_reference(*dest_trans, dest_el.type(), composite_ir,  dest_ir);
//                
//                //            src.dof_map().dof_indices(&src_el,  master_dofs);
//                //            dest.dof_map().dof_indices(&dest_el, slave_dofs);
//                
//                
//       
//                
//                assert(!master_dofs.empty());
//                assert(!slave_dofs.empty());
//                //composite_ir.print_info();
//                
//                
//                master_fe->attach_quadrature_rule(&src_ir);
//                master_fe->reinit(&src_el);
//                
//                slave_fe->attach_quadrature_rule(&dest_ir);
//                slave_fe->reinit(&dest_el);
//                
//                elemmat.zero();
//                
//                
//                
//                if(use_biorth_) {
//                    mortar_assemble_biorth(*master_fe, *slave_fe, dest_el.type(), elemmat);
//                    
//                } else {
//                    mortar_assemble(*master_fe, *slave_fe, elemmat);
//                }
//                
//                // std::cout << "-----------------------------------------\n";
//                // std::cout << src_index << ", " << dest_index << "\n";
//                // elemmat.print(std::cout);
//                // for(auto i : slave_dofs) {
//                // 	std::cout << i << " ";
//                // }
//                // std::cout << "\n";
//                
//                // for(auto i : master_dofs) {
//                // 	std::cout << i << " ";
//                // }
//                // std::cout << "\n";
//                // std::cout << "-----------------------------------------\n";
//                
//                auto partial_sum = std::accumulate(elemmat.get_values().begin(), elemmat.get_values().end(), libMesh::Real(0.0));
//                // std::cout << src_index << ", " << dest_index << ": " << partial_sum << std::endl;
//                // dest_ir.print_info();
//                
//                local_element_matrices_sum += partial_sum;
//                
//                intersected = true;
//                
//                ++n_intersections;
//                
//                
//                if(slave_dofs.size() != elemmat.m()) {
//                    std::cout << slave_dofs.size() << " != " <<  elemmat.m() << std::endl;
//                }
//                
//                assert(slave_dofs.size() == elemmat.m());
//                assert(master_dofs.size() == elemmat.n());
//                
//               // std::cout<<"slave_dofs.size()"<<slave_dofs.size()<<std::endl;
//               // std::cout<<"master_dofs.size()"<<master_dofs.size()<<std::endl;
//                
//                for(int i = 0; i < slave_dofs.size(); ++i) {
//                    
//                    const long dof_I = slave_dofs[i];
//                    
//                    for(int j = 0; j < master_dofs.size(); ++j) {
//                        
//                        const long dof_J = master_dofs[j];
//                        
//                        mat_buffer.add(dof_I, dof_J, elemmat(i, j));
//                    }
//                }
//                
//                return true;
//                
//            } else {
//                
//                return false;
//            }
//            
//        };
//
//        
//        
//        // comm.barrier();
//        // utopia::Chrono c2;
//        // c2.start();
//        
//        
//        if(!Assemble<Dimensions>(comm, master, slave, dof_master, dof_slave, _from_var_num, _to_var_num, fun, settings, use_biorth_)) {
//            std::cout << "n_intersections: false2" <<std::endl;
//            return false;
//        }
//        
//        // c2.stop();
//        // std::cout << "Local mortars" << std::endl;
//        // c2.describe(std::cout);
//        
//        
//        // comm.barrier();
//        // std::stringstream ss;
//        // ss << "Setup_time: " << element_setup_time << "\n";
//        // ss << "intersection_time: " << intersection_time << "\n";
//        // ss << "assembly_time: " << assembly_time << std::endl;
//        
//        //express::SynchDescribe(ss.str(), comm, std::cout);
//        // comm.barrier();
//        // c2.start();
//        
//        std::cout << "n_intersections: " <<std::endl;
//        
//        double volumes[2] = { local_element_matrices_sum,  total_intersection_volume };
//        
//        comm.allReduce(volumes, 2, express::MPISum());
//        
//        const processor_id_type master_proc_id  = master->processor_id();
//        
//        const dof_id_type n_dofs_on_proc_master = dof_master->n_local_dofs();
//        
//        const processor_id_type slave_proc_id   = slave->processor_id();
//        
//        const dof_id_type n_dofs_on_proc_slave  =dof_slave->n_local_dofs();
//        
//        const int n_dofs_on_proc_print  = dof_slave->n_local_dofs();
//        
//        if(comm.isRoot()) {
//            std::cout << "sum(B): " << volumes[0] << ", vol(I): " << volumes[1] << std::endl;
//        }
//        
//        express::Array<express::SizeType>  ownershipRangesMaster(comm.size()+1);
//        ownershipRangesMaster.allSet(0);
//        
//        express::Array<express::SizeType>  ownershipRangesSlave(comm.size()+1);
//        ownershipRangesSlave.allSet(0);
//        
//        
//        ownershipRangesMaster[comm.rank()+1] += static_cast<unsigned int>(n_dofs_on_proc_master);
//        
//        ownershipRangesSlave[comm.rank()+1]  += static_cast<unsigned int>(n_dofs_on_proc_slave);
//        
//        
//        
//        comm.allReduce(&ownershipRangesMaster[0], ownershipRangesMaster.size(), express::MPISum());
//        
//        comm.allReduce(&ownershipRangesSlave[0],  ownershipRangesSlave.size(),  express::MPISum());
//        
//        
//        std::partial_sum(ownershipRangesMaster.begin(), ownershipRangesMaster.end(), ownershipRangesMaster.begin());
//        std::partial_sum(ownershipRangesSlave.begin(), ownershipRangesSlave.end(), ownershipRangesSlave.begin());
//        
//        
//        if(comm.isRoot()) {
//            std::cout <<ownershipRangesMaster << std::endl;
//            std::cout<<"prova"<<n_dofs_on_proc_print<<std::endl;
//            
//        }
//        
//        
//        
//        
//        express::Redistribute< express::MapSparseMatrix<double> > redist(comm.getMPIComm());
//        
//        redist.apply(ownershipRangesSlave, mat_buffer, express::AddAssign<double>());
//        
//        assert(ownershipRangesSlave.empty() == ownershipRangesMaster.empty() || ownershipRangesMaster.empty());
//        
//        express::RootDescribe("petsc assembly begin", comm, std::cout);
//        
//        SizeType  mMaxRowEntries = mat_buffer.maxEntriesXCol();
//        
//        comm.allReduce(&mMaxRowEntries, 1, express::MPIMax());
//        
//        const SizeType local_range_slave_range  = ownershipRangesSlave [comm.rank()+1] - ownershipRangesSlave [comm.rank()];
//        const SizeType local_range_master_range = ownershipRangesMaster[comm.rank()+1] - ownershipRangesMaster[comm.rank()];
//        
//        B = utopia::local_sparse(local_range_slave_range, local_range_master_range, mMaxRowEntries);
//        
//        {
//            utopia::Write<utopia::DSMatrixd> write(B);
//            for (auto it = mat_buffer.iter(); it; ++it) {
//                B.set(it.row(), it.col(), *it);
//                
//            }
//        }
//        
//        express::RootDescribe("petsc assembly end", comm, std::cout);
//        
//        // c2.stop();
//        // std::cout << "Global stuff\n";
//        // c2.describe(std::cout);
//        return true;
//        
//    }
//
//
//
//inline bool AssembleMOOSE(express::Communicator &comm,
//                       const std::shared_ptr<MeshBase> &master,
//                       const std::shared_ptr<MeshBase> &slave,
//                       const std::shared_ptr<DofMap> &dof_master,
//                       const std::shared_ptr<DofMap> &dof_slave,
//                       const std::shared_ptr<const unsigned int> & _from_var_num,
//                       const std::shared_ptr<const unsigned int> & _to_var_num,
//                       bool  use_biorth_,
//                       DSMatrixd &B)
//    {
//        cutk::Settings settings;
//        
//        if(master->mesh_dimension() == 2) {
//            std::cout<<"Assemble_matrix::I am in assemble"<<std::endl;
//            return utopia::Assemble<2>(comm, master, slave, dof_master, dof_slave, _from_var_num,  _to_var_num, B, settings,use_biorth_);
//        }
//        
//        if(master->mesh_dimension() == 3) {
//            std::cout<<"Assemble_matrix::I am in assemble"<<std::endl;
//            return utopia::Assemble<3>(comm, master, slave, dof_master, dof_slave, _from_var_num,  _to_var_num, B, settings,use_biorth_);
//        }
//        
//        assert(false && "Dimension not supported!");
//        return false;
//    }
////
////
//////    bool AssembleMOOSE(express::Communicator &comm,
//////                       const std::shared_ptr<MeshBase> &mesh_master,
//////                       const std::shared_ptr<MeshBase> &mesh_slave,
//////                       libMesh::Order master_order,
//////                       libMesh::Order slave_order,
//////                       DSMatrixd &B)
//////    {
//////        cutk::Settings settings;
//////        
//////    
//////        
//////        LibMeshFEContext<LinearImplicitSystem> master_context(mesh_master);
//////        auto master_space = fe_space(LAGRANGE, master_order, master_context);
//////        master_context.equation_systems.init();
//////        
//////        LibMeshFEContext<LinearImplicitSystem> slave_context(mesh_slave);
//////        auto slave_space = fe_space(LAGRANGE, slave_order, slave_context);
//////        slave_context.equation_systems.init();
//////
//////        
//////        if(mesh_master->mesh_dimension() == 2) {
//////            std::cout<<"Assemble_matrix::I am in assemble"<<std::endl;
//////            return utopia::Assemble<2>(comm, make_ref(master_space), make_ref(slave_space), B, settings);
//////        }
//////        
//////        if(mesh_master->mesh_dimension() == 3) {
//////            std::cout<<"Assemble_matrix::I am in assemble"<<std::endl;
//////            return utopia::Assemble<3>(comm, make_ref(master_space), make_ref(slave_space), B, settings);
//////        }
//////        
//////        assert(false && "Dimension not supported!");
//////        return false;
//////    }
//////    
////
////
////
////   
////
}

#endif //LIBMESH_CUTLIBPP_ADAPTERS_HPP

    

