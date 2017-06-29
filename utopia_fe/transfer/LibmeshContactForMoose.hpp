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

#include "utopia_Socket.hpp"

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
    
    
    //Hacky but it is the only way with libmesh that we know of without having to clone the boundary info!!!!
    inline static void nodes_are_boundary_hack(
                                               const libMesh::DenseMatrix<libMesh::Real> &mat,
                                               std::vector<bool> &rows,
                                               std::vector<bool> &cols)
    {
        
        rows.resize(mat.m());
        cols.resize(mat.n());
        
        std::fill(rows.begin(), rows.end(), 0);
        std::fill(cols.begin(), cols.end(), 0);
        
        std::vector<double> sum_rows(mat.m(), 0.);
        std::vector<double> sum_cols(mat.n(), 0.);
        double sum_all = 0.0;
        double vol_check = 0.0;
        
#ifndef NDEBUG
        std::vector<double> sum_rows_with_sign(mat.m(), 0.);
        
#endif //NDEBUG
        
        for(int i = 0; i < sum_rows.size(); ++i) {
            for(int j = 0; j < sum_cols.size(); ++j) {
                const double abs_ij = std::abs(mat(i, j));
                sum_rows[i] += abs_ij;
                sum_cols[j] += abs_ij;
                sum_all += abs_ij;
                
                vol_check += mat(i, j);
                
#ifndef NDEBUG
                sum_rows_with_sign[i] += mat(i, j);
#endif //NDEBUG
            }
        }
        
        for(int i = 0; i < sum_rows.size(); ++i) {
            if(std::abs(sum_rows[i]/sum_all) > 1e-8) {
                rows[i] = true;
            }
        }
        
        for(int i = 0; i < sum_cols.size(); ++i) {
            if(std::abs(sum_cols[i]/sum_all) > 1e-8) {
                cols[i] = true;
            }
        }
        
        
        assert(vol_check > 0);
#ifndef NDEBUG
        for(int i = 0; i < sum_rows_with_sign.size(); ++i) {
            assert(sum_rows_with_sign[i]/sum_all >= -1e-8);
        }
#endif //NDEBUG
    }

        inline static bool check_node_is_boundary(const ElemType &type,
                                                  const std::vector<bool> &is_boundary)
        {

            int n_bound = 0;
            for(std::size_t i = 0; i < is_boundary.size(); ++i) {
                n_bound += is_boundary[i];
            }
            
            switch(type) {
                case TRI3:
                case QUAD4:
                case QUADSHELL4:
                {
                    assert(n_bound == 2);
                    return n_bound == 2;   
                }
                    
                case TET4:
                {
                    assert(n_bound == 3);
                    return n_bound == 3;   
                }
                    
                default:
                {
                    std::cerr << "Missing implementation for ElemType " << type << std::endl;
                    assert(false && "implement me!");
                    break;
                }
            }
        }
    
    
    
    
    inline static void assemble_trace_biorth_weights_from_space(const ElemType &type,
                                                                const std::vector<bool> &is_boundary,
                                                                libMesh::DenseMatrix<libMesh::Real> &weights)
    {
        switch(type) {
            case TRI3:
            {
                weights.resize(3, 3);
                weights.zero();
                
#ifndef NDEBUG
                int n_bound = 0;
                for(std::size_t i = 0; i < is_boundary.size(); ++i) {
                    n_bound += is_boundary[i];
                }
                
                assert(n_bound == 2);
#endif //NDEBUG
                
                for(std::size_t i = 0; i < is_boundary.size(); ++i) {
                    if(!is_boundary[i]) { continue; }
                    
                    weights(i, i) = 2;
                    
                    for(std::size_t j = 0; j < is_boundary.size(); ++j) {
                        if(is_boundary[j] && i != j) {
                            weights(i, j) = -1;
                        }
                    }
                }
                
                break;
            }
                
            case QUAD4:
            case QUADSHELL4:
            {
                
                weights.resize(4, 4);
                weights.zero();
#ifndef NDEBUG
                int n_bound = 0;
                for(std::size_t i = 0; i < is_boundary.size(); ++i) {
                    n_bound += is_boundary[i];
                }
                
                assert(n_bound == 2);
#endif //NDEBUG
                
                for(std::size_t i = 0; i < is_boundary.size(); ++i) {
                    if(!is_boundary[i]) { continue; }
                    
                    weights(i, i) = 2;
                    
                    for(std::size_t j = 0; j < is_boundary.size(); ++j) {
                        if(is_boundary[j] && i != j) {
                            weights(i, j) = -1;
                        }
                    }
                }
                
                break;
            }
                
            case TET4:
            {
                weights.resize(4, 4);
                weights.zero();
                
#ifndef NDEBUG
                int n_bound = 0;
                for(std::size_t i = 0; i < is_boundary.size(); ++i) {
                    n_bound += is_boundary[i];
                }
                
                assert(n_bound == 3);
#endif //NDEBUG
                
                for(std::size_t i = 0; i < is_boundary.size(); ++i) {
                    if(!is_boundary[i]) { continue; }
                    
                    weights(i, i) = 3;
                    
                    for(std::size_t j = 0; j < is_boundary.size(); ++j) {
                        if(is_boundary[j] && i != j) {
                            weights(i, j) = -1;
                        }
                    }
                }
                
                break;
            }
                
            default:
            {
                std::cerr << "Missing implementation for ElemType " << type << std::endl;
                assert(false && "implement me!");
                break;
            }
        }
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
        //for(int side_number = 0; side_number < e.n_sides(); ++side_number)
            //boundary_face_dof_map[side_number] gives the dofs associated to 
            //if(boundary_face_dof_map[side_number].empty()) { means is not a boundary face }
        //the boundary face
        //std::vector< std::vector<long> > boundary_face_dof_map;

        //get_side_dof_query(side, std::vector<bool> &is_boundary)
        
        
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
                   const std::shared_ptr<const unsigned int> &_var_num,  const std::vector< std::pair<int, int>> &tags)
        {
            
            utopiamesh_.reserve(1);
            
            utopiamesh_.push_back(master_slave);
            
            must_destroy_attached[0] = false;
            
            express::Communicator comm = master_slave->comm().get();
            
            const int n_elements = master_slave->n_elem();
            
            copy_global_dofs(*master_slave, original_dofmap, _var_num,
                             dof_maps_[0], var_type_[0], n_elements, var_number_[0],
                             subdomain_id_[0], side_set_id_[0], side_set_id_tag_[0],
                             face_set_id_global_[0],ownershipRangesFaceID_[0], tags); 
            
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
        
        
        
        inline const std::vector<ElementDofMap> &side_set_id_tag() const
        {
            return side_set_id_tag_[0];
        }
        
        inline express::Array<express::SizeType> &ownershipRangesFaceID()
        {
            
            return ownershipRangesFaceID_[0];
        }
        
        
        
        inline const express::Array<express::SizeType> &ownershipRangesFaceID() const
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
        
        
       inline static bool is_tagged_contact_boundary(const MeshBase &mesh, 
                                                const Elem *elem,
                                                const int side_elem,
                                                const std::vector< std::pair<int, int> > &tags)
       {
           for(auto t : tags) {
               if (mesh.get_boundary_info().has_boundary_id(elem, side_elem, t.first) || 
                   mesh.get_boundary_info().has_boundary_id(elem, side_elem, t.second)) {
                    return true;
                }
            }

            return false;
        }

        
        
        
        inline static void copy_global_dofs(MeshBase &space, const std::shared_ptr<libMesh::DofMap>  &original_dofmap,
                                            const std::shared_ptr<const unsigned int>  &var_num,
                                            std::vector<ElementDofMap> &dof_map, std::vector<ElementDofMap> &variable_type,
                                            const int n_elements, std::vector<ElementDofMap> &variable_number,
                                            std::vector<ElementDofMap> &subdomain_id, std::vector<ElementDofMap> &side_set_id,
                                            std::vector<ElementDofMap> &side_set_id_tag, std::vector<ElementDofMap> &face_set_id_global,
                                            express::Array<express::SizeType>  & ownershipRangesFaceID, const std::vector< std::pair<int, int> >  &tags)
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
                    if (check_side_id_one) {
                        // side_set_id[elem->id()].global.insert(side_set_id[elem->id()].global.end(),-1);
                        // side_set_id[elem->id()].global.push_back(-1);
                        check_side_id_one=false;
                        jj_side_id_one++;
                    }
                }
                
                if (elem->on_boundary()){
                    for (int side_elem=0; side_elem<elem->n_sides(); side_elem++){
                        {
                            if (is_tagged_contact_boundary(mesh, elem, side_elem, tags) &&
                                check_side_id_one_tag){
                                // side_set_id[elem->id()].global.insert(side_set_id[elem->id()].global.end()-1,mesh.get_boundary_info().boundary_id(elem,side_elem));
                                side_set_id[elem->id()].global.push_back(mesh.get_boundary_info().boundary_id(elem, side_elem));
                                check_side_id_one_tag = false;
                                jj_side_id_one_tag++;
                            }
                        }
                    }
                }
                
                
                
                for (int side_elem=0; side_elem<elem->n_sides(); side_elem++){
                    if (check_side_id_one_check){
                        check_side_id_one_check=false;
                        jj_side_id_one_check++;
                    }
                }
                
                if (elem->on_boundary()){
                    
                    for(uint side_elem = 0; side_elem < elem->n_sides(); ++side_elem){
                        {
                            if (is_tagged_contact_boundary(mesh, elem, side_elem, tags)) {
                                
                                face_set_id[elem->id()].global.push_back(f_id++);
                                
                                offset++;
                            } else {
                                face_set_id[elem->id()].global.push_back(-1);
                            }
                        }
                    }
                }
                else
                {
                    
                    face_set_id[elem->id()].global.insert(face_set_id[elem->id()].global.end(), -1);
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
    static void write_space(const Iterator &begin, const Iterator &end,MeshBase &space, const std::vector<ElementDofMap> &dof_map, const std::vector<ElementDofMap> &variable_number, const std::vector<ElementDofMap> &variable_order, const std::vector<ElementDofMap> &subdomain_id, const std::vector<ElementDofMap> &side_set_id, const std::vector<ElementDofMap> &face_set_id_global,cutk::OutputStream &os)
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
        write_space(begin, end, *m, utopiamesh.dof_map(), utopiamesh.variable_number(), utopiamesh.variable_number(), utopiamesh.subdomain_id(), utopiamesh.side_set_id(), utopiamesh.face_set_id_global(), os); //FIXME
    }
    
    
    static void read_space(cutk::InputStream &is, cutk::shared_ptr<MeshBase> & space,
                           std::vector<ElementDofMap> &dof_map,
                           std::vector<ElementDofMap> &variable_number,
                           std::vector<ElementDofMap> &variable_order,
                           std::vector<ElementDofMap> &subdomain_id,
                           std::vector<ElementDofMap> &side_set_id,
                           std::vector<ElementDofMap> &face_set_id_global,
                           const libMesh::Parallel::Communicator &comm)
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
    
        read_space(is, utopiamesh.utopiamesh()[0], utopiamesh.dof_map(), utopiamesh.variable_number(), utopiamesh.variable_order(), utopiamesh.subdomain_id(), utopiamesh.side_set_id(), utopiamesh.face_set_id_global(), comm_mesh); //FIXME
        
        utopiamesh.set_must_destroy_attached(0,true);
        
    }
    
    
    template<int Dimensions, class Fun>
    static bool SurfaceAssemble(express::Communicator &comm,
                                const std::shared_ptr<MeshBase> &master_slave,
                                const std::shared_ptr<DofMap> &dof_map,
                                const std::shared_ptr<const unsigned int> &_var_num,
                                Fun process_fun,
                                const cutk::Settings &settings,
                                const libMesh::Real search_radius,
                                const std::vector< std::pair<int, int> >  &tags,
                                const bool use_biorth)
    {
        
        std::shared_ptr<UtopiaMesh> local_fun_spaces_new = cutk::make_shared<UtopiaMesh>(master_slave, dof_map, _var_num, tags);
        auto predicate = std::make_shared<cutlibpp::MasterAndSlave>();
        
        for(auto t : tags)
            predicate->add(t.first, t.second);
        
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
        
        const int dim_master = master_slave->mesh_dimension();
        const int dim_slave = master_slave->mesh_dimension();
        
        MeshBase::const_element_iterator e_it = mesh->active_elements_begin();
        const MeshBase::const_element_iterator e_end = mesh->active_elements_end();
        std::vector<int> block_id;
        std::vector<int> block_id_def;
                
        EXPRESS_EVENT_BEGIN("create_adapters");
        ////////////////////////////////////////////////////////////////////////////////////////////////////
        cutk::shared_ptr<NTreeT> tree = NTreeT::New(predicate, maxNElements, maxDepth);
        tree->reserve(n_elements);
        
        std::cout << "nElements = tree->memory().nData()_inside " << n_elements << std::endl;
        
        std::shared_ptr<UtopiaMesh> local_spaces = make_shared<UtopiaMesh>(master_slave, dof_map, _var_num, tags);
        
        int jj=0;
        
        for (auto it = master_slave->active_local_elements_begin(); 
             it != master_slave->active_local_elements_end(); ++it) {
            
            auto elem = *it;
            
            if(!elem->on_boundary()) {
                continue;
            }
            
            bool check_size=false;
            
            for(uint side_elem = 0; side_elem < elem->n_sides(); ++side_elem){
                if ((predicate->select(master_slave->get_boundary_info().boundary_id(elem, side_elem))) && check_size==false){
                    SurfaceAdapter a(*master_slave, elem->id(), elem->id(), master_slave->get_boundary_info().boundary_id(elem, side_elem), search_radius);
                    assert(!local_spaces->dof_map()[elem->id()].empty());
                    a.set_dof_map(&local_spaces->dof_map()[elem->id()].global);
                    a.set_face_id(&local_spaces->face_set_id_global()[elem->id()].global);
                    tree->insert(a);
                    check_size=true;
                    //                    jj++;
                }
            }
        }
                
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
                         const cutk::Settings &settings,
                         const libMesh::Real search_radius,
                         const std::vector< std::pair<int, int> > &tags,
                         const bool use_biorth)
    {
        std::shared_ptr<UtopiaMesh> local_fun_spaces_new = cutk::make_shared<UtopiaMesh>(master_slave, dof_map, _var_num, tags);
        
        libMesh::DenseMatrix<libMesh::Real> src_pts;
        libMesh::DenseMatrix<libMesh::Real> dest_pts;
        libMesh::DenseMatrix<libMesh::Real> intersection2;
        Polyhedron src_poly, dest_poly;
        Polyhedron  intersection3,temp_poly;
        Intersector isector;
        
        std::shared_ptr<MeshBase> master_slave_space = master_slave;
        
        auto predicate = std::make_shared<cutlibpp::MasterAndSlave>();
        
        for(auto t : tags){
            predicate->add(t.first, t.second);
            std::cout << "[Status] Added master-slave pair = "<< t.first << ", " << t.second << std::endl;
        }
        
        static const double tol = 1e-8;
        
        std::vector<libMesh::dof_id_type> master_dofs, slave_dofs;

        libMesh::DenseMatrix<libMesh::Real> elemmat;
        libMesh::DenseMatrix<libMesh::Real> cumulative_elemmat;

        DenseMatrix<Real> side_polygon_master, side_polygon_slave;
        DenseMatrix<Real> isect_polygon_master, isect_polygon_slave;
        
        std::shared_ptr<Transform> trafo_master;
        std::shared_ptr<Transform> trafo_slave;
        
        Point n_master, n_slave;
        
        const int dim = master_slave->mesh_dimension();
        
        libMesh::Real total_intersection_volume = 0.0;
        libMesh::Real local_element_matrices_sum = 0.0;
               
       //all face dof buffers
        express::MapSparseMatrix<double> B_buffer; 
        express::MapSparseMatrix<double> P_buffer;
        express::MapSparseMatrix<double> Q_buffer;
        express::MapSparseMatrix<double> rel_area_buffer;
        express::MapSparseMatrix<double> gap_buffer;
        express::MapSparseMatrix<double> normal_buffer;
        
        std::cout<<"*********** master_slave->dof_map().n_dofs() = "<<  dof_map->n_dofs() <<std::endl;
                
        DenseMatrix<Real> biorth_weights;
        
        auto fun = [&](const SurfaceElementAdapter<Dimensions> &master,
                       const SurfaceElementAdapter<Dimensions> &slave) -> bool {
            
            long n_intersections = 0;
            
            using namespace cutlibpp;
            using namespace express;
            using namespace cutk;
            
            const auto &master_mesh = master.space();
            const auto &slave_mesh  = slave.space();
                        
            libMesh::DenseMatrix<libMesh::Real> elemmat;
            
            const int src_index  = master.element();
            const int dest_index = slave.element();
            
            auto &el_master  = *master_mesh.elem(src_index);
            auto &dest_el = *slave_mesh.elem(dest_index);
            
            const int dim_master = master_mesh.mesh_dimension();
            const int dim_slave = slave_mesh.mesh_dimension();
            
            Box box_master(dim_master), box_slave(dim_slave);
            
            QMortar src_ir_ref(dim_master);
            QMortar src_ir(dim_master);
            QMortar dest_ir(dim_slave);
            QMortar dest_ir_ref(dim_slave);
            
            std::vector<long> order_master = local_fun_spaces_new->variable_order()[0].global;
            const int approx_order = order_master[0];
            
            std::shared_ptr<Contact> surface_assemble;
            
            const auto &side_id_master = master.dof_map_face();
            const auto &face_id_slave  = slave.dof_map_face();
            
            std::unique_ptr<libMesh::FEBase> master_fe, slave_fe;
           
            //
            master_fe = libMesh::FEBase::build(master_mesh.mesh_dimension(), FIRST);
            slave_fe  = libMesh::FEBase::build(slave_mesh.mesh_dimension(),  FIRST);
            //
            
            typedef Intersector::Scalar Scalar;
            
            if(dim_slave == 2)  {
                make_polygon(el_master,   src_pts);
                make_polygon(dest_el,  dest_pts);
                trafo_master  = std::make_shared<Transform2>(el_master);
                trafo_slave = std::make_shared<Transform2>(dest_el);
                
            }
            
            else if(dim_slave == 3) {
                make_polyhedron(el_master,  src_poly);
                make_polyhedron(dest_el, dest_poly);
                trafo_master  = std::make_shared<Transform3>(el_master);
                trafo_slave = std::make_shared<Transform3>(dest_el);
                
            }
            
            bool intersected = false;
            
            for(uint side_index_master = 0; 
                side_index_master < el_master.n_sides(); 
                ++side_index_master) {

                if(side_id_master[side_index_master] < 0) continue;
                
                if(el_master.neighbor_ptr(side_index_master) != nullptr) continue;
                
                auto side_master = el_master.build_side_ptr(side_index_master);
                
                compute_side_normal(dim_master, *side_master, n_master);
                fix_normal_orientation(el_master, side_index_master, n_master);
                
                box_master.reset();
                enlarge_box_from_side(dim_master, *side_master, box_master, search_radius);
                
                if(dim_slave == 2) {
                    make_polygon(*side_master, side_polygon_master);
                } else if(dim_master == 3) {
                    make_polygon_3(*side_master, side_polygon_master);
                } else {
                    assert(false);
                }
                
                for(uint side_2 = 0; side_2 < dest_el.n_sides(); ++side_2) {
                    if(face_id_slave[side_2] < 0) continue;
                    if(dest_el.neighbor_ptr(side_2) != nullptr) continue;
                    // if (!predicate->tagsAreRelated(tag_1, tag_2)) continue;
                    
                    auto side_slave = dest_el.build_side_ptr(side_2);
                    
                    compute_side_normal(dim_slave, *side_slave, n_slave);
                    
                    fix_normal_orientation(dest_el, side_2, n_slave);
                    
                    const Real cos_angle = n_master.contract(n_slave);
                    
                    //if the angle is more than 60 degrees ( cos(60/180*pi) == 0.5 ) or has same orientation skip
                    if(cos_angle >= -0.5) {
                        continue;
                    }
                    
                    box_slave.reset();
                    enlarge_box_from_side(dim_slave, *side_slave, box_slave, search_radius);
                    
                    if(!box_master.intersects(box_slave, tol)) {
                        continue;
                    }
                    
                    
                    bool pair_intersected = false;
                    if(dim_slave == 2){
                        make_polygon(*side_slave, side_polygon_slave);
                        
                        //plot_lines(2, 2, &side_polygon_master.get_values()[0], "in_master/" + std::to_string(master_facq) + "_" + std::to_string(cos_angle));
                        //plot_lines(2, 2, &side_polygon_slave.get_values()[0], "in_slave/" + std::to_string(face_id_slave[0]) + "_" + std::to_string(cos_angle));
                        
                        if(!project_2D(side_polygon_master, side_polygon_slave, isect_polygon_master, isect_polygon_slave)){
                            continue;
                        }


                        const Scalar dx = side_polygon_slave(0, 0) - side_polygon_slave(1, 0);
                        const Scalar dy = side_polygon_slave(0, 1) - side_polygon_slave(1, 1);
                        
                        const Scalar isect_dx = isect_polygon_slave(0, 0) - isect_polygon_slave(1, 0);
                        const Scalar isect_dy = isect_polygon_slave(0, 1) - isect_polygon_slave(1, 1);
                        
                        const Scalar area   = std::sqrt(isect_dx*isect_dx + isect_dy*isect_dy);
                        const Scalar weight = area/std::sqrt(dx*dx + dy*dy);

                        if(weight < 1e-15) continue;
                        
                        const int order = order_for_l2_integral(dim_master, el_master, approx_order, dest_el, approx_order);
                        
                        make_composite_quadrature_on_surf_2D(isect_polygon_master, weight, order, src_ir);
                        
                        make_composite_quadrature_on_surf_2D(isect_polygon_slave, weight, order, dest_ir);
                        
                        pair_intersected = true;
                        
                        surface_assemble = std::make_shared<Contact>();
                        surface_assemble->isect_area	   = area;
                        surface_assemble->relative_area    = weight;
                        
                        // plot_polygon(2, 2, &side_polygon_master.get_values()[0], "master");
                        // plot_polygon(2, 2, &side_polygon_slave.get_values()[0], "slave");
                        
                        
                    } else if(dim_slave == 3) {
                        make_polygon_3(*side_slave, side_polygon_slave);
                        
                        if(!project_3D(
                                       side_polygon_master,
                                       side_polygon_slave,
                                       isect_polygon_master,
                                       isect_polygon_slave))
                        {
                            continue;
                        }
                        
                        const Scalar area_slave = isector.polygon_area_3(side_polygon_slave.m(),  &side_polygon_slave.get_values()[0]);
                        const Scalar area   	= isector.polygon_area_3(isect_polygon_slave.m(), &isect_polygon_slave.get_values()[0]);
                        const Scalar weight 	= area/area_slave;
                        
                        assert(area_slave > 0);
                        assert(area > 0);
                        assert(weight > 0);
                        
                        const int order = order_for_l2_integral(dim_master, el_master, approx_order, dest_el, approx_order);
                        
                        make_composite_quadrature_on_surf_3D(isect_polygon_master, weight, order, src_ir);
                        make_composite_quadrature_on_surf_3D(isect_polygon_slave, weight, order, dest_ir);
                        
                        pair_intersected = true;
                        
                        surface_assemble = std::make_shared<Contact>();
                        surface_assemble->isect_area	= area;
                        surface_assemble->relative_area = weight;
                        
                    } else {
                        assert(false);
                        return false;
                    }
                    
                    
                    if(pair_intersected) {
                        
                        // plot_polygon(3, isect_polygon_master.get_values().size()/3, &isect_polygon_master.get_values()[0], "master");
                        // plot_polygon(3, isect_polygon_slave.get_values().size()/3, &isect_polygon_slave.get_values()[0], "slave");
                        
                        // std::cout << "isect: " << master.handle() << " -> " << slave.handle() << std::endl;
                        
                        //////////////////////////////////ASSEMBLY ////////////////////////////////////////
                        //////////////////////////////////////////////////////////////////////////////////////
                        transform_to_reference_surf(*trafo_master,  el_master.type(),  src_ir, src_ir_ref);
                        transform_to_reference_surf(*trafo_slave, dest_el.type(), dest_ir, dest_ir_ref);
                        
                        master_fe->attach_quadrature_rule(&src_ir_ref);
                        
                        master_fe->get_phi();
                        master_fe->reinit(&el_master);
                        
                        slave_fe->attach_quadrature_rule(&dest_ir_ref);
                        
                        slave_fe->get_xyz();
                        slave_fe->reinit(&dest_el);
                        
                        
                        surface_assemble->parent_element_master  = src_index;
                        
                        surface_assemble->id_master 			 = el_master.id();
                        
                        surface_assemble->parent_element_slave   = dest_index;
                        
                        surface_assemble->id_slave 			     = dest_el.id();
                        
                        surface_assemble->coupling.zero();
                        
                        elemmat.zero();
                        
                        mortar_assemble(*master_fe, *slave_fe, elemmat);
                        
                        std::vector<bool> node_is_boundary_slave;
                        std::vector<bool> node_is_boundary_master;
                        
                        nodes_are_boundary_hack(elemmat, node_is_boundary_slave, node_is_boundary_master);

                        assert( check_node_is_boundary((*master_slave->active_local_elements_begin())->type(), node_is_boundary_master) );
                        
                        if(use_biorth) {
                            assemble_trace_biorth_weights_from_space((*master_slave->active_local_elements_begin())->type(),
                                                                     node_is_boundary_slave,
                                                                     biorth_weights);
                            elemmat.zero();
                            mortar_assemble_weighted_biorth(*master_fe, *slave_fe, biorth_weights, elemmat);
                        }
                        
                        
                        const libMesh::Point pp = side_master->point(0);
                        const Real plane_offset = n_master.contract(pp);
                        
                        
                        if(use_biorth) {
                            mortar_normal_and_gap_assemble_weighted_biorth(
                                                                           *slave_fe,
                                                                           dim,
                                                                           n_slave,
                                                                           n_master,
                                                                           plane_offset,
                                                                           biorth_weights,
                                                                           surface_assemble->normals,
                                                                           surface_assemble->gap);
                        } else {
                            mortar_normal_and_gap_assemble(
                                                           dim,
                                                           *slave_fe,
                                                           n_slave,
                                                           n_master,
                                                           plane_offset,
                                                           surface_assemble->normals,
                                                           surface_assemble->gap);
                        }
                        //////////////////////////////////////////////////////////////////////////////////////
                        
                        rel_area_buffer.add(face_id_slave[side_2], 0, surface_assemble->relative_area);
                        
                        const auto &master_dofs = master.dof_map();
                        const auto &slave_dofs  = slave.dof_map();
                    
                        int n_nodes_face_slave  = side_slave->n_nodes();
                        int n_nodes_face_master = side_master->n_nodes();
                        
                        
                        //plot_lines(dim_master,  side_master->n_nodes(), &side_polygon_master.get_values()[0] , "master/" + std::to_string(side_id_master[0]));
                        //plot_lines(dim_slave, side_slave->n_nodes(), &side_polygon_slave.get_values()[0] , "slave/" + std::to_string(face_id_slave[0]));
                        
                        std::vector<dof_id_type> face_node_id_slave(slave_dofs.size(),   -1);
                        std::vector<dof_id_type> face_node_id_master(master_dofs.size(), -1);
                        
                        //generate face-node ids for slave
                        int offset = 0;
                        for(uint i = 0; i < node_is_boundary_slave.size(); ++i){
                            if (node_is_boundary_slave[i]) {
                                face_node_id_slave[i] =  face_id_slave[side_2] * n_nodes_face_slave + offset++;
                            }
                        }
                        
                        //generate face-node ids for master
                        offset = 0;
                        for(uint i = 0; i < node_is_boundary_master.size(); ++i) {
                            if (node_is_boundary_master[i]) {
                                face_node_id_master[i] = side_id_master[side_index_master] * n_nodes_face_master + offset++;
                            }
                        }
                        
                        //fill-up slave permutation
                        for(int i = 0; i <  slave_dofs.size(); ++i) {
                            const long dof_I = slave_dofs[i];
                            const long dof_J = face_node_id_slave[i];
                            
                            if(node_is_boundary_slave[i]) {
                                P_buffer.setAt(dof_I, dof_J, 1.);
                            }
                        }
                        
                        //fill-up master permutation
                        for(int i = 0; i <  face_node_id_master.size(); ++i) {
                            const long dof_I = master_dofs[i];
                            const long dof_J = face_node_id_master[i];
                            
                            if(node_is_boundary_master[i]) {
                                Q_buffer.setAt(dof_I, dof_J, 1.);
                            }
                        }
                        
                        auto partial_sum = std::accumulate(elemmat.get_values().begin(), elemmat.get_values().end(), libMesh::Real(0.0));
                        assert(partial_sum > 0);

                        local_element_matrices_sum += partial_sum;
                        
                        assert(slave_dofs.size() == elemmat.m());
                        assert(master_dofs.size() == elemmat.n());
                        
                        for(int i = 0; i < face_node_id_slave.size(); ++i) {
                            if(!node_is_boundary_slave[i]) continue;
                            
                            const long dof_I = face_node_id_slave[i];
                            
                            gap_buffer.add(dof_I, 0, surface_assemble->gap(i));
                            
                            for (int k = 0; k < dim_slave; ++k) {
                                normal_buffer.add(dof_I, k, surface_assemble->normals(i,k));
                            }
                            
                            for(int j = 0; j <  face_node_id_master.size(); ++j) {
                                if(!node_is_boundary_master[j]) continue;

                                const long dof_J = face_node_id_master[j];
                                                            
                                B_buffer.add(dof_I, dof_J, elemmat(i, j));
                            }
                        }
                        
                        intersected = true;
                    }
                }
            }
            
            return intersected;
            
        };
        
        if(!SurfaceAssemble<Dimensions>(comm, master_slave, dof_map, _var_num, fun, settings, search_radius, tags, use_biorth)) {
            return false;
        }
        
        
        double volumes[1] = { local_element_matrices_sum };
        
        comm.allReduce(volumes, 1, express::MPISum());
        
        const processor_id_type master_proc_id  = master_slave->processor_id();
        
        const dof_id_type n_dofs_on_proc_master = dof_map->n_dofs_on_processor(master_proc_id);
        
        const processor_id_type slave_proc_id   = master_slave->processor_id();
        
        const dof_id_type n_dofs_on_proc_slave  = dof_map->n_dofs_on_processor(slave_proc_id);
        
        if(comm.isRoot()) {
            std::cout << "sum(B_tilde): " << volumes[0] <<std::endl;
        }

        express::Array<express::SizeType>  ownership_ranges_master(comm.size()+1);
        ownership_ranges_master.allSet(0);
                
        express::Array<express::SizeType>  ownership_ranges_slave(comm.size()+1);
        ownership_ranges_slave.allSet(0);
        
        const int n_nodes_x_face = master_slave->elem(0)->build_side_ptr(0)->n_nodes();
        express::Array<express::SizeType>  side_node_ownership_ranges = local_fun_spaces_new->ownershipRangesFaceID();

        for(SizeType i = 0; i < side_node_ownership_ranges.size(); ++i) {
            side_node_ownership_ranges[i] *= n_nodes_x_face;
        }
        
        ownership_ranges_master[comm.rank() + 1] += static_cast<unsigned int>(n_dofs_on_proc_master);
        ownership_ranges_slave[comm.rank()  + 1] += static_cast<unsigned int>(n_dofs_on_proc_slave);
        
        comm.allReduce(&ownership_ranges_master[0], ownership_ranges_master.size(), express::MPISum());
        comm.allReduce(&ownership_ranges_slave[0],  ownership_ranges_slave.size(),  express::MPISum());
        
        std::partial_sum(ownership_ranges_master.begin(), ownership_ranges_master.end(),
                         ownership_ranges_master.begin());
        
        std::partial_sum(ownership_ranges_slave.begin(), ownership_ranges_slave.end(),
                         ownership_ranges_slave.begin());
        
        
        const SizeType n_local_dofs_slave     = ownership_ranges_slave [comm.rank() + 1] - ownership_ranges_slave [comm.rank()];
        const SizeType n_local_dofs_master    = ownership_ranges_master[comm.rank() + 1] - ownership_ranges_master[comm.rank()];
        const SizeType n_local_side_node_dofs = side_node_ownership_ranges[comm.rank() + 1] - side_node_ownership_ranges[comm.rank()];
        
        B_buffer.finalizeStructure();
        P_buffer.finalizeStructure();
        Q_buffer.finalizeStructure();
        rel_area_buffer.finalizeStructure();
        gap_buffer.finalizeStructure();
        normal_buffer.finalizeStructure();
        
        SizeType sizes[6] = {
            B_buffer.rows(), B_buffer.columns(),
            P_buffer.rows(), P_buffer.columns(),
            Q_buffer.rows(), Q_buffer.columns()
        };
        
        comm.allReduce(sizes, 6, express::MPIMax());
        
        const SizeType n_side_node_dofs = express::Math<SizeType>::Max(sizes[0], sizes[1], sizes[3], sizes[5]);
        std::cout << "n_side_node_dofs = " << n_side_node_dofs << std::endl;
        
        B_buffer.setSize(n_side_node_dofs, n_side_node_dofs);
        P_buffer.setSize(sizes[2], n_side_node_dofs);
        Q_buffer.setSize(sizes[4], n_side_node_dofs);
        rel_area_buffer.setSize(n_side_node_dofs, 1);
        gap_buffer.setSize(n_side_node_dofs, 1);
        normal_buffer.setSize(n_side_node_dofs, dim);
        
        express::Redistribute< express::MapSparseMatrix<double> > redist(comm.getMPIComm());
        redist.apply(side_node_ownership_ranges, B_buffer,     express::AddAssign<double>());
        redist.apply(side_node_ownership_ranges, gap_buffer,       express::AddAssign<double>());
        redist.apply(side_node_ownership_ranges, normal_buffer,    express::AddAssign<double>());
        
        redist.apply(local_fun_spaces_new->ownershipRangesFaceID(), rel_area_buffer, express::AddAssign<double>());
        
        redist.apply(ownership_ranges_slave, P_buffer, express::Assign<double>());
        redist.apply(ownership_ranges_master, Q_buffer, express::Assign<double>());
        
        std::cout << "n_local_side_node_dofs: " << n_local_side_node_dofs << std::endl;
       
        express::RootDescribe("petsc rel_area_buffer assembly begin", comm, std::cout);
        
        express::Array<bool> remove_row(n_local_side_node_dofs);
                
        if(!remove_row.isNull()) {
            long n_remove_rows = 0;
            remove_row.allSet(false);
            
            {
                for (auto it = rel_area_buffer.iter(); it; ++it) {
                    if(*it < 1 - 1e-8) {
                        const SizeType faceId = it.row();

                        for(int k = 0; k < n_nodes_x_face; ++k) {
                            const SizeType nodeId = faceId * n_nodes_x_face + k;
                            const SizeType index  = nodeId - side_node_ownership_ranges[comm.rank()];
                            assert(index < remove_row.size());
                            remove_row[index] = true;
                            ++n_remove_rows;
                        }
                    }
                }
            }
            
            std::cout << "n_remove_rows: " <<n_remove_rows << std::endl;
        }
        
        express::RootDescribe("petsc B_buffer assembly begin", comm, std::cout);
        
        DVectord is_contact_node_tilde = local_zeros(n_local_side_node_dofs);
        DVectord gap_tilde = local_zeros(n_local_side_node_dofs);
        {
            utopia::Write<utopia::DVectord> write(gap_tilde);
            for (auto it = gap_buffer.iter(); it; ++it) {
                
                const SizeType index = it.row() - side_node_ownership_ranges[comm.rank()];
                assert(index < remove_row.size());
                                
                if(!remove_row[index]) {
                    gap_tilde.set(it.row(), *it);
                    is_contact_node_tilde.set(it.row(), 1.0);
                }
            }
        }
        
        DSMatrixd normal_tilde = utopia::local_sparse(n_local_side_node_dofs, dim, dim);
        {
            utopia::Write<utopia::DSMatrixd> write(normal_tilde);
            for (auto it = normal_buffer.iter(); it; ++it) {
                
                const SizeType index = it.row() - side_node_ownership_ranges[comm.rank()];
                assert(index < remove_row.size());
                            
                if(!remove_row[index]) {
                    normal_tilde.set(it.row(), it.col(), *it);
                }
            }
        }

        SizeType n_max_row_entries_bpq[3] = { B_buffer.maxEntriesXCol(), P_buffer.maxEntriesXCol(), Q_buffer.maxEntriesXCol() };
        comm.allReduce(n_max_row_entries_bpq, 3, express::MPIMax());

        const SizeType n_max_row_entries_b = n_max_row_entries_bpq[0];
        const SizeType n_max_row_entries_p = n_max_row_entries_bpq[1];
        const SizeType n_max_row_entries_q = n_max_row_entries_bpq[2];

        DSMatrixd B_tilde = utopia::local_sparse(n_local_side_node_dofs, n_local_side_node_dofs, n_max_row_entries_b);
        {
            utopia::Write<utopia::DSMatrixd> write(B_tilde);
            for (auto it = B_buffer.iter(); it; ++it) {
                
                const SizeType index = it.row() - side_node_ownership_ranges[comm.rank()];
                assert(index < remove_row.size());
                                
                if(!remove_row[index]) {
                    B_tilde.set(it.row(), it.col(), *it);
                }
            }
        }
        
        comm.barrier();
        express::RootDescribe("petsc P_buffer assembly begin", comm, std::cout);
         
        DSMatrixd P = utopia::local_sparse(n_local_dofs_slave, n_local_side_node_dofs, n_max_row_entries_p);
        {
            utopia::Write<utopia::DSMatrixd> write(P);
            for (auto it = P_buffer.iter(); it; ++it) {
                P.set(it.row(), it.col(), *it);
            }
        }
        
        comm.barrier();
        express::RootDescribe("petsc Q_buffer assembly begin", comm, std::cout);
        
        DSMatrixd Q = utopia::local_sparse(n_local_dofs_master, n_local_side_node_dofs, n_max_row_entries_q);
        {
            utopia::Write<utopia::DSMatrixd> write(Q);
            for (auto it = Q_buffer.iter(); it; ++it) {
                Q.set(it.row(), it.col(), *it);
            }
        }
                
        DSMatrixd Q_t = transpose(Q);
        DSMatrixd B_x = P * B_tilde * Q_t;
        
        normals = P * normal_tilde;
        
        DVectord gap_x = P * gap_tilde;
        
        is_contact_node = P * is_contact_node_tilde;
        

        DVectord normals_vec = local_zeros(n_local_dofs_slave);
        {
            Write<DVectord> w(normals_vec);
            
            each_read(normals, [&](const SizeType i, const SizeType j, const double value){
                normals_vec.set(i + j, value);
            });
        }
        
        
        bool has_contact = false;
        
        each_read(is_contact_node, [&](const SizeType i , const double value){
            if (value > 0)
            {
                is_contact_node.set(i, 1);
                has_contact = true;
            }
        });
        
        orthogonal_trafos = local_sparse(n_local_dofs_slave , n_local_dofs_slave , dim);        
        {
            typedef Intersector::Scalar Scalar;
            std::vector<Scalar> normal(dim, 0);
            std::vector<Scalar> H(dim * dim, 0);
            
            Read<DVectord>  r_n(normals_vec);
            Write<DSMatrixd> w_o(orthogonal_trafos);
            Read<DVectord> r_icn(is_contact_node);
            
            bool check_has_contact = false;
            
            utopia::Range r = utopia::range(normals_vec);
            for(uint i = r.begin(); i < r.end(); i += dim) {
                bool use_identity = true;
                bool is_cn_i = is_contact_node.get(i) > 0;
                
                if(is_cn_i) {                    
                    check_has_contact = true;
                    
                    for(uint d = 0; d < dim; ++d) {
                        normal[d] = normals_vec.get(i + d);
                    }
                    
                    normalize(normal);
                    
                    if(std::abs(normal[0] - 1.) > 1e-8) {
                        use_identity = false;
                        
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
                                orthogonal_trafos.set((i + di), (i + dj), H[di * dim + dj]);
                            }
                        }
                    }
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
        
        auto size_B_x = local_size(B_x);
        B = local_sparse(size_B_x.get(0), size_B_x.get(1), n_max_row_entries_b * dim);
        
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
        
        comm.barrier();
        express::RootDescribe("petsc assembly end", comm, std::cout);
        return true;
    }
    
    

    
// }




   inline bool MooseSurfaceAssemble(express::Communicator &comm,
                                    const std::shared_ptr<MeshBase> &mesh,
                                    const std::shared_ptr<DofMap> &dof_map,
                                    const std::shared_ptr<const unsigned int> &_var_num,
                                    DSMatrixd &B,
                                    DSMatrixd &orthogonal_trafos,
                                    DVectord &gap,
                                    DSMatrixd &normals,
                                    DVectord &is_contact_node,
                                    const libMesh::Real search_radius,
                                    const std::vector< std::pair<int, int> > &tags,
                                    const bool use_biorth = true)
   {
       cutk::Settings settings;
       if(mesh->mesh_dimension() == 2) {
           return utopia::SurfaceAssemble<2>(comm, mesh, dof_map, _var_num, B,  orthogonal_trafos, gap, normals, is_contact_node, settings, search_radius, tags, use_biorth);
       }


       if(mesh->mesh_dimension() == 3) {
           return utopia::SurfaceAssemble<3>(comm, mesh, dof_map, _var_num, B,  orthogonal_trafos, gap, normals,  is_contact_node, settings, search_radius, tags, use_biorth);
       }

       assert(false && "Dimension not supported!");
       return false;
   }


   inline bool MooseSurfaceAssemble(express::Communicator &comm,
                                    const std::shared_ptr<MeshBase> &mesh,
                                    const std::shared_ptr<DofMap> &dof_map,
                                    const std::shared_ptr<const unsigned int> &_var_num,
                                    DSMatrixd &B,
                                    DSMatrixd &orthogonal_trafos,
                                    DVectord &gap,
                                    DSMatrixd &normals,
                                    DVectord &is_contact_node,
                                    const libMesh::Real search_radius,
                                    const int tag_1, 
                                    const int tag_2,
                                    const bool use_biorth = true)
   {
        return MooseSurfaceAssemble(
            comm, mesh, dof_map, _var_num, 
            B, orthogonal_trafos, 
            gap, normals, is_contact_node, 
            search_radius,
            { {tag_1, tag_2} },
            use_biorth);
    }
}

#endif //LIBMESH_CUTLIBPP_ADAPTERS_HPP



