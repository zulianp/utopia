#ifndef LibmeshTransferForMoose_HPP
#define LibmeshTransferForMoose_HPP

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
#include "libmesh/serial_mesh.h"

#include "express_Profiler.hpp"
#include "express_Redistribute.hpp"

#include "Box.hpp"
#include "MapSparseMatrix.hpp"
#include "utopia_fe.hpp"
#include "MortarAssemble.hpp"
#include "libmesh/serial_mesh.h"

#include <cmath>
#include <queue>



namespace utopia {
    using namespace libMesh;
    
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
    };
    
    
    template<int _Dimension>
    class TreeTraits {
    public:
        enum {
            Dimension = _Dimension
        };
        
        typedef utopia::BoxBoxAdapter<Dimension> Bound;
        typedef utopia::ElementAdapter<Dimension> DataType;
        
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
    
    
    class Spaces {
    public:
        
        
        explicit Spaces(const express::Communicator &comm) : comm(comm)
        {
            must_destroy_attached[0] = false;
            must_destroy_attached[1] = false;
        }
        
        Spaces(const std::shared_ptr<MeshBase> &master,
               const std::shared_ptr<MeshBase> &slave,
               const std::shared_ptr<libMesh::DofMap>  &dof_map_master,
               const std::shared_ptr<libMesh::DofMap>  &dof_map_slave,
               const std::shared_ptr<const unsigned int> &_from_var_num,
               const std::shared_ptr<const unsigned int> &_to_var_num)
        {
            
            spaces_.reserve(2);
            spaces_.push_back(master);
            spaces_.push_back(slave);
            
            must_destroy_attached[0] = false;
            must_destroy_attached[1] = false;
            
            const int n_elements_master = master->n_elem();
            const int n_elements_slave  = slave->n_elem();
            
            const int n_elements = n_elements_master + n_elements_slave;
            
            std::cout<<"MASTER DOF"<<std::endl;
            copy_global_dofs(*master, dof_map_master, _from_var_num, dof_maps_[0], var_type_[0], n_elements);
            
            std::cout<<"SLAVE DOF"<<std::endl;
            copy_global_dofs(*slave,  dof_map_slave, _to_var_num, dof_maps_[1], var_type_[1], n_elements);
            
            //            copy_var_number(*master, var_number_[0]);
            //            copy_var_number(*slave,  var_number_[1]);
            //
            copy_var_order(*dof_map_master, var_order_[0]);
            copy_var_order(*dof_map_slave,  var_order_[1]);
            
        }
        
        
        inline std::vector< std::shared_ptr<MeshBase> > &spaces()
        {
            return spaces_;
            
        }
        
        inline const std::vector< std::shared_ptr<MeshBase> > &spaces() const
        {
            return spaces_;
            
        }
        
        inline long n_elements() const
        {
            long ret = 0;
            for(auto s : spaces_) {
                if(s) {
                    ret += s->n_elem();
                }
            }
            
            return ret;
        }
        
        
        inline std::vector<ElementDofMap> &dof_map(const int i)
        {
            assert(i < 2);
            assert(i >= 0);
            return dof_maps_[i];
        }
        
        inline const std::vector<ElementDofMap> &dof_map(const int i) const
        {
            assert(i < 2);
            assert(i >= 0);
            return dof_maps_[i];
        }
        
        inline void set_must_destroy_attached(const int index, const bool value)
        {
            assert(index < 2);
            assert(index >= 0);
            must_destroy_attached[index] = value;
        }
        
        
        inline  std::vector<ElementDofMap> &variable_number(const int i)
        {
            assert(i < 2);
            assert(i >= 0);
            return var_number_[i];
        }
        
        inline const std::vector<ElementDofMap> &variable_number(const int i) const
        {
            assert(i < 2);
            assert(i >= 0);
            return var_number_[i];
        }
        
        
        
        inline std::vector<ElementDofMap> &variable_order(const int i)
        {
            assert(i < 2);
            assert(i >= 0);
            return var_order_[i];
        }
        
        
        
        inline const std::vector<ElementDofMap> &variable_order(const int i) const
        {
            assert(i < 2);
            assert(i >= 0);
            return var_order_[i];
        }
        
        
        inline std::vector<ElementDofMap> &variable_type(const int i)
        {
            assert(i < 2);
            assert(i >= 0);
            return var_type_[i];
        }
        
        
        
        inline const std::vector<ElementDofMap> &variable_type(const int i) const
        {
            assert(i < 2);
            assert(i >= 0);
            return var_type_[i];
        }
        
        
    private:
        express::Communicator comm;
        std::vector<std::shared_ptr< MeshBase>> spaces_;
        std::vector<ElementDofMap> dof_maps_[2];
        std::vector<ElementDofMap> var_number_[2];
        std::vector<ElementDofMap> var_order_[2];
        std::vector<ElementDofMap> var_type_[2];
        bool must_destroy_attached[2];
        
        
        
        
        inline static void copy_global_dofs(MeshBase &space, const std::shared_ptr<libMesh::DofMap>  &original_dof_map, const std::shared_ptr<const unsigned int>  &var_num,
                                            std::vector<ElementDofMap> &dof_map, std::vector<ElementDofMap> &variable_type, const int n_elements)
        {
            
            //            auto &mesh = space.get_mesh();
            //            auto &original_dof_map = space.get_dof_map();
            std::vector<dof_id_type> temp;
            std::shared_ptr<ElementDofMap> temp_ptr;
            dof_map.resize(n_elements);
            
            
            //        std::cout<<"______________________________COPY_DOF_BEGIN____________________________"<<std::endl;
            
            MeshBase::const_element_iterator e_it        =space.active_local_elements_begin();
            const MeshBase::const_element_iterator e_end =space.active_local_elements_end();
            
            int i=0;
            
            variable_type.resize(1);
            
            bool first=true;
            
            
            for (; e_it != e_end; ++e_it){
                
                Elem *elem = *e_it;
                original_dof_map->dof_indices(elem, temp, *var_num);
                
                dof_map[elem->id()].global.insert(dof_map[elem->id()].global.end(), temp.begin(), temp.end());
                
                if (first)
                {
                    variable_type[0].global.push_back(elem->type());
                    first=false;
                    
                }
                i++;
            }
            
            
            
            // std::cout<<"______________________________COPY_DOF_END____________________________"<<std::endl;
            
        }
        
        
        
        //        inline static void copy_var_number(MeshBase &space, std::vector<ElementDofMap> &variable_number)
        //        {
        //            variable_number.resize(1);
        //            variable_number[0].global.push_back(space.var_num());
        //        }
        
        inline static void copy_var_order(DofMap &dofmap, std::vector<ElementDofMap> &variable_order)
        {
            variable_order.resize(1);
            FEType fe_type =  dofmap.variable(0).type();
            variable_order[0].global.push_back(fe_type.order);
        }
        
        
    };
    
    
    template<class Iterator>
    static void write_space(const Iterator &begin, const Iterator &end,MeshBase &space,
                            const std::vector<ElementDofMap> &dof_map,/* const std::vector<ElementDofMap> &variable_number,*/
                            const std::vector<ElementDofMap> &variable_order,const int role, cutk::OutputStream &os)
    {
        const int dim 		  = space.mesh_dimension();
        const long n_elements = std::distance(begin, end);
        
        std::set<long> nodeIds;
        std::map<long, long> mapping;
        
        //        std::cout<<"------------------------------------write_space-IN-------------------------------------------"<<std::endl;
        
        
        
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
        os << dim << role;
        
        
        int index = 0;
        for (auto nodeId : nodeIds) {
            mapping[nodeId] = index++;
        }
        
        //WRITE 2
        os << n_nodes;
        
        //WRITE 6
        os << n_elements;
        
        //        std::cout<<"write_n_el = "<<n_elements<<std::endl;
        
        for(auto node_id : nodeIds){
            
            const Point &p = space.node(node_id);
            
            for(int i = 0; i < dim; ++i) {
                
                //WRITE 3
                os << p(i);
                
            }
            
            //std::cout<<"write_point"<<p<<std::endl;
            
        }
        
        std::vector<dof_id_type> indices_vector;
        
        CHECK_STREAM_WRITE_BEGIN("elements", os);
        
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
            
            
            
        }
        
        CHECK_STREAM_WRITE_END("elements", os);
        
        //       std::cout<<"------------------------------------write_space-OUT------------------------------------------"<<std::endl;
        
        
        //WRITE 10
        //        os << variable_number.at(0);
        
        //WRITE 11
        os << variable_order.at(0);
        
        
        
    }
    
    
    
    
    template<class Iterator>
    static void write_element_selection(const Iterator &begin, const Iterator &end, const Spaces &spaces, cutk::OutputStream &os)
    {
        //       std::cout<<"------------------------------------WRITE_ELEM_SELEC-IN------------------------------------------"<<std::endl;
        
        
        if(spaces.spaces().empty()){
            assert(false);
            return;
        }
        
        
        
        auto m = spaces.spaces()[0];
        assert(m);
        std::shared_ptr<MeshBase> s =nullptr;
        
        if(spaces.spaces().size()>1) {
            s=spaces.spaces()[1];
        }
        
        
        
        unsigned int dim_master=m->mesh_dimension();
        
        std::vector<long> master_selection;
        std::vector<long> slave_selection;
        
        bool met_slave_selection = false;
        
        //       std::cout<<"------------------------------------WRITE_ELEM_SELEC-MID---------------------------------------"<<std::endl;
        
        for(Iterator it = begin; it != end; ++it) {
            int index =*it;
            
            if(m && index >= m->n_elem()) {
                index -= m->n_elem();
                slave_selection.push_back(index);
            }
            
            else if(!m) {
                met_slave_selection = true;
                slave_selection.push_back(index);
            }
            
            else {
                assert(!met_slave_selection);
                assert(index < m->n_elem());
                master_selection.push_back(index);
            }
        }
        
        
        
        const bool has_master = !master_selection.empty();
        const bool has_slave  = !slave_selection.empty();
        
        os << has_master << has_slave;
        
        
        if(has_master) {
            
            //            std::cout<<"I am in master"<<std::endl;
            write_space(master_selection.begin(), master_selection.end(), *m, spaces.dof_map(0),
                        /*spaces.variable_number(0),*/ spaces.variable_order(0), 0, os);
        }
        
        if(has_slave) {
            //            std::cout<<"I am in slave"<<std::endl;
            write_space(slave_selection.begin(), slave_selection.end(), *s, spaces.dof_map(1),
                        /*spaces.variable_number(1),*/ spaces.variable_order(1),1, os);
        }
        
        //   std::cout<<"------------------------------------WRITE_ELEM_SELEC--OUT----------------------------------------"<<std::endl;
        
        
    }
    //
    
    static void read_space(cutk::InputStream &is, cutk::shared_ptr<MeshBase> & space,
                           std::vector<ElementDofMap> &dof_map, /*std::vector<ElementDofMap> &variable_number,*/
                           std::vector<ElementDofMap> &variable_order, const libMesh::Parallel::Communicator &comm)
    {
        
        //   std::cout<<"------------------------------------READ-SPACE-IN-------------------------------------------"<<std::endl;
        
        
        using namespace std;
        
        
        
        //READ 1
        int dim, role;
        is >> dim >> role;
        
        //READ 2
        long n_nodes;
        is >> n_nodes;
        
        //READ 6
        long n_elements;
        is >> n_elements;
        
        auto mesh_ptr = std::make_shared<SerialMesh>(comm, dim);
        
        //      std::cout<<"read_n_el = "<<n_elements<<std::endl;
        //
        //      EquationSystems equation_systems (*mesh_ptr);
        //
        //      LinearImplicitSystem & system = equation_systems.add_system<LinearImplicitSystem> ("Serial");
        //
        //      std::cout<<"ciao2 = "<<n_elements<<std::endl;
        
        //      equation_systems.init();
        
        //      equation_systems.print_info();
        
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
        
        
        CHECK_STREAM_READ_BEGIN("elements", is);
        
        for(long i = 0; i !=n_elements; ++i) {
            
            //READ 7
            
            int type, e_n_nodes;
            
            is >> type >> e_n_nodes;
            
            auto elem = Elem::build(ElemType(type)).release();
            
            
            int index;
            
            for (int ii = 0; ii != e_n_nodes; ++ii) {
                
                //READ 8
                is >> index;
                
                elem->set_node(ii) = & mesh_ptr->node(index);
                
            }
            
            //READ 9
            is >> dof_map.at(i);
            
            mesh_ptr->add_elem(elem);
            
            libmesh_assert(elem);
            
            //std::cout<<"Local Elem"<<mesh_ptr->n_local_elem()<<std::endl;
            
            
            
            
        }
        
        CHECK_STREAM_READ_END("elements", is);
        
        //READ 10
        //  variable_number.resize(1);
        //  is >> variable_number.at(0);
        
        //READ 11
        variable_order.resize(1);
        is >> variable_order.at(0);
        
        
        //!!!! dummy parameters
        
        
        
        
        //  space = make_shared<System>(equation_systems,"Serial",0);
        
        space = mesh_ptr;
        
        
        
        //        std::cout<<"ciao= "<<n_elements<<std::endl;
        
        //   std::cout<<"------------------------------------READ-SPACE-OUT-------------------------------------------"<<std::endl;
        
    }
    
    static void read_spaces(cutk::InputStream &is, Spaces &spaces, const libMesh::Parallel::Communicator &comm_master, const libMesh::Parallel::Communicator &comm_slave)
    {
        //    std::cout<<"------------------------------------READ-SPACES-IN-------------------------------------------"<<std::endl;
        
        bool has_master, has_slave;
        is >> has_master >> has_slave;
        
        spaces.spaces().resize(2);
        
        
        if(has_master) {
            read_space(is, spaces.spaces()[0], spaces.dof_map(0),
                       /* spaces.variable_number(0),*/spaces.variable_order(0), comm_master);
            spaces.set_must_destroy_attached(0,true);
        } else {
            spaces.spaces()[0] = nullptr;
            
        }
        
        if(has_slave) {
            read_space(is, spaces.spaces()[1], spaces.dof_map(1),
                       /*spaces.variable_number(1),*/spaces.variable_order(1),comm_slave);
            spaces.set_must_destroy_attached(1,true);
        } else {
            spaces.spaces()[1] = nullptr;
            
        }
        
        //       std::cout<<"------------------------------------READ-SPACES-OUT------------------------------------------"<<std::endl;
        
    }
    
    template<int Dimensions, class Fun>
    static bool Assemble(express::Communicator &comm,
                         const std::shared_ptr<MeshBase> &master,
                         const std::shared_ptr<MeshBase> &slave,
                         const std::shared_ptr<DofMap> &dof_master,
                         const std::shared_ptr<DofMap> &dof_slave,
                         const std::shared_ptr<const unsigned int> &_from_var_num,
                         const std::shared_ptr<const unsigned int> &_to_var_num,
                         Fun process_fun,
                         const cutk::Settings &settings, bool use_biorth_, bool impact)
    {
        
        
        using namespace cutlibpp;
        using namespace express;
        using namespace cutk;
        
        typedef LibMeshTree<Dimensions> NTreeT;
        typedef typename NTreeT::DataContainer DataContainer;
        typedef typename NTreeT::DataType Adapter;
        
        long maxNElements = 100;
        long maxDepth = 5;
        
        
        if (!settings.get("max_depth").isNull()) {
            maxDepth = settings.get("max_depth").toInt();
            std::cout<<"max_depth  = "<< maxDepth  <<std::endl;
        }
        
        const auto &master_mesh = master;
        const auto &slave_mesh  = slave;
        const int n_elements_master = master_mesh->n_elem();
        const int n_elements_slave  = slave_mesh->n_elem();
        const int n_elements 		= n_elements_master + n_elements_slave;
        
        
        const Parallel::Communicator &libmesh_comm_master = master_mesh->comm();
        const Parallel::Communicator &libmesh_comm_slave = slave_mesh->comm();
        
        
        auto predicate = make_shared<MasterAndSlave>();
        predicate->add(0, 1);
        
        EXPRESS_EVENT_BEGIN("create_adapters");
        ////////////////////////////////////////////////////////////////////////////////////////////////////
        cutk::shared_ptr<NTreeT> tree = NTreeT::New(predicate, maxNElements, maxDepth);
        std::cout<<"n_elements"<<n_elements<<std::endl;
        tree->reserve(n_elements);
        
        
        std::shared_ptr<Spaces> local_spaces = make_shared<Spaces>(master, slave, dof_master, dof_slave, _from_var_num, _to_var_num);
        
        int offset = 0;
        int space_num = 0;
        
        for(auto s : local_spaces->spaces()) {
            if(s) {
                
                bool first = true;
                for (auto it = s->active_local_elements_begin(); it != s->active_local_elements_end(); ++it) {
                    auto elem=*it;
                    Adapter a(*s, elem->id(), offset+elem->id(), space_num);
                    assert(!local_spaces->dof_map(space_num)[elem->id()].empty());
                    a.set_dof_map(&local_spaces->dof_map(space_num)[elem->id()].global);
                    tree->insert(a);
                }
                
                offset += s->n_elem(); //s->mesh().n_active_local_elem();//(*s->mesh().active_local_elements_end())->id();
                
            }
            
            
            ++space_num;
            
            
        }
        
        //std::cout<<" --------------------------- tree->memory().nData()=" << tree->memory().nData()<<std::endl; ;
        
        tree->getRoot()->getBound().staticBound().enlarge(1e-8);
        
        ////////////////////////////////////////////////////////////////////////////////////////////////////
        EXPRESS_EVENT_END("create_adapters");
        
        // std::cout<<"-----------------------------------ADAPTERS-------------------------------------------------"<<std::endl;
        //Just to have an indexed-storage
        std::map<long, cutk::shared_ptr<Spaces> > spaces;
        std::map<long, std::vector<cutk::shared_ptr<Spaces> > > migrated_spaces;
        
        
        auto read = [&spaces, &migrated_spaces, comm, &libmesh_comm_master, &libmesh_comm_slave ]
        (
         const long ownerrank,
         const long senderrank,
         bool is_forwarding, DataContainer &data,
         InputStream &in
         ) {
            
            
            //   std::cout<<"------------------------------------AUTO-READ-IN---------------------------------------------"<<std::endl;
            
            
            CHECK_STREAM_READ_BEGIN("vol_proj", in);
            
            cutk::shared_ptr<Spaces> proc_space = cutk::make_shared<Spaces>(comm);
            
            read_spaces(in, *proc_space, libmesh_comm_master, libmesh_comm_slave);
            
            if (!is_forwarding) {
                assert(!spaces[ownerrank]);
                spaces[ownerrank] = proc_space;
            } else {
                migrated_spaces[ownerrank].push_back(proc_space);
            }
            
            data.reserve(data.size() + 3000);
            //std::cout<<proc_space->dof_map(0)[0].global<<std::endl;
            
            int space_num = 0;
            long offset = 0;
            for(auto s : proc_space->spaces()) {
                if(s) {
                    for (int i=0; i<s->n_elem(); i++) {
                        data.push_back(Adapter(*s, i, offset + i, space_num) );
                        assert(!proc_space->dof_map(space_num)[i].empty());
                        data.back().set_dof_map(&proc_space->dof_map(space_num)[i].global);
                    }
                    
                    offset += s->n_elem();
                    
                }
                
                ++space_num;
                
            }
            
            
            //    std::cout<<"------------------------------------AUTO-READ-OUT---------------------------------------------"<<std::endl;
            
            CHECK_STREAM_READ_END("vol_proj", in);
            
            
            
        };
        
        
        auto write = [&local_spaces, &spaces, &comm]
        (
         const long ownerrank, const long recvrank,
         const std::vector<long>::const_iterator &begin,
         const std::vector<long>::const_iterator &end,
         const DataContainer &data,
         OutputStream &out) {
            
            CHECK_STREAM_WRITE_BEGIN("vol_proj", out);
            
            //       std::cout<<"------------------------------------AUTO-WRITE-IN-------------------------------------------"<<std::endl;
            
            if (ownerrank == comm.rank()) {
                
                write_element_selection(begin, end, *local_spaces, out);
                
                
            } else {
                
                auto it = spaces.find(ownerrank);
                assert(it != spaces.end());
                cutk::shared_ptr<Spaces> spaceptr = it->second;
                assert(std::distance(begin, end) > 0);
                write_element_selection(begin, end, *spaceptr, out);
                
            }
            
            //      std::cout<<"------------------------------------AUTO-WRITE-OUT-------------------------------------------"<<std::endl;
            //      comm.barrier();
            
            CHECK_STREAM_WRITE_END("vol_proj", out);
            
        };
        
        
        long n_false_positives = 0, n_intersections = 0;
        
        //       std::cout<<"------------------------------------SEARCH-COMPUTE-------------------------------------------"<<std::endl;
        
        auto fun = [&n_false_positives, &n_intersections, &process_fun](
                                                                        
                                                                        Adapter &master, Adapter &slave) -> bool {
            
            bool ok = process_fun(master, slave);
            
            if(ok) {
                n_intersections++;
                
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
        
        
        
        long n_total_candidates = n_intersections + n_false_positives;
        
        long n_collection[3] = {n_intersections, n_total_candidates, n_false_positives};
        
        
        comm.allReduce(n_collection, 3, express::MPISum());
        
        
        if (comm.isRoot()) {
            std::cout << "n_intersections: " << n_collection[0]
            << ", n_total_candidates: " 	 << n_collection[1]
            << ", n_false_positives: " 	     << n_collection[2] << std::endl;
        }
        
        return true;
    }
    
    static void assemble_biorth_weights_from_space(const std::shared_ptr<MeshBase> &mesh,
                                                   const std::shared_ptr<DofMap> &dof_map,
                                                   const int var_num,
                                                   libMesh::DenseMatrix<libMesh::Real> &weights)
    {
        const int dim = mesh->mesh_dimension();
        std::unique_ptr<libMesh::FEBase> biorth_elem =
        libMesh::FEBase::build(dim,
                               dof_map->variable_type(var_num));
        
        auto &el = **mesh->active_local_elements_begin();
        
        const int order = order_for_l2_integral(dim,
                                                el, dof_map->variable(var_num).type().order,
                                                el, dof_map->variable(var_num).type().order);
        
        libMesh::QGauss qg(dim, libMesh::Order(order));
        biorth_elem->attach_quadrature_rule(&qg);
        biorth_elem->reinit(&el);
        mortar_assemble_weights(*biorth_elem, weights);
    }
    
    template<int Dimensions>
    bool Assemble(
                  express::Communicator &comm,
                  const std::shared_ptr<MeshBase> &master,
                  const std::shared_ptr<MeshBase> &slave,
                  const std::shared_ptr<DofMap> &dof_master,
                  const std::shared_ptr<DofMap> &dof_slave,
                  const std::shared_ptr<const unsigned int> &_from_var_num,
                  const std::shared_ptr<const unsigned int> &_to_var_num,
                  DSMatrixd &B1,/* DSMatrixd &B2,*/
                  const cutk::Settings &settings,bool  use_biorth_, bool  impact)
    {
        
        const int var_num_slave = *_to_var_num;
        
        std::shared_ptr<Spaces> local_fun_spaces = cutk::make_shared<Spaces>(master, slave, dof_master, dof_slave,_from_var_num,_to_var_num);
        
        libMesh::DenseMatrix<libMesh::Real> src_pts;
        libMesh::DenseMatrix<libMesh::Real> dest_pts;
        libMesh::DenseMatrix<libMesh::Real> intersection2;
        Polyhedron src_poly, dest_poly;
        Polyhedron  intersection3,temp_poly;
        Intersector isector;
        
        std::shared_ptr<MeshBase> master_space = master;
        std::shared_ptr<MeshBase> slave_space  = slave;
        
        
        // std::vector<libMesh::dof_id_type> master_dofs, slave_dofs;
        libMesh::DenseMatrix<libMesh::Real> elemmat;
        libMesh::DenseMatrix<libMesh::Real> cumulative_elemmat;
        
        std::shared_ptr<Transform> src_trans;
        std::shared_ptr<Transform> dest_trans;
        
        
        int skip_zeros = 1;
        
        
        libMesh::Real total_intersection_volume = 0.0;
        libMesh::Real local_element_matrices_sum = 0.0;
        
        
        
        express::MapSparseMatrix<double> mat_buffer(dof_slave->n_dofs(), dof_master->n_dofs());
        
        //        std::cout<<"dof_slave->n_dofs()"<<dof_slave->n_dofs()<<std::endl;
        //        std::cout<<"dof_master->n_dofs()"<<dof_master->n_dofs()<<std::endl;
        
        bool intersected = false;
        
        double element_setup_time = 0.0;
        double intersection_time = 0.0;
        double assembly_time     = 0.0;
        
        utopia::Chrono c;
        
        auto fun = [&](const ElementAdapter<Dimensions> &master,
                       const ElementAdapter<Dimensions> &slave) -> bool {
            
            c.start();
            
            libMesh::DenseMatrix<libMesh::Real> biorth_weights;
            
            if(use_biorth_) {
                assemble_biorth_weights_from_space(slave_space,
                                                   dof_slave,
                                                   var_num_slave,
                                                   biorth_weights);
            }
            
            long n_intersections = 0;
            
            bool pair_intersected = false;
            
            const auto &src  = master.space();
            
            const auto &dest = slave.space();
            
            const auto &src_mesh  = src;
            
            const auto &dest_mesh = dest;
            
            const int src_index  = master.element();
            
            const int dest_index = slave.element();
            
            auto &src_el  = *src_mesh.elem(src_index);
            
            auto &dest_el = *dest_mesh.elem(dest_index);
            
            const int dim = src_mesh.mesh_dimension();
            
            
            std::unique_ptr<libMesh::FEBase> master_fe, slave_fe;
            
            master_fe = libMesh::FEBase::build(src_mesh.mesh_dimension(),  dof_master->variable_type(0));
            slave_fe  = libMesh::FEBase::build(dest_mesh.mesh_dimension(), dof_slave->variable_type(0));
            
            QMortar composite_ir(dim);
            QMortar src_ir(dim);
            QMortar dest_ir(dim);
            
            
            const int order = order_for_l2_integral(dim, src_el, dof_master->variable(0).type().order , dest_el,dof_slave->variable(0).type().order);
            
            c.stop();
            element_setup_time += c.get_seconds();
            c.start();
            
            if(dim == 2)  {
                make_polygon(src_el,   src_pts);
                make_polygon(dest_el, dest_pts);
                
                if(intersect_2D(src_pts, dest_pts, intersection2)) {
                    total_intersection_volume += fabs(isector.polygon_area_2(intersection2.m(), &intersection2.get_values()[0]));
                    
                    const libMesh::Real weight=isector.polygon_area_2(dest_pts.m(), &dest_pts.get_values()[0]);
                    
                    make_composite_quadrature_2D(intersection2, weight, order, composite_ir);
                    pair_intersected = true;
                    
                    src_trans  = std::make_shared<Transform2>(src_el);
                    dest_trans = std::make_shared<Transform2>(dest_el);
                    pair_intersected = true;
                }
            }
            else if(dim == 3) {
                make_polyhedron(src_el,  src_poly);
                make_polyhedron(dest_el, dest_poly);
                
                
                if(intersect_3D(src_poly, dest_poly, intersection3)) {
                    
                    total_intersection_volume += isector.p_mesh_volume_3(intersection3);
                    
                    const libMesh::Real weight = isector.p_mesh_volume_3(dest_poly);
                    
                    make_composite_quadrature_3D(intersection3, weight, order, composite_ir);
                    src_trans  = std::make_shared<Transform3>(src_el);
                    dest_trans = std::make_shared<Transform3>(dest_el);
                    pair_intersected = true;
                }
                
            } else {
                assert(false);
                return false;
            }
            
            c.stop();
            intersection_time += c.get_seconds();
            c.start();
            
            const auto &master_dofs = master.dof_map();
            const auto &slave_dofs  = slave.dof_map();
            
            if(pair_intersected) {
                
                
                transform_to_reference(*src_trans,  src_el.type(),  composite_ir,  src_ir);
                transform_to_reference(*dest_trans, dest_el.type(), composite_ir,  dest_ir);
                
                //            src.dof_map().dof_indices(&src_el,  master_dofs);
                //            dest.dof_map().dof_indices(&dest_el, slave_dofs);
                
                
                
                
                assert(!master_dofs.empty());
                assert(!slave_dofs.empty());
                //composite_ir.print_info();
                
                
                master_fe->attach_quadrature_rule(&src_ir);
                master_fe->reinit(&src_el);
                
                slave_fe->attach_quadrature_rule(&dest_ir);
                slave_fe->reinit(&dest_el);
                
                elemmat.zero();
                
                
                
                if(use_biorth_) {
                    // mortar_assemble_biorth(*master_fe, *slave_fe, dest_el.type(), elemmat);
                    //std::cout<<"I am here"<<std::endl;
                    mortar_assemble_weighted_biorth(*master_fe, *slave_fe, biorth_weights, elemmat);
                    
                } else {
                    mortar_assemble(*master_fe, *slave_fe, elemmat);
                }
                
                // std::cout << "-----------------------------------------\n";
                // std::cout << src_index << ", " << dest_index << "\n";
                // elemmat.print(std::cout);
                // for(auto i : slave_dofs) {
                // 	std::cout << i << " ";
                // }
                // std::cout << "\n";
                
                // for(auto i : master_dofs) {
                // 	std::cout << i << " ";
                // }
                // std::cout << "\n";
                // std::cout << "-----------------------------------------\n";
                
                auto partial_sum = std::accumulate(elemmat.get_values().begin(), elemmat.get_values().end(), libMesh::Real(0.0));
                // std::cout << src_index << ", " << dest_index << ": " << partial_sum << std::endl;
                // dest_ir.print_info();
                
                local_element_matrices_sum += partial_sum;
                
                intersected = true;
                
                ++n_intersections;
                
                
                //                if(slave_dofs.size() != elemmat.m()) {
                //                    std::cout << slave_dofs.size() << " != " <<  elemmat.m() << std::endl;
                //                }
                
                assert(slave_dofs.size() == elemmat.m());
                assert(master_dofs.size() == elemmat.n());
                
                // std::cout<<"slave_dofs.size()"<<slave_dofs.size()<<std::endl;
                // std::cout<<"master_dofs.size()"<<master_dofs.size()<<std::endl;
                
                for(int i = 0; i < slave_dofs.size(); ++i) {
                    
                    const long dof_I = slave_dofs[i];
                    
                    for(int j = 0; j < master_dofs.size(); ++j) {
                        
                        const long dof_J = master_dofs[j];
                        
                        mat_buffer.add(dof_I, dof_J, elemmat(i, j));
                    }
                }
                
                return true;
                
            } else {
                
                return false;
            }
            
        };
        
        
        
        // comm.barrier();
        // utopia::Chrono c2;
        // c2.start();
        
        
        if(!Assemble<Dimensions>(comm, master, slave, dof_master, dof_slave, _from_var_num, _to_var_num, fun, settings, use_biorth_,impact)) {
            std::cout << "n_intersections: false2" <<std::endl;
            return false;
        }
        
        // c2.stop();
        // std::cout << "Local mortars" << std::endl;
        // c2.describe(std::cout);
        
        
        // comm.barrier();
        // std::stringstream ss;
        // ss << "Setup_time: " << element_setup_time << "\n";
        // ss << "intersection_time: " << intersection_time << "\n";
        // ss << "assembly_time: " << assembly_time << std::endl;
        
        //express::SynchDescribe(ss.str(), comm, std::cout);
        // comm.barrier();
        // c2.start();
        
        std::cout << "n_intersections: " <<std::endl;
        
        double volumes[2] = { local_element_matrices_sum,  total_intersection_volume };
        
        comm.allReduce(volumes, 2, express::MPISum());
        
        const processor_id_type master_proc_id  = master->processor_id();
        
        const dof_id_type n_dofs_on_proc_master = dof_master->n_local_dofs();
        
        const processor_id_type slave_proc_id   = slave->processor_id();
        
        const dof_id_type n_dofs_on_proc_slave  =dof_slave->n_local_dofs();
        
        const int n_dofs_on_proc_print  = dof_slave->n_local_dofs();
        
        if(comm.isRoot()) {
            std::cout << "sum(B): " << volumes[0] << ", vol(I): " << volumes[1] << std::endl;
        }
        
        express::Array<express::SizeType>  ownershipRangesMaster(comm.size()+1);
        ownershipRangesMaster.allSet(0);
        
        express::Array<express::SizeType>  ownershipRangesSlave(comm.size()+1);
        ownershipRangesSlave.allSet(0);
        
        
        ownershipRangesMaster[comm.rank()+1] += static_cast<unsigned int>(n_dofs_on_proc_master);
        
        ownershipRangesSlave[comm.rank()+1]  += static_cast<unsigned int>(n_dofs_on_proc_slave);
        
        
        
        comm.allReduce(&ownershipRangesMaster[0], ownershipRangesMaster.size(), express::MPISum());
        
        comm.allReduce(&ownershipRangesSlave[0],  ownershipRangesSlave.size(),  express::MPISum());
        
        
        std::partial_sum(ownershipRangesMaster.begin(), ownershipRangesMaster.end(), ownershipRangesMaster.begin());
        std::partial_sum(ownershipRangesSlave.begin(), ownershipRangesSlave.end(), ownershipRangesSlave.begin());
        
        
        //        if(comm.isRoot()) {
        //            std::cout <<ownershipRangesMaster << std::endl;
        //            std::cout<<"prova"<<n_dofs_on_proc_print<<std::endl;
        //
        //        }
        
        
        
        int dim = master->mesh_dimension();
        
        express::Redistribute< express::MapSparseMatrix<double> > redist(comm.getMPIComm());
        
        redist.apply(ownershipRangesSlave, mat_buffer, express::AddAssign<double>());
        
        assert(ownershipRangesSlave.empty() == ownershipRangesMaster.empty() || ownershipRangesMaster.empty());
        
        express::RootDescribe("petsc assembly begin", comm, std::cout);
        
        SizeType  mMaxRowEntries = mat_buffer.maxEntriesXCol();
        
        comm.allReduce(&mMaxRowEntries, 1, express::MPIMax());
        
        const SizeType local_range_slave_range  = ownershipRangesSlave [comm.rank()+1] - ownershipRangesSlave [comm.rank()];
        const SizeType local_range_master_range = ownershipRangesMaster[comm.rank()+1] - ownershipRangesMaster[comm.rank()];
        
        B1 = utopia::local_sparse(local_range_slave_range, local_range_master_range, mMaxRowEntries);
        
        {
            utopia::Write<utopia::DSMatrixd> write(B1);
            for (auto it = mat_buffer.iter(); it; ++it) {
                B1.set(it.row(), it.col(), *it);
                
            }
        }

        
        
//        if (impact)
//        {
//            auto s_B_x = local_size(B1);
//            B2 = local_sparse(s_B_x.get(0), s_B_x.get(1), mMaxRowEntries);
//            
//            {
//                std::cout<< "i am modifying the impact  matrix"<<std::endl;
//                utopia::Write<DSMatrixd> w_B(B2);
//                utopia::each_read(B1, [&](const utopia::SizeType i, const utopia::SizeType j, const double value) {
//                    for(utopia::SizeType d = 0; d < dim; ++d) {
//                        //std::cout<< "i am modifying the matrix with value"<<value<<std::endl;
//                        B2.set(i, j, 0.0);
//                        B2.set(i, j+d, value);
//                    }
//                });
//                
//                
//            }
//        }
//        else
//        {
//            auto s_B_x = local_size(B1);
//          
//            B2 = local_sparse(s_B_x.get(0), s_B_x.get(1), mMaxRowEntries);
//            
//            
//            std::cout<< "i am modifying the matrix"<<std::endl;
//            utopia::Write<DSMatrixd> w_B(B2);
//            utopia::each_read(B1, [&](const utopia::SizeType i, const utopia::SizeType j, const double value) {
//                for(utopia::SizeType d = 0; d < dim; ++d) {
//                   // std::cout<< "i am modifying the matrix with value"<<value<<std::endl;
//                    B2.set(i, j, 0.0);
//                    B2.set(i+d, j, value);
//                }
//            });
//            
//            
//            
//            
//        }
        
        
        
        express::RootDescribe("petsc assembly end", comm, std::cout);
        
//        write("_B_1.m", B1);
//        write("_B_2.m", B2);
        
        // c2.stop();
        // std::cout << "Global stuff\n";
        // c2.describe(std::cout);
        return true;
        
    }
    
    
    
    inline bool AssembleMOOSE(express::Communicator &comm,
                              const std::shared_ptr<MeshBase> &master,
                              const std::shared_ptr<MeshBase> &slave,
                              const std::shared_ptr<DofMap> &dof_master,
                              const std::shared_ptr<DofMap> &dof_slave,
                              const std::shared_ptr<const unsigned int> & _from_var_num,
                              const std::shared_ptr<const unsigned int> & _to_var_num,
                              bool  use_biorth_, bool impact,
                              DSMatrixd &B1/*, DSMatrixd &B2*/)
    {
        cutk::Settings settings;
        
        if(master->mesh_dimension() == 2) {
            std::cout<<"Assemble_matrix::I am in assemble"<<std::endl;
            return utopia::Assemble<2>(comm, master, slave, dof_master, dof_slave, _from_var_num,  _to_var_num, B1, /*B2,*/ settings,use_biorth_, impact);
        }
        
        if(master->mesh_dimension() == 3) {
            std::cout<<"Assemble_matrix::I am in assemble"<<std::endl;
            return utopia::Assemble<3>(comm, master, slave, dof_master, dof_slave, _from_var_num,  _to_var_num, B1, /*B2,*/ settings,use_biorth_, impact);
        }
        
        assert(false && "Dimension not supported!");
        return false;
    }
    //
    //
    ////    bool AssembleMOOSE(express::Communicator &comm,
    ////                       const std::shared_ptr<MeshBase> &mesh_master,
    ////                       const std::shared_ptr<MeshBase> &mesh_slave,
    ////                       libMesh::Order master_order,
    ////                       libMesh::Order slave_order,
    ////                       DSMatrixd &B)
    ////    {
    ////        cutk::Settings settings;
    ////
    ////
    ////
    ////        LibMeshFEContext<LinearImplicitSystem> master_context(mesh_master);
    ////        auto master_space = fe_space(LAGRANGE, master_order, master_context);
    ////        master_context.equation_systems.init();
    ////
    ////        LibMeshFEContext<LinearImplicitSystem> slave_context(mesh_slave);
    ////        auto slave_space = fe_space(LAGRANGE, slave_order, slave_context);
    ////        slave_context.equation_systems.init();
    ////
    ////
    ////        if(mesh_master->mesh_dimension() == 2) {
    ////            std::cout<<"Assemble_matrix::I am in assemble"<<std::endl;
    ////            return utopia::Assemble<2>(comm, make_ref(master_space), make_ref(slave_space), B, settings);
    ////        }
    ////
    ////        if(mesh_master->mesh_dimension() == 3) {
    ////            std::cout<<"Assemble_matrix::I am in assemble"<<std::endl;
    ////            return utopia::Assemble<3>(comm, make_ref(master_space), make_ref(slave_space), B, settings);
    ////        }
    ////
    ////        assert(false && "Dimension not supported!");
    ////        return false;
    ////    }
    ////
    //
    //
    //
    //
    //
}

#endif //LIBMESH_CUTLIBPP_ADAPTERS_HPP



