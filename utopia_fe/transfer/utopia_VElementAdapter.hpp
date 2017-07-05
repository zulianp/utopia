#ifndef UTOPIA_V_ELEMENT_ADAPTER_HPP
#define UTOPIA_V_ELEMENT_ADAPTER_HPP 

#include "cutk_Serializable.hpp"
#include "cutk_InputStream.hpp"
#include "cutk_OutputStream.hpp"

#include "utopia_BoxAdapter.hpp"

#include "libmesh/serial_mesh.h"

namespace utopia {

template<int Dimension>
    class VElementAdapter : public cutk::Serializable {
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
        
        VElementAdapter(libMesh::MeshBase &fe, const libMesh::dof_id_type &element, const long element_handle, const int tag)
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
        
        
        
        VElementAdapter()
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
        
        
        inline const  libMesh::MeshBase  &space() const
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
        

         void set_dof_map_reverse(std::vector<long> * ptr)
        {
            
            dof_map_reverse_ = ptr;
        }

        inline const std::vector<long> &dof_map_reverse() const
        {
            assert(dof_map_reverse_);
            return *dof_map_reverse_;
        }
           
    private:
        libMesh::MeshBase * fe_;
        libMesh::dof_id_type element_;
        long element_handle_;
        int tag_;
        BoxBoxAdapter<Dimension> bound_;
        std::vector<long> * dof_map_;
        std::vector<long> * dof_map_reverse_;
    };
}
#endif //UTOPIA_V_ELEMENT_ADAPTER_HPP