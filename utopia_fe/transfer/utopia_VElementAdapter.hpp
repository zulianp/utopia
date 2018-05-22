#ifndef UTOPIA_V_ELEMENT_ADAPTER_HPP
#define UTOPIA_V_ELEMENT_ADAPTER_HPP 

#include "moonolith_serializable.hpp"
#include "moonolith_input_stream.hpp"
#include "moonolith_output_stream.hpp"

#include "utopia_BoxAdapter.hpp"
#include "MortarAssemble.hpp"

#include "libmesh/serial_mesh.h"

namespace utopia {

template<int Dimension>
    class VElementAdapter : public moonolith::Serializable {
    public:
        inline int tag() const
        {
            return tag_;
        }
        
        const BoxBoxAdapter<Dimension> &bound() const
        {
            return bound_;
        }
        
        BoxBoxAdapter<Dimension> &bound()
        {
            return bound_;
        }
        
        void apply_read_write(moonolith::Stream &stream) override
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

            bool is_s = is_shell();

            libMesh::Point n;

            if(is_s) {
                compute_side_normal(Dimension, *e, n);
                assert(n.size()> 0.99);
            }
            
            std::array<double, Dimension> p_a;
            for (libMesh::dof_id_type i = 0; i < e->n_nodes(); ++i) {
                const libMesh::Point &p = fe.node(e->node(i));
                for(int d = 0; d < Dimension; ++d) {
                    p_a[d] = p(d);
                }

                if(is_s) {
                    std::array<double, Dimension> p_a_plus;
                    std::array<double, Dimension> p_a_minus;

                    for(int d = 0; d < Dimension; ++d) {
                        p_a_plus[d]  = p_a[d] + n(d) * 1e-10;
                        p_a_minus[d] = p_a[d] - n(d) * 1e-10;
                    }

                    bound_.static_bound()  += p_a_minus;
                    bound_.dynamic_bound() += p_a_minus; 

                    bound_.static_bound()  += p_a_plus;
                    bound_.dynamic_bound() += p_a_plus; 

                } else {
                    bound_.static_bound()  += p_a;
                    bound_.dynamic_bound() += p_a; 
                }
            }

            assert(!bound_.static_bound().empty());
            assert(!bound_.dynamic_bound().empty());
        }
        
        VElementAdapter()
        : fe_(nullptr) , element_(-1), element_handle_(-1), tag_(-1), dof_map_(nullptr) {}
        
        
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


        bool is_shell() const {
            assert(fe_);

            libMesh::Elem &e = *fe_->elem(element_);
            
            if(Dimension == 3) {
                return is_tri(e.type()) || is_quad(e.type());

            } else if(Dimension == 2) {
                return !is_tri(e.type()) && !is_quad(e.type());
            }

            return false;
        }
    };
}

#endif //UTOPIA_V_ELEMENT_ADAPTER_HPP
