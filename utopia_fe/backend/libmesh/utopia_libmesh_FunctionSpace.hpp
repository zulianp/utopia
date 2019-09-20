#ifndef UTOPIA_LIBMESH_TAG_FUNCTION_SPACE_HPP
#define UTOPIA_LIBMESH_TAG_FUNCTION_SPACE_HPP

#include "utopia_libmesh_FEForwardDeclarations.hpp"
#include "utopia_FunctionSpace.hpp"
#include "utopia_Traits.hpp"
#include "utopia_libmesh_Types.hpp"

#include "libmesh/equation_systems.h"
#include "libmesh/dof_map.h"
#include "libmesh/reference_counter.h"
#include "libmesh/linear_implicit_system.h"

#include <memory>

namespace utopia {

    class LibMeshFunctionSpace final : public FunctionSpace<LibMeshFunctionSpace> {
    public:
        inline explicit LibMeshFunctionSpace(
        libMesh::System &equation_system,
        const int var_num)
        : equation_systems_(
            make_ref(equation_system.get_equation_systems())
          ),
         system_num_(equation_system.number())
         {
             this->set_subspace_id(var_num);
         }


        inline explicit LibMeshFunctionSpace(
            const std::shared_ptr<libMesh::EquationSystems> &equation_systems,
            const int system_num,
            const int var_num)
        : equation_systems_(equation_systems), system_num_(system_num)
        {
            this->set_subspace_id(var_num);
        }

        inline explicit LibMeshFunctionSpace(
            const std::shared_ptr<libMesh::EquationSystems> &equation_systems,
            const libMesh::FEFamily &type = libMesh::LAGRANGE,
            const libMesh::Order &order = libMesh::FIRST,
            const std::string &var_name = "",
            const int system_num = 0)
        : equation_systems_(equation_systems),
          system_num_(system_num)
        {
            std::string var_name_copy = var_name;

            if(var_name_copy.empty()) {
                var_name_copy = "var_" + std::to_string(equation_system().n_vars());
            }

            const int var_num = equation_system().add_variable(var_name_copy, order, type);

            assert(equation_system().n_vars() > 0);

            this->set_subspace_id(var_num);
        }


        inline explicit LibMeshFunctionSpace(
            libMesh::MeshBase &mesh,
            const libMesh::FEFamily &type = libMesh::LAGRANGE,
            const libMesh::Order &order = libMesh::FIRST,
            const std::string &var_name = "")
        : equation_systems_(std::make_shared<libMesh::EquationSystems>(mesh)), system_num_(0)
        {
            equation_systems_->add_system<libMesh::LinearImplicitSystem>("main");

            std::string var_name_copy = var_name;

            if(var_name_copy.empty()) {
                var_name_copy = "var_" + std::to_string(equation_system().n_vars());
            }

            const int var_num = equation_system().add_variable(var_name_copy, order, type);

            assert(equation_system().n_vars() > 0);

            this->set_subspace_id(var_num);
        }

        // inline bool is_initialized() const
        // {
        // 	return equation_system().is_initialized();
        // }

        inline void initialize()
        {
            if(!equation_system().is_initialized()) {
                equation_system().init();
                // equation_systems_->init();
            }
        }

        inline libMesh::Order order(const int) const
        {
            return dof_map().variable_order(this->subspace_id());
        }

        inline libMesh::Order order() const
        {
            return dof_map().variable_order(this->subspace_id());
        }


        inline libMesh::FEType type()
        {
            return dof_map().variable_type(this->subspace_id());
        }

        inline libMesh::DofMap &dof_map() {
            return equation_system().get_dof_map();
        }

        inline const libMesh::DofMap &dof_map() const {
            return equation_system().get_dof_map();
        }

        inline libMesh::System &equation_system()
        {
            return equation_systems_->get_system(system_num_);
        }

        inline const libMesh::System &equation_system() const
        {
            return equation_systems_->get_system(system_num_);
        }

        inline libMesh::EquationSystems &equation_systems() {
            return *equation_systems_;
        }

        inline const libMesh::EquationSystems &equation_systems() const
        {
            return *equation_systems_;
        }

        inline const std::shared_ptr<libMesh::EquationSystems> &equation_systems_ptr() const
        {
            return equation_systems_;
        }

        inline libMesh::MeshBase &mesh() { return equation_systems_->get_mesh(); }
        inline const libMesh::MeshBase &mesh() const { return equation_systems_->get_mesh(); }

        inline std::string get_class() const override {
            return "LibMeshFunctionSpace";
        }

    private:
        std::shared_ptr<libMesh::EquationSystems> equation_systems_;
        int system_num_;
    };

    template<>
    class Traits<LibMeshFunctionSpace> //: public LibMeshAlgebraTraits<double> 
    {
    public:
        static const int Backend = LIBMESH_TAG;
        static const int Order = 1;
        static const int FILL_TYPE = FillType::DENSE;

        typedef double Scalar;
        typedef utopia::LMDenseVector Vector;
        typedef utopia::LMDenseMatrix Matrix;
        // typedef libMesh::TensorValue<Scalar> TensorValueT;
        // typedef libMesh::VectorValue<Scalar> VectorValueT;

        typedef utopia::LibMeshFunctionSpace Implementation;

        typedef libMesh::FEBase FE;

        typedef std::vector<std::vector<double>> FunctionType;
        typedef std::vector<std::vector<utopia::LMDenseVector>> GradientType;
        typedef std::vector<std::vector<double>> DivergenceType;
        typedef std::vector<std::vector<utopia::LMDenseMatrix>> JacobianType;
        typedef std::vector<std::vector<utopia::LMDenseVector>> CurlType;

        typedef std::vector<double> DXType;
        typedef libMesh::MeshBase MeshType;
        typedef utopia::LibMeshAssemblyValues AssemblyValues;
    };

    typedef utopia::Traits<LibMeshFunctionSpace> LibMeshTraits;

    inline auto elements_begin(const libMesh::MeshBase &m) -> decltype(m.active_local_elements_begin())
    {
        return m.active_local_elements_begin();
    }

    inline auto elements_end(const libMesh::MeshBase &m) -> decltype(m.active_local_elements_end())
    {
        return m.active_local_elements_end();
    }


    void write(const Path &path, LibMeshFunctionSpace &space, UVector &x);

    std::size_t max_nnz_x_row(const libMesh::DofMap &dof_map);
    std::size_t max_nnz_x_row(const LibMeshFunctionSpace &space);
}

#endif //UTOPIA_LIBMESH_TAG_FUNCTION_SPACE_HPP
