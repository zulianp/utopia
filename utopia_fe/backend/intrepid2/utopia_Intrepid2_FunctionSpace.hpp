#ifndef UTOPIA_INTREPID2_TAG_FUNCTION_SPACE_HPP
#define UTOPIA_INTREPID2_TAG_FUNCTION_SPACE_HPP

#include "utopia_FunctionSpace.hpp"
#include "utopia_Intrepid2_FEForwardDeclarations.hpp"
#include "utopia_Intrepid2_Types.hpp"
#include "utopia_Traits.hpp"

//#include "libmesh/dof_map.h"
/*#include "libmesh/equation_systems.h"
#include "libmesh/linear_implicit_system.h"*/
#include "libmesh/reference_counter.h"

#include <memory>

namespace utopia {

    class Intrepid2FunctionSpace final : public FunctionSpace<Intrepid2FunctionSpace> {
    public:
        inline explicit Intrepid2FunctionSpace(
            // libMesh::System &equation_system,
            const int var_num)
        //: equation_systems_(
        //    make_ref(equation_system.get_equation_systems())
        //  ),
        // system_num_(equation_system.number())
        {
            this->set_subspace_id(var_num);
        }

        inline explicit Intrepid2FunctionSpace(
            //   const std::shared_ptr<libMesh::EquationSystems> &equation_systems,
            const int system_num,
            const int var_num)
        //  : equation_systems_(equation_systems), system_num_(system_num)
        {
            this->set_subspace_id(var_num);
        }

        inline explicit Intrepid2FunctionSpace(
            //  const std::shared_ptr<libMesh::EquationSystems> &equation_systems,
            //  const libMesh::FEFamily &type = libMesh::LAGRANGE,
            //  const libMesh::Order &order = libMesh::FIRST,
            const std::string &var_name = "",
            const int system_num = 0)
            //: equation_systems_(equation_systems),
            : system_num_(system_num) {
            std::string var_name_copy = var_name;

            if (var_name_copy.empty()) {
                // var_name_copy = "var_" + std::to_string(equation_system().n_vars());
            }

            const int var_num = 4;  // equation_system().add_variable(var_name_copy, order, type);

            // assert(equation_system().n_vars() > 0);

            this->set_subspace_id(var_num);
        }

        inline explicit Intrepid2FunctionSpace(
            //  libMesh::MeshBase &mesh,
            //  const libMesh::FEFamily &type = libMesh::LAGRANGE,
            //  const libMesh::Order &order = libMesh::FIRST,
            const std::string &var_name = "")
        // : equation_systems_(std::make_shared<libMesh::EquationSystems>(mesh)), system_num_(0)
        {
            //   equation_systems_->add_system<libMesh::LinearImplicitSystem>("main");

            std::string var_name_copy = var_name;

            if (var_name_copy.empty()) {
                //   var_name_copy = "var_" + std::to_string(equation_system().n_vars());
            }

            const int var_num = 4;  // equation_system().add_variable(var_name_copy, order, type);

            // assert(equation_system().n_vars() > 0);

            this->set_subspace_id(var_num);
        }

        // inline bool is_initialized() const
        // {
        // 	return equation_system().is_initialized();
        // }

        inline void initialize() {
            //   if(!equation_system().is_initialized()) {
            //       equation_system().init();
            // equation_systems_->init();
            //   }
        }

        /*inline libMesh::Order order(const int) const
        {
            return dof_map().variable_order(this->subspace_id());
        }

        inline libMesh::FEType type()
        {
            return dof_map().variable_type(this->subspace_id());
        }*/

        //     inline libMesh::DofMap &dof_map() {
        //         return equation_system().get_dof_map();
        //     }

        //     inline const libMesh::DofMap &dof_map() const {
        //         return equation_system().get_dof_map();
        //     }

        /*        inline libMesh::System &equation_system()
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
                }*/

        //   inline libMesh::MeshBase &mesh() { return equation_systems_->get_mesh(); }
        //   inline const libMesh::MeshBase &mesh() const { return equation_systems_->get_mesh(); }

        inline std::string get_class() const override { return "Intrepid2FunctionSpace"; }

    private:
        // std::shared_ptr<libMesh::EquationSystems> equation_systems_; //TODO
        int system_num_;
    };

    /*template<>
    class Traits<Intrepid2FunctionSpace> : public LibMeshAlgebraTraits<double> { //TODO
    public:
        static const int Backend = LIBMESH_TAG;
        static const int Order = 1;
        static const int FILL_TYPE = FillType::DENSE;

        typedef double Scalar;
        typedef utopia::LMDenseVector Vector;
        typedef utopia::LMDenseMatrix Matrix;
        typedef libMesh::TensorValue<Scalar> TensorValueT;
        typedef libMesh::VectorValue<Scalar> VectorValueT;

        typedef utopia::LibMeshFunctionSpace Implementation;

        typedef libMesh::FEBase FE;

        typedef std::vector<std::vector<libMesh::FEBase::OutputShape>> FunctionType;
        typedef std::vector<std::vector<libMesh::FEBase::OutputGradient>> GradientType;
        typedef std::vector<std::vector<libMesh::FEBase::OutputDivergence>> DivergenceType;
        typedef std::vector<std::vector<TensorValueT>> JacobianType;
        typedef std::vector<std::vector<VectorValueT>> CurlType;

        typedef std::vector<libMesh::Real> DXType;
        typedef libMesh::MeshBase MeshType;







    };*/

    typedef utopia::Traits<Intrepid2FunctionSpace> Intrepid2Traits;

    inline auto elements_begin(const libMesh::MeshBase &m) -> decltype(m.active_local_elements_begin()) {
        return m.active_local_elements_begin();
    }

    inline auto elements_end(const libMesh::MeshBase &m) -> decltype(m.active_local_elements_end()) {
        return m.active_local_elements_end();
    }

    void write(const Path &path, Intrepid2FunctionSpace &space, UVector &x);
}  // namespace utopia

#endif  // UTOPIA_INTREPID2_TAG_FUNCTION_SPACE_HPP
