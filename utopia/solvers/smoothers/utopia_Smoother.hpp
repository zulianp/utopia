// #ifndef UTOPIA_SMOOTHER_HPP
// #define UTOPIA_SMOOTHER_HPP
// #include "utopia_Core.hpp"
// #include "utopia_Clonable.hpp"
// #include "utopia_Input.hpp"

// #include <iomanip>

//      namespace utopia
//      {
//         template<class Vector>
//         class Smoother : public virtual Clonable, public virtual Configurable
//         {
//             typedef UTOPIA_SCALAR(Vector)           Scalar;
//             typedef UTOPIA_SIZE_TYPE(Vector)        SizeType;

//         public:

//         /**
//          * @brief      Base class for smoothers.
//          */
//         Smoother()
//         : _sweeps(1), _relaxation_parameter(1.)
//         {

//         }

//         virtual ~Smoother() {}

//         /**
//          * @brief      Single sweep. Function needs to be provided by actual
//          smoothers.
//          * @return
//          */
//         virtual bool smooth(const Vector &rhs, Vector &x) = 0;

//         // virtual void update(const std::shared_ptr<const Operator<Vector> >
//         &) = 0;

//         /**
//          * @brief      Quick interface for smoothing with projecting
//          constraints.
//          */
//         // virtual bool nonlinear_smooth(const Matrix &/*A*/, const Vector
//         &/*rhs*/, const Vector& /*ub*/, const Vector& /*lb*/, Vector &/*x*/,
//         std::vector<SizeType>& /*zero_rows*/){ return 0; }

//         /**
//          * @brief      Get number of sweeps.
//          *
//          * @return
//          */
//         virtual SizeType sweeps()
//         {
//             return _sweeps;
//         }

//         /**
//          * @brief      Set the sweeps.
//          *
//          * @param[in]  sweeps   The number of sweeps.
//          *
//          * @return
//          */
//         virtual void sweeps(const SizeType & sweeps_in)
//         {
//             _sweeps = sweeps_in;
//         }

//         /**
//          * @brief      Set omega.
//          *
//          * @return     omega  The relaxation parameter.
//          */
//         virtual Scalar relaxation_parameter()
//         {
//             return _relaxation_parameter;
//         }

//         /**
//          * @brief      Set omega.
//          *
//          * @return     omega  The relaxation parameter.
//          */
//         virtual void relaxation_parameter(const Scalar &
//         relaxation_parameter)
//         {
//             _relaxation_parameter = relaxation_parameter;
//         }

//         virtual Smoother * clone() const override = 0;

//         virtual void read(Input &in) override
//         {
//             in.get("sweeps", _sweeps);
//             in.get("relaxation_parameter", _relaxation_parameter);
//         }

//         virtual void print_usage(std::ostream &os) const override
//         {
//             this->print_param_usage(os, "sweeps", "int", "Number of smoothing
//             steps.", "1.0"); this->print_param_usage(os,
//             "relaxation_parameter", "real", "Relaxation parameter.", "1.0");
//         }

//     private:
//         SizeType     _sweeps;
//         Scalar       _relaxation_parameter;
// };

// }

// #endif //UTOPIA_SMOOTHER_HPP
