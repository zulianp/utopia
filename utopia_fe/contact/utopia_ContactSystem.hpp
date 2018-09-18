#ifndef UTOPIA_CONTACT_SYSTEM_HPP
#define UTOPIA_CONTACT_SYSTEM_HPP 

#include "utopia_Mechanics.hpp"
#include "utopia_Contact.hpp"
#include "utopia.hpp"

#include <memory>


namespace libMesh {
	class EquationSystems;
}

namespace utopia {

	class ContactSystem {
	public:
		ContactSystem(const std::shared_ptr<libMesh::EquationSystems> &equation_systems, const int main_system_num);
		~ContactSystem();

		void update(
			const MechanicsState &state,
			const Contact &contact,
			const double dt = 1.
			);

		void set_wear_coefficient(const double wear_coefficient)
		{
			wear_coefficient_ = wear_coefficient;
		}

	private:
		//init aux system for plotting
		void init();
		std::shared_ptr<libMesh::EquationSystems> equation_systems_;
		int system_num_;
		std::vector<int> var_num_;
		int main_system_num_;

		//additional quantities
		UVector wear_;
		double wear_coefficient_;
		std::vector<double> total_wear_;


		static void update_wear(
			const double dt,
			const double wear_coefficient,
			const UVector &sliding_distance,
			const UVector &normal_stress,
			UVector &wear
			);
	};
}

#endif //UTOPIA_CONTACT_SYSTEM_HPP
