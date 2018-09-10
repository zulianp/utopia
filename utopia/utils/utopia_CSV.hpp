#ifndef UTOPIA_CSV_HPP
#define UTOPIA_CSV_HPP

#include "utopia_Base.hpp"
#include "utopia_Path.hpp"

#include <ostream>
#include <istream>
#include <vector>
#include <string>

namespace utopia {

	class CSV {
	public:

		CSV(const char separator = ',')
		: separator_(separator)
		{}

		bool read(const Path &path);
		bool read(std::istream &is);
		bool write(std::ostream &os) const;
		bool write(const Path &path) const;

		const std::vector<std::string> & get_row(const SizeType i) const;
		void get(const SizeType i, const SizeType j, double &val) const;
		void get(const SizeType i, const SizeType j, int &val) const;
		void get(const SizeType i, const SizeType j, std::string &val) const;


		void clear();
	private:
		char separator_;
		std::vector<std::vector<std::string>> table_;
	};
}

#endif //UTOPIA_CSV_HPP
