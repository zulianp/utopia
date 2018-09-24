#ifndef UTOPIA_CSV_HPP
#define UTOPIA_CSV_HPP

#include "utopia_Base.hpp"
#include "utopia_Path.hpp"

#include <ostream>
#include <istream>
#include <vector>
#include <string>

namespace utopia {

	class CSV final {
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
		inline std::size_t n_rows() const
		{
			return table_.size();
		}

		inline std::size_t n_cols() const
		{
			if(table_.empty()) {
				return 0;
			}

			return table_[0].size();
		}


		void clear();
	private:
		char separator_;
		std::vector<std::vector<std::string>> table_;
	};
}

#endif //UTOPIA_CSV_HPP
