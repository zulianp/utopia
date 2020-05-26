#include "utopia_CSV.hpp"

#include <fstream>
#include <iostream>
#include <string>

namespace utopia {
    bool CSV::read(const Path &path) {
        std::ifstream is(path.c_str());

        if (!is.good()) {
            is.close();
            return false;
        }

        bool ok = read(is);
        is.close();
        return ok;
    }

    bool CSV::read(std::istream &is) {
        std::size_t cols = 0;
        while (is.good()) {
            std::string line;
            std::getline(is, line, '\n');

            if (line.empty()) continue;

            std::stringstream ss(line);

            table_.emplace_back();

            auto &row = table_.back();

            if (cols > 0) {
                row.reserve(cols);
            }

            while (ss.good()) {
                std::string value;
                std::getline(ss, value, separator_);
                row.push_back(value);
            }

            if (cols == 0) {
                cols = row.size();
            }
        }

        return true;
    }

    bool CSV::write(std::ostream &os) const {
        for (const auto &row : table_) {
            auto n = row.size();
            for (std::size_t i = 0; i < n; ++i) {
                os << row[i];

                if (i + 1 < n) {
                    os << separator_;
                }
            }

            os << "\n";
        }

        return true;
    }

    bool CSV::write(const Path &path) const {
        std::ofstream os(path.c_str());

        if (!os.good()) {
            os.close();
            return false;
        }

        bool ok = write(os);
        os.close();
        return ok;
    }

    const std::vector<std::string> &CSV::get_row(const SizeType i) const {
        static const std::vector<std::string> null;

        SizeType rows = table_.size();
        if (i < rows) {
            return table_[i];
        } else {
            return null;
        }
    }

    void CSV::get(const SizeType i, const SizeType j, double &val) const {
        std::string str;
        get(i, j, str);
        val = atof(str.c_str());
    }

    void CSV::get(const SizeType i, const SizeType j, int &val) const {
        std::string str;
        get(i, j, str);
        val = atoi(str.c_str());
    }

    void CSV::get(const SizeType i, const SizeType j, std::string &val) const {
        SizeType rows = table_.size();
        if (i < rows) {
            SizeType cols = table_[i].size();
            if (j < cols) {
                val = table_[i][j];
                return;
            }
        }

        std::cerr << "[Error] no value at " << i << ", " << j << std::endl;
    }

    void CSV::clear() { table_.clear(); }

}  // namespace utopia
