
#ifndef UTOPIA_PATH_HPP
#define UTOPIA_PATH_HPP

#include <memory>
#include <ostream>
#include <string>

#include "utopia_Describable.hpp"

// forward declaration of dirent included in cpp
struct dirent;

namespace utopia {
    class Path;
    class PathIterator;

    class PathIterator {
    private:
        struct DirHandle {
            // DIR * hidden
            void *dir;
            DirHandle(const std::string &path);
            ~DirHandle();
            struct dirent *next();
        };

        const Path &path_;
        std::shared_ptr<DirHandle> dir_;
        struct dirent *it_;
        bool skip_hidden_;

    public:
        PathIterator(const Path &path);
        operator bool() const;
        PathIterator &operator++();
        Path operator*() const;
    };

    class Path : public Describable {
    private:
        std::string path_;
        static void resolve_path_separators(std::string &path);

    public:
        using Iterator = PathIterator;

        Path(std::string path = "");
        Path(const char *path);

        Path operator+(const Path &other) const;
        Path operator/(const Path &other) const;

        Path &operator+=(const Path &other);
        Path &operator/=(const Path &other);

        inline const std::string &to_string() const { return path_; }
        inline operator const std::string &() const { return path_; }
        inline const char *c_str() const { return path_.c_str(); }
        friend std::ostream &operator<<(std::ostream &os, const Path &path) {
            os << path.to_string();
            return os;
        }

        inline void describe(std::ostream &os) const override { os << this->to_string(); }

        inline bool empty() const { return path_.empty(); }
        std::string file_name() const;
        std::string extension() const;
        Path without_extension() const;
        Path parent() const;

        PathIterator iter() const;

        virtual ~Path();

        bool is_dir() const;
        bool exists() const;
        bool make_dir(const int permissions = 0777) const;

        inline std::string &raw_type() { return path_; }
    };

}  // namespace utopia

#endif  // UTOPIA_PATH_HPP
