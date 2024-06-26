
#include "utopia_Path.hpp"

#include <algorithm>
#include <fstream>
#include <iostream>
#include <string>
#include <utility>
#include <vector>
#include <filesystem>


// for windows looks for windows/dirent.h
#include <dirent.h>
#include <sys/stat.h>

#ifdef WIN32
static const char PATH_SEPARATOR = '\\';
#include <filesystem>
#else
#include <unistd.h>
static const char PATH_SEPARATOR = '/';
#endif

namespace utopia {

    Path::Path(std::string path) : path_(std::move(path)) { resolve_path_separators(path_); }

    Path::Path(const char *path) : path_(path) { resolve_path_separators(path_); }

    Path::~Path() = default;

    void Path::resolve_path_separators(std::string &path) {
#ifdef WIN32
        static const char toReplace = '/';
        static const char replacement = '\\';
#else
        static const char toReplace = '\\';
        static const char replacement = '/';
#endif
        std::replace(path.begin(), path.end(), toReplace, replacement);

        if (!path.empty() && path[path.size() - 1] == PATH_SEPARATOR) path.resize(path.size() - 1);
    }

    PathIterator Path::iter() const { return PathIterator(*this); }

    Path Path::without_extension() const {
        size_t index = path_.find_last_of(".");
        if (index >= path_.size()) {
            return *this;
        }

        return Path(path_.substr(0, index));
    }

    std::string Path::extension() const {
        std::size_t index = path_.find_last_of('.');

        if (std::string::npos == index) {
            return "";
        } else {
            return path_.substr(index + 1);
        }
    }

    Path Path::operator+(const Path &other) const {
        if (other.empty()) return *this;
        if (empty()) return other;

        return Path(path_ + other.path_);
    }

    Path Path::operator/(const Path &other) const {
        if (other.empty()) return *this;
        if (empty()) return other;

        std::string tmp = other.path_;
        if (tmp[0] != PATH_SEPARATOR) tmp = PATH_SEPARATOR + tmp;

        return Path(path_ + tmp);
    }

    std::string Path::file_name() const {
        const std::size_t found = path_.find_last_of(PATH_SEPARATOR);
        std::string name = path_.substr(found + 1, path_.size());
        const std::size_t dot = name.find_last_of(".");
        return name.substr(0, dot);
    }

    Path Path::parent() const {
        const size_t found = path_.find_last_of(PATH_SEPARATOR);
        if (found == std::string::npos) return Path(".");
        return path_.substr(0, found);
    }

    Path &Path::operator+=(const Path &other) {
        (*this) = (*this) + other;
        return *this;
    }

    Path &Path::operator/=(const Path &other) {
        (*this) = (*this) / other;
        return *this;
    }

    bool Path::is_dir() const {
        auto dir = opendir(path_.c_str());
        if (dir) {
            closedir((DIR *)dir);
            return true;
        } else {
            return false;
        }
    }

    bool Path::exists() const
    {
        const std::filesystem::path p{path_};
        return std::filesystem::exists(p);
    }

    bool Path::make_dir(const int permissions) const {
#ifdef WIN32
        auto ok = std::filesystem::create_directory(path_.c_str());
        if(!ok) return false;

        std::filesystem::permissions(path_.c_str(),
                          std::filesystem::perms(permissions));
#else
        int result = mkdir(path_.c_str(), permissions);

        return result == 0;
#endif
    }

    PathIterator::DirHandle::DirHandle(const std::string &path) { dir = opendir(path.c_str()); }

    PathIterator::DirHandle::~DirHandle() {
        if (dir) {
            closedir((DIR *)dir);
        }
    }

    struct dirent *PathIterator::DirHandle::next() {
        if (dir) return readdir((DIR *)dir);
        return nullptr;
    }

    PathIterator::PathIterator(const Path &path)
        : path_(path), dir_(new DirHandle(path.to_string())), it_(nullptr), skip_hidden_(true) {
        ++(*this);  // init and skip hidden
    }

    PathIterator::operator bool() const { return it_ != nullptr; }

    PathIterator &PathIterator::operator++() {
        it_ = dir_->next();
        if (!it_) return *this;

        std::string name = it_->d_name;
        while (skip_hidden_ && name[0] == '.') {
            it_ = dir_->next();
            if (!it_) break;
            name = it_->d_name;
        }
        return *this;
    }

    Path PathIterator::operator*() const { return Path(path_.to_string() + PATH_SEPARATOR + it_->d_name); }
}  // namespace utopia
