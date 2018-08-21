
#include "utopia_Path.hpp"

#include <algorithm>
#include <vector>
#include <string>
#include <iostream>
#include <fstream>

//for windows looks for windows/dirent.h
#include <dirent.h>

#ifdef WIN32
static const char PATH_SEPARATOR='\\';
#else
static const char PATH_SEPARATOR='/';
#endif

namespace utopia {

	Path::Path(const std::string &path)
		: path_(path)
	{
		resolve_path_separators(path_);
	}

	Path::Path(const char *path)
		: path_(path)
	{
		resolve_path_separators(path_);
	}

	Path::~Path() {}

	void Path::resolve_path_separators(std::string &path)
	{
#ifdef WIN32
		static const char toReplace = '/';
		static const char replacement = '\\';
#else
		static const char toReplace = '\\';
		static const char replacement = '/';
#endif
		std::replace(path.begin(), path.end(), toReplace, replacement);

		if(!path.empty() && path[path.size()-1]==PATH_SEPARATOR)
			path.resize(path.size()-1);
	}

	PathIterator Path::iter() const { return PathIterator(*this); }

	Path Path::without_extension() const
	{
		size_t index = path_.find_last_of(".");
		if(index >= path_.size()) {
			return *this;
		}

		return Path(path_.substr(0, index));
	}

	std::string Path::extension() const
	{
		std::size_t index = path_.find_last_of('.');
	
		if(std::string::npos == index) {
			return "";
		} else {
			return path_.substr(index + 1);
		}

	}

	Path Path::operator+(const Path &other) const
	{
		if(other.empty())
			return *this;
		if(empty()) 
			return other;

		return Path(path_ + other.path_);
	} 

	Path Path::operator/(const Path &other) const
	{
		if(other.empty())
			return *this;
		if(empty()) 
			return other;


		std::string tmp=other.path_;
		if(tmp[0]!=PATH_SEPARATOR)
			tmp=PATH_SEPARATOR+tmp;

		return Path(path_ + tmp);
	} 

	std::string Path::file_name() const
	{
		const std::size_t found = path_.find_last_of(PATH_SEPARATOR);
		std::string name=path_.substr(found+1,path_.size());
		const std::size_t dot = name.find_last_of(".");
		return name.substr(0,dot);
	}

	Path Path::parent() const 
	{
		const size_t found = path_.find_last_of(PATH_SEPARATOR);
		if(found == std::string::npos) return Path(".");
		return path_.substr(0, found);
	}

	Path & Path::operator+=(const Path &other)
	{
		(*this)=(*this)+other;
		return *this;
	}

	Path & Path::operator/=(const Path &other)
	{
		(*this)=(*this)/other;
		return *this;
	}

	PathIterator::DirHandle::DirHandle(const std::string &path)
	{
		dir = opendir(path.c_str());
	}	

	PathIterator::DirHandle::~DirHandle()
	{
		if(dir)
			closedir((DIR *)dir);
	}

	struct dirent * PathIterator::DirHandle::next() 
	{ 
		if(dir)
			return readdir((DIR *)dir); 
		return NULL;
	}


	PathIterator::PathIterator(const Path &path)
		: path_(path), dir_(new DirHandle(path.to_string())), it_(NULL), skip_hidden_(true)
	{
		++(*this); //init and skip hidden
	}

	PathIterator::operator bool() const 
	{
		return it_ != NULL;
	} 

	PathIterator & PathIterator::operator ++()
	{
		it_ = dir_->next();
		if(!it_)
			return *this;


		std::string name = it_->d_name;
		while(skip_hidden_ && name[0] == '.') {
			it_ = dir_->next();
			if(!it_)
				break;
			name = it_->d_name;
		}
		return *this;
	}

	Path PathIterator::operator *() const
	{
		return Path(path_.to_string() + PATH_SEPARATOR + it_->d_name);
	}
}
