#ifndef UTOPIA_UTILS_HPP
#define UTOPIA_UTILS_HPP

#include <sys/stat.h>
#include <unistd.h>
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <memory>
#include <string>
#include <vector>

#include "utopia_MPI.hpp"
#include "utopia_Chrono.hpp"
#include "utopia_Path.hpp"
#include "utopia_CSV.hpp"

#ifdef UTOPIA_WITH_BRANCH_PREDICTION
#define UTOPIA_LIKELY(x)       __builtin_expect((x),1)
#define UTOPIA_UNLIKELY(x)     __builtin_expect((x),0)
#else
#define UTOPIA_LIKELY(x)     (x)
#define UTOPIA_UNLIKELY(x)   (x)
#endif

namespace utopia
{
    bool is_matlab_file(const std::string &path);


    inline std::string str(const char *char_array)
    {
        return std::string(char_array);
    }

    template<typename T>
    inline std::string str(const T &val)
    {
        return std::to_string(val);
    }
    template<typename T>
    inline bool is_nan_or_inf(const T &v)
    {
        return std::isinf(v) || std::isnan(v);
    }

    bool read(const std::string &path, std::string &str);
    bool write(const std::string &path, const std::string &str);
    void disp(const std::string &str);
    void disp(const char * str);

    inline void disp(const double value, std::ostream &os = std::cout)
    {
        os << value << "\n";
    }

    inline void disp(const int value, std::ostream &os = std::cout)
    {
        os << value << "\n";
    }

    inline void disp(const long value, std::ostream &os = std::cout)
    {
        os << value << "\n";
    }

    inline void disp(const long long int value, std::ostream &os = std::cout)
    {
        os << value << "\n";
    }

    template<typename T>
    void disp(const T & val, const std::string & name)
    {
        if(mpi_world_rank() == 0)
            std::cout<<"---  "<< name <<"  --- \n";
        return disp(val);
    }

    std::vector<std::string> split_string(const std::string &, char);

    /** \addtogroup Utilities
     * @ingroup ptr_utilities
     * @brief Provide basic functionalities for pointer manipulations.
     *  @{
     */

    template <typename T>
    struct EmptyDeleter
    {
        EmptyDeleter() /* noexcept */
            = default;

        template <typename U>
        EmptyDeleter(const EmptyDeleter<U>&,
            typename std::enable_if<
                std::is_convertible<U*, T*>::value
            >::type* = nullptr) /* noexcept */
        {
        }

        void operator()(T* const) const /* noexcept */
        {
            // do nothing
        }
    };


    /*!
     * @fn std::shared_ptr<T> make_ref(T &obj)
     * @brief use only when you are sure that the scope of the pointer is limited to the life of obj
     * @return a shared_ptr that does not delete the pointer
     */
    template<typename T>
    std::shared_ptr<T> make_ref(T &obj)
    {
        return std::shared_ptr<T>(&obj, EmptyDeleter<T>());
    }

    /*!
     * @fn std::shared_ptr<const T> make_ref(const T &obj)
     * @brief use only when you are sure that the scope of the pointer is limited to the life of obj
     * @return a shared_ptr that does not delete the pointer
     */
    template<typename T>
    std::shared_ptr<const T> make_ref(const T &obj)
    {
        return std::shared_ptr<const T>(&obj, EmptyDeleter<const T>());
    }



    enum ColorCode {
        FG_RED      		= 31,
        FG_GREEN    		= 32,
        FG_YELLOW 			= 33,
        FG_BLUE     		= 34,
        FG_MAGENTA 			= 35,
        FG_CYAN 			= 36,
        FG_LIGHT_GRAY 		= 37,
        FG_DEFAULT  		= 39,
        FG_DARK_GRAY 		= 90,
        FG_LIGHT_RED 		= 91,
        FG_LIGHT_GREEN 		= 92,
        FG_LIGHT_YELLOW		= 93,
        FG_LIGHT_BLUE 		= 94,
        FG_LIGHT_MAGENTA	= 95,
        FG_LIGHT_CYAN 		= 96,
        FG_WHITE	 		= 97,

        BG_RED      = 41,
        BG_GREEN    = 42,
        BG_BLUE     = 44,
        BG_DEFAULT  = 49
    };

    class ColorModifier
    {
        ColorCode code;

    public:

        void set_color_code(ColorCode new_code){ code = new_code;  }

        ColorModifier(ColorCode pCode) : code(pCode) {}
        friend std::ostream&
        operator<<(std::ostream& os, const ColorModifier& mod) {
            return os << "\033[" << mod.code << "m";
        }
    };


    class CSVWriter
    {
         public:
             inline constexpr static const char * append_binary() { return "ab+"; }
             inline constexpr static const char * read() { return "r"; }
             inline constexpr static const char * write() { return "w"; }
             inline constexpr static const char * append() { return "a+"; }

             void open_file(const std::string & file_name, const std::string &mode = append_binary())
             {
                 pFile = fopen(file_name.c_str(), mode.c_str());
             }

             template<class T>
             void write_table_row(const std::vector<T> vars);

             void close_file()
             {
                 if (pFile != nullptr && mpi_world_rank() == 0) fclose(pFile);
             }


             inline bool file_exists(const std::string& file_name)
             {
                 struct stat buffer {};
                 return (stat(file_name.c_str(), &buffer) == 0);
            }


         private:
             FILE *pFile;

    };

    template<typename data_type>
    bool in_array(const std::string &value, const std::vector<data_type> &array)
    {
        return std::find(array.begin(), array.end(), value) != array.end();
    }

      /** @}*/
}

#endif //UTOPIA_UTILS_HPP
