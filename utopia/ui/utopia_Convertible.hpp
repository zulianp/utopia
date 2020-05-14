#ifndef UTOPIA_CONVERTIBLE_HPP
#define UTOPIA_CONVERTIBLE_HPP

#include <string>
#include <type_traits>

#include "utopia_Clonable.hpp"

namespace utopia {
    class IConvertible : public Clonable {
    public:
        ~IConvertible() override {}
        virtual void get(double &) const = 0;
        virtual void get(float &) const = 0;
        virtual void get(int &) const = 0;
        virtual void get(long &) const = 0;
        // virtual void get(long long &) const = 0;
        virtual void get(unsigned long &) const = 0;
        virtual void get(long long int &) const = 0;
        virtual void get(bool &) const = 0;
        virtual void get(std::string &) const = 0;

        virtual bool is_double() const { return false; }
        virtual bool is_float() const { return false; }
        virtual bool is_int() const { return false; }
        virtual bool is_long() const { return false; }
        virtual bool is_longlong() const { return false; }
        virtual bool is_ulong() const { return false; }
        virtual bool is_bool() const { return false; }
        virtual bool is_string() const { return false; }
        virtual bool is_long_long_int() const { return false; }
        IConvertible *clone() const override = 0;
    };

    //does not do anything (no defaults)
    class NullConvertible final : public IConvertible {
    public:
        void get(double &) const override {}
        void get(float &) const override {}
        void get(int &) const override {}
        void get(long &) const override {}
        // void get(long long &) const override {}
        void get(unsigned long &) const override {}
        void get(long long int &) const override {}
        void get(bool &) const override {}
        void get(std::string &) const override {}

        NullConvertible * clone() const override
        {
            return new NullConvertible();
        }
    };

    template<typename In, typename Out>
    class Convert {
    public:
        static void apply(const In &in, Out &out)
        {
            out = static_cast<Out>(in);
        }
    };

    template<typename InOut>
    class Convert<InOut, InOut> {
    public:
        static void apply(const InOut &in, InOut &out)
        {
            out = in;
        }
    };

    template<typename In>
    class Convert<In, std::string> {
    public:
        static void apply(const In &in, std::string &out)
        {
            out = std::to_string(in);
        }
    };

    template<>
    class Convert<std::string, std::string> {
    public:
        static void apply(const std::string &in, std::string &out)
        {
            out = in;
        }
    };

    template<>
    class Convert<std::string, bool> {
    public:
        static void apply(const std::string &in, bool &out)
        {
            out = (in == "1" || in == "true" || in == "yes");
        }
    };

    template<>
    class Convert<std::string, double> {
    public:
        static void apply(const std::string &in, double &out)
        {
            out = atof(in.c_str());
        }
    };

    template<>
    class Convert<std::string, float> {
    public:
        static void apply(const std::string &in, float &out)
        {
            out = atof(in.c_str());
        }
    };

    template<>
    class Convert<std::string, int> {
    public:
        static void apply(const std::string &in, int &out)
        {
            out = atoi(in.c_str());
        }
    };

    template<>
    class Convert<std::string, long> {
    public:
        static void apply(const std::string &in, long &out)
        {
            out = atol(in.c_str());
        }
    };

    // template<>
    // class Convert<std::string, long long> {
    // public:
    //     static void apply(const std::string &in, long long &out)
    //     {
    //         out = atol(in.c_str());
    //     }
    // };

    template<>
    class Convert<std::string, unsigned long> {
    public:
        static void apply(const std::string &in, unsigned long &out)
        {
            //FIXME
            out = atol(in.c_str());
        }
    };

    template<>
    class Convert<std::string, long long int> {
    public:
        static void apply(const std::string &in, long long int &out)
        {
            //FIXME
            out = atoll(in.c_str());
        }
    };

    template<typename T>
    class Convertible final : public IConvertible {
    public:
        Convertible(const T &value) : value_(value) {}

        inline void get(double &in_out) const override
        {
            Convert<T, double>::apply(value_, in_out);
        }

        inline void get(float &in_out) const override
        {
            Convert<T, float>::apply(value_, in_out);
        }

        inline void get(int &in_out) const override
        {
            Convert<T, int>::apply(value_, in_out);
        }

        inline void get(long &in_out) const override
        {
            Convert<T, long>::apply(value_, in_out);
        }

        // inline void get(long long &in_out) const override
        // {
        //     Convert<T, long long>::apply(value_, in_out);
        // }

        inline void get(unsigned long &in_out) const override
        {
            Convert<T, unsigned long>::apply(value_, in_out);
        }

        inline void get(long long int &in_out) const override
        {
            Convert<T, long long int>::apply(value_, in_out);
        }

        inline void get(bool &in_out) const override
        {
            Convert<T, bool>::apply(value_, in_out);
        }

        inline void get(std::string &in_out) const override
        {
            Convert<T, std::string>::apply(value_, in_out);
        }

        ////////////////////////////////////////

        inline void set(const double &in)
        {
            Convert<double, T>::apply(in, value_);
        }

        inline void set(const float &in)
        {
            Convert<float, T>::apply(in, value_);
        }

        inline void set(const int &in)
        {
            Convert<int, T>::apply(in, value_);
        }

        inline void set(const long &in)
        {
            Convert<long, T>::apply(in, value_);
        }

        // inline void set(const long long &in)
        // {
        //     Convert<long long, T>::apply(in, value_);
        // }

        inline void set(const unsigned long &in)
        {
            Convert<long, T>::apply(in, value_);
        }

        inline void set(const long long int &in)
        {
            Convert<long long int, T>::apply(in, value_);
        }

        inline void set(const bool &in)
        {
            Convert<bool, T>::apply(in, value_);
        }

        inline void set(const std::string &in)
        {
            Convert<std::string, T>::apply(in, value_);
        }

        inline operator const T &() const
        {
            return value_;
        }

        inline operator T &()
        {
            return value_;
        }

        inline bool is_double() const override  { return std::is_same<T, double>::value; }
        inline bool is_float() const override 	{ return std::is_same<T, float>::value; }
        inline bool is_int() const override 	{ return std::is_same<T, int>::value; }
        inline bool is_long() const override 	{ return std::is_same<T, long>::value; }
        inline bool is_ulong() const override 	{ return std::is_same<T, long>::value; }
        inline bool is_bool() const override 	{ return std::is_same<T, bool>::value; }
        inline bool is_string() const override  { return std::is_same<T, std::string>::value; }
        inline bool is_long_long_int() const override { return std::is_same<T, long long int>::value; }

        inline Convertible * clone() const override
        {
            return new Convertible(value_);
        }

    private:
        T value_;
    };

    using Double = Convertible<double>;
    using Float  = Convertible<float>;
    using Int    = Convertible<int>;
    using Long   = Convertible<long>;
    using LLong  = Convertible<long long int>;
    using ULong  = Convertible<unsigned long>;
    using Bool   = Convertible<bool>;
    using String = Convertible<std::string>;
}


#endif //UTOPIA_CONVERTIBLE_HPP
