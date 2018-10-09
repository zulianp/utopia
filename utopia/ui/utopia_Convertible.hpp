#ifndef UTOPIA_CONVERTIBLE_HPP
#define UTOPIA_CONVERTIBLE_HPP

#include <string>

namespace utopia {
	class IConvertible {
	public:
		virtual ~IConvertible() {}
		virtual void get(double &) const = 0;
		virtual void get(float &) const = 0;
		virtual void get(int &) const = 0;
		virtual void get(long &) const = 0;
		virtual void get(bool &) const = 0;
		virtual void get(std::string &) const = 0;
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
			out = in == "1" || in == "true" || in == "yes";
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

	template<typename T>
	class Convertible final : public IConvertible {
	public:
	
		void get(double &in_out) const override
		{
			Convert<T, double>::apply(value_, in_out);
		}

		void get(float &in_out) const override
		{
			Convert<T, float>::apply(value_, in_out);
		}

		void get(int &in_out) const override
		{
			Convert<T, int>::apply(value_, in_out);
		}

		void get(long &in_out) const override
		{
			Convert<T, long>::apply(value_, in_out);
		}

		void get(bool &in_out) const override
		{
			Convert<T, bool>::apply(value_, in_out);
		}

		void get(std::string &in_out) const override
		{
			Convert<T, std::string>::apply(value_, in_out);
		}

	private:
		T value_;
	};
}


#endif //UTOPIA_CONVERTIBLE_HPP
