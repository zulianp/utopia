#ifndef UTOPIA_OBJECT_FACTORY_HPP
#define UTOPIA_OBJECT_FACTORY_HPP

#include <string>
#include <map>
#include "utopia_Base.hpp"
#include "utopia_FactoryMethod.hpp"
#include "utopia_AbstractVector.hpp"
#include "utopia_AbstractMatrix.hpp"
#include "utopia_Traits.hpp"
#include "utopia_Input.hpp"

namespace utopia {
    template<typename Scalar, typename SizeType>
    class AlgebraFactory : public Configurable {
    public:
        using BackendType    = std::string;
        using AbstractVector = utopia::AbstractVector<Scalar, SizeType>;
        using AbstractMatrix = utopia::AbstractMatrix<Scalar, SizeType>;

        inline void read(Input &in) override
        {
            in.get("default-backend", default_backend_);
        }

        static void init(Input &in)
        {
            instance().read(in);
        }

        static AlgebraFactory &instance()
        {
            static AlgebraFactory instance_;
            return instance_;
        }

        static void set_default_backend(const BackendType &backend)
        {
            instance().default_backend_ = backend;
        }

        static void print_info()
        {
            const auto &self = instance();

            std::cout << "Available implementations:\n";
            for(const auto &v : self.vector_factory_) {
                std::cout << "  - " << v.first << std::endl;
            }
        }

        //////////////////////////////////////////////////////////////////////////////////////////////////

        template<class Vector>
        static int register_vector()
        {
            static_assert(std::is_same<Scalar,  typename Traits<Vector>::Scalar>::value, "registered vector and factory scalar types have to be the same");
            static_assert(std::is_same<SizeType, typename Traits<Vector>::SizeType>::value, "registered vector and factory scalar types have to be the same");

            return instance().template register_vector_impl<Vector>();
        }

        inline static std::unique_ptr<AbstractVector> new_vector()
        {
            return instance().new_vector_impl(instance().default_backend_);
        }

        inline static std::unique_ptr<AbstractVector> new_vector(const BackendType &backend)
        {
            return instance().new_vector_impl(backend);
        }


        //////////////////////////////////////////////////////////////////////////////////////////////////

        template<class Matrix>
        static int register_matrix()
        {
            static_assert(std::is_same<Scalar,  typename Traits<Matrix>::Scalar>::value, "registered matrix and factory scalar types have to be the same");
            static_assert(std::is_same<SizeType, typename Traits<Matrix>::SizeType>::value, "registered matrix and factory scalar types have to be the same");

            return instance().template register_matrix_impl<Matrix>();
        }

        inline static std::unique_ptr<AbstractMatrix> new_matrix()
        {
            return instance().new_matrix_impl(instance().default_backend_);
        }

        inline static std::unique_ptr<AbstractMatrix> new_matrix(const BackendType &backend)
        {
            return instance().new_matrix_impl(backend);
        }

        //////////////////////////////////////////////////////////////////////////////////////////////////

    private:
        BackendType default_backend_;

        std::map<
            BackendType,
            std::unique_ptr< IFactoryMethod<AbstractVector> >
        > vector_factory_;

        std::map<
            BackendType,
            std::unique_ptr< IFactoryMethod<AbstractMatrix> >
        > matrix_factory_;


          AlgebraFactory()
        : default_backend_("")
        {}

        ////////////////////////////////////////////////////////////////////////////////

        template<class Vector>
        int register_vector_impl()
        {
            vector_factory_[Traits<Vector>::backend_info().get_name()] =
                utopia::make_unique< FactoryMethod<AbstractVector, Wrapper<Vector>> >();

            return Traits<Vector>::Backend;
        }

        std::unique_ptr<AbstractVector> new_vector_impl(const BackendType &backend)
        {
            if(vector_factory_.empty()) {
                std::cerr << "AlgebraFactory::new_vector(): no backend have registered to factory" << std::endl;
                return nullptr;
            }

            if(backend.empty()) {
                return vector_factory_.begin()->second->make();
            }

            auto it = vector_factory_.find(backend);
            if(it == vector_factory_.end()) {
                std::cerr << "AlgebraFactory::new_vector(): no backend with id " << backend << " has registered to factory" << std::endl;
                return nullptr;
            } else {
                return it->second->make();
            }
        }

        ////////////////////////////////////////////////////////////////////////////////

        template<class Matrix>
        int register_matrix_impl()
        {
            matrix_factory_[Traits<Matrix>::backend_info().get_name()] =
                utopia::make_unique< FactoryMethod<AbstractMatrix, Wrapper<Matrix>> >();

            return Traits<Matrix>::Backend;
        }

        std::unique_ptr<AbstractMatrix> new_matrix_impl(const BackendType &backend)
        {
            if(matrix_factory_.empty()) {
                std::cerr << "AlgebraFactory::new_matrix(): no backend have registered to factory" << std::endl;
                return nullptr;
            }

            if(backend.empty()) {
                return matrix_factory_.begin()->second->make();
            }

            auto it = matrix_factory_.find(backend);
            if(it == matrix_factory_.end()) {
                std::cerr << "AlgebraFactory::new_matrix(): no backend with id " << backend << " has registered to factory" << std::endl;
                return nullptr;
            } else {
                return it->second->make();
            }
        }
    };
}


#define UTOPIA_DEFINE_FACTORY_VAR(macro_in) dummy_factory_variable_ ## macro_in ## name
#define UTOPIA_FACTORY_REGISTER_VECTOR(macro_Vector_) \
 static int UTOPIA_DEFINE_FACTORY_VAR(macro_Vector_) = \
utopia::AlgebraFactory<typename utopia::Traits<macro_Vector_>::Scalar, typename utopia::Traits<macro_Vector_>::SizeType>::register_vector<macro_Vector_>();

#endif //UTOPIA_OBJECT_FACTORY_HPP
