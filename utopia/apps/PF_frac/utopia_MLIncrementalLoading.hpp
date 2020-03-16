#ifndef UTOPIA_DM_RMTR_SETUP_HPP
#define UTOPIA_DM_RMTR_SETUP_HPP

#include "utopia_Input.hpp"
#include "utopia_Multigrid.hpp"
#include "utopia_IPTransfer.hpp"
#include "utopia_make_unique.hpp"
#include "utopia_IsotropicPhaseField.hpp"

#include <memory>

namespace utopia {

    //FIXME complete the overriding process
    template<class FunctionSpace, class ProblemType, class BCType, class ICType>
    class MLIncrementalLoading final : public Configurable {
    public:
        using Matrix   = typename FunctionSpace::Matrix;
        using Vector   = typename FunctionSpace::Vector;
        using Scalar   = typename FunctionSpace::Scalar;
        using SizeType = typename FunctionSpace::SizeType;


        MLIncrementalLoading(FunctionSpace &space_coarse, const SizeType & n_levels) : n_levels_(n_levels){
            init(space_coarse); 
        }   

        void read(Input &in) override {

            for (auto l=0; l < level_functions_.size(); l++){
                level_functions_[l]->read(in);    
            }

            
        }


        bool init(FunctionSpace &space){
            return init(make_ref(space));
        }

        bool init(const std::shared_ptr<FunctionSpace> &space)
        {
            if(n_levels_ < 2) {
                std::cerr << "n_levels must be at least 2" << std::endl;
                return false;
            }

            spaces_.resize(n_levels_);
            spaces_[0] = space;

            level_functions_.resize(n_levels_); 

            auto fun = std::make_shared<ProblemType>(*spaces_[0]);
            level_functions_[0] = fun; 


            transfers_.resize(n_levels_ - 1);

            for(SizeType i = 1; i < n_levels_; ++i) {
                spaces_[i] = spaces_[i-1]->uniform_refine();


                auto fun = std::make_shared<ProblemType>(*spaces_[i]);
                level_functions_[i] = fun; 


                auto I = std::make_shared<Matrix>();
                spaces_[i-1]->create_interpolation(*spaces_[i], *I);
                assert(!empty(*I));
                transfers_[i-1] = std::make_shared<IPTransfer<Matrix, Vector> >(I);
            }

            return true;
        }

        FunctionSpace &fine_space()
        {
            return *spaces_.back();
        }

        const FunctionSpace &fine_space() const
        {
            return *spaces_.back();
        }

        std::shared_ptr<FunctionSpace> fine_space_ptr()
        {
            return spaces_.back();
        }

        std::shared_ptr<const FunctionSpace> fine_space_ptr() const
        {
            return spaces_.back();
        }

    private:
        SizeType n_levels_; 

        std::vector<std::shared_ptr<FunctionSpace>> spaces_;
        std::vector<std::shared_ptr<Transfer<Matrix, Vector> > > transfers_;
        std::vector<std::shared_ptr<ExtendedFunction<Matrix, Vector> > >  level_functions_;

    };

}

#endif //UTOPIA_DM_RMTR_SETUP_HPP
