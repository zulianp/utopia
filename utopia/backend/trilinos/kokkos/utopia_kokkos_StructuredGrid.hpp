#ifndef UTOPIA_KOKKOS_StructuredGrid_hpp
#define UTOPIA_KOKKOS_StructuredGrid_hpp

#include "utopia_Base.hpp"
#include "utopia_Core.hpp"
#include "utopia_kokkos_Traits.hpp"
#include "utopia_Views.hpp"
//"



namespace utopia {

  template<int Dim,typename Vector,int Backend = Traits<Vector>::Backend>
  class StructuredGrid {};

  template<typename Vector>
  class StructuredGrid<2,Vector,TRILINOS>  {
    public:

        using Scalar   = typename Traits<Vector>::Scalar;
        using SizeType = typename Traits<Vector>::SizeType;

        //FIXME
        typedef Kokkos::TeamPolicy<>               TeamPolicy;
        typedef Kokkos::TeamPolicy<>::member_type  MemberType;


        using DeviceVector      = utopia::VectorView<Kokkos::View<Scalar *>>;

        //FIXME

        using DofView           = Kokkos::View<SizeType **>;
        using PointView         = Kokkos::View<Scalar **>;
        

        template<typename... Args>
        inline static void parallel_for(Args&&... args)
        {
            using ForLoop  = utopia::ParallelFor<Traits<Vector>::Backend>;
            return ForLoop::apply(std::forward<Args>(args)...);
        }

        // using VectorD = utopia::VectorView<Kokkos::View<Scalar[Dim]>>;

        static SizeType compute_n_elements(const SizeType &n)
        {
            SizeType ret = n;
            for(int i = 1; i < 2; ++i) {
                ret *= n;
            }

            return ret;
        }

        static SizeType compute_n_points(const SizeType &n)
        {
            SizeType ret = n + 1;
            for(int i = 1; i < 2; ++i) {
                ret *= (n + 1);
            }

            return ret;
        }






      StructuredGrid(const SizeType &n) :
        n_(n),
        n_elements_(compute_n_elements(n)),
        n_points_(compute_n_points(n)),
        dof_("dof", n_elements_, n),
        point_("point", n_points_, 2)
        {

          init_mesh();
        }

        void init_mesh()
        {
   

    
            const Scalar h = 1./n_;
            parallel_for(
                n_,
                UTOPIA_LAMBDA(const SizeType &i)
                {
                    // Kokkos::parallel_for(Kokkos::TeamThreadRange(team_member, n_), [&] (const SizeType j) {
                    for(SizeType j = 0; j < n_; ++j) {

                        const SizeType e_id = i * n_ + j;
                        DeviceVector p(Kokkos::subview(point_, e_id, Kokkos::ALL()));

                        p.set(0, i * h);
                        p.set(1, j * h);

                        const SizeType n_p = n_ + 1;

                        dof_(e_id, 0) = i * n_p + j;
                        dof_(e_id, 1) = i * n_p + (j + 1);
                        //flipped for dof-consistency and ccw local
                        dof_(e_id, 2) = (i + 1) * n_p + (j + 1);
                        dof_(e_id, 3) = (i + 1) * n_p + j;
                      }
                }
            );
        }

        void dof_indices(const SizeType &e_id, DofView &dof_indices)
        {
           int size = 4;

                    // Kokkos::parallel_for(Kokkos::TeamThreadRange(team_member, n_), [&] (const SizeType j) {
            for(SizeType j = 0; j < size; ++j){

                dof_indices[j] = dof_(e_id, j);

            }
           
        }


      private:
        SizeType n_;
        SizeType n_elements_;
        SizeType n_points_;
        DofView dof_;
        PointView point_;

      

  };


  

}

#endif

