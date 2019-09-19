
#ifndef UTOPIA_TPETRA_VECTOR_HPP
#define UTOPIA_TPETRA_VECTOR_HPP

#include "utopia_Range.hpp"
#include "utopia_make_unique.hpp"
#include "utopia_Base.hpp"
#include "utopia_Size.hpp"
#include "utopia_Tensor.hpp"
#include "utopia_Writable.hpp"
#include "utopia_Constructible.hpp"
#include "utopia_Reducible.hpp"
#include "utopia_ElementWiseOperand.hpp"
#include "utopia_Transformable.hpp"
#include "utopia_Comparable.hpp"
#include "utopia_Vector.hpp"

#include "utopia_kokkos_Eval_Binary.hpp"
#include "utopia_kokkos_Eval_Unary.hpp"

#include <Tpetra_Map_decl.hpp>
#include <Tpetra_Vector_decl.hpp>

#include <Kokkos_DefaultNode.hpp>
#include <Tpetra_CrsMatrix.hpp>

#include "utopia_trilinos_Traits.hpp"
#include "utopia_trilinos_Communicator.hpp"

#include <memory>

namespace utopia {

    class TpetraVector : 
    public DistributedVector<TpetraScalar, TpetraSizeType>,
    public ElementWiseOperand<TpetraVector>,
    public ElementWiseOperand<TpetraScalar>,
    public Transformable<TpetraScalar>,
    public Tensor<TpetraVector, 1>,
    public Constructible<TpetraScalar, TpetraSizeType, 1>,
    public Reducible<TpetraScalar>,
    public Comparable<TpetraVector>
    {
    public:

        using SizeType = utopia::TpetraSizeType;
        using Scalar   = utopia::TpetraScalar;


    typedef Tpetra::Operator<>::scalar_type SC;
    typedef Tpetra::Operator<SC>::local_ordinal_type LO;
    typedef Tpetra::Operator<SC, LO>::global_ordinal_type GO;

    //types of Kokkos Parallel Nodes
    typedef Kokkos::Compat::KokkosSerialWrapperNode serial_node;

#ifdef KOKKOS_ENABLE_CUDA
    typedef Kokkos::Compat::KokkosCudaWrapperNode cuda_node;
    typedef cuda_node NT;
#elif defined KOKKOS_ENABLE_ROCM //Kokkos::Compat::KokkosROCmWrapperNode doesn't exist
    typedef Kokkos::Compat::KokkosDeviceWrapperNode<Kokkos::ROCm> rocm_node;
    typedef rocm_node NT;
#elif defined KOKKOS_ENABLE_OPENMP
    typedef Kokkos::Compat::KokkosOpenMPWrapperNode openmp_node;
    typedef openmp_node NT;
#else
    typedef serial_node NT;
#endif

    typedef Tpetra::Operator<SC, LO, GO, NT> OP;

    typedef Tpetra::Map<LO, GO, NT>                   map_type;
    typedef Tpetra::Vector<SC, LO, GO, NT>            vector_type;
    typedef Tpetra::MultiVector<SC, LO, GO, NT>       multi_vector_type;
    typedef Teuchos::RCP<vector_type>                 rcpvector_type;
    typedef Teuchos::RCP<const Teuchos::Comm<int> >   rcp_comm_type;
    typedef Teuchos::RCP<const map_type>              rcp_map_type;

//        typedef Tpetra::Vector<>::mag_type                magnitude_type;
    // typedef vector_type::scalar_type                  Scalar;

        using IndexSet = Traits<TpetraVector>::IndexSet;


        ////////////////////////////////////////////////////////////////////
        ///////////////////////// BOILERPLATE CODE FOR EDSL ////////////////
        ////////////////////////////////////////////////////////////////////

        using Super         = utopia::Tensor<TpetraVector, 1>;
        using Constructible = utopia::Constructible<Scalar, GO, 1>;

        using Super::Super;
        // using Constructible::values;
        // using Constructible::local_values;
        // using Constructible::local_zeros;
        // using Constructible::zeros;

         template<class Expr>
         TpetraVector(const Expression<Expr> &expr)
         {
             //THIS HAS TO BE HERE IN EVERY UTOPIA TENSOR CLASS
             Super::construct_eval(expr.derived());
         }

         template<class Expr>
         inline TpetraVector &operator=(const Expression<Expr> &expr)
         {
             Super::assign_eval(expr.derived());
             return *this;
         }

        void assign(const TpetraVector &other) override;
        void assign(TpetraVector &&other) override;


        ////////////////////////////////
        TpetraVector(const TrilinosCommunicator &comm = Tpetra::getDefaultComm())
        : comm_(comm)
        {}

        ~TpetraVector()
        {}

        TpetraVector(const TpetraVector &other);


        TpetraVector(TpetraVector &&other)
        : vec_(std::move(other.vec_))
        { }

        rcp_comm_type communicator()
        {
            return implementation().getMap()->getComm();
        }

        const rcp_comm_type communicator() const
        {
            return implementation().getMap()->getComm();
        }

        TpetraVector &operator=(const TpetraVector &other);

        TpetraVector &operator=(TpetraVector &&other)
        {
            if(this == &other) return *this;
            vec_ = std::move(other.vec_);
            ghosted_vec_ = std::move(other.ghosted_vec_);
            return *this;
        }

        void select(const IndexSet &index, TpetraVector &result) const override
        {
            generic_select(index, result);
        }

        void copy(const TpetraVector &other);

        /////////////////////////////////////////////////////////////////////////////////
        //////////////////////////////// OVERRIDES for Constructible ////////////////////
        /////////////////////////////////////////////////////////////////////////////////

        void values(const SizeType &s, const Scalar &val) override
        {
            assert(false && "IMPLEMENT");
        } 


        /////////////////////////////////////////////////////////////////////////////////
        //////////////////////////////// OVERRIDES for Reducible ////////////////////
        /////////////////////////////////////////////////////////////////////////////////

        inline Scalar reduce(const Plus &) const override
        {
            return sum();
        }

        inline Scalar reduce(const Min &)  const override
        {
            return min();
        }

        inline Scalar reduce(const Max &)  const override
        {
            return max();
        }

        Scalar sum() const override;
        Scalar min() const override;
        Scalar max() const override;

        /////////////////////////////////////////////////////////////////////////////////
        //////////////////////////////// OVERRIDES for ElementWiseOperand ////////////////////
        /////////////////////////////////////////////////////////////////////////////////

        void e_mul(const TpetraVector &other) override
        {
            raw_type()->elementWiseMultiply(
                1.,
                *raw_type(),
                *other.raw_type(),
                0.
            );
        }

        void e_div(const TpetraVector &other) override
        {
            KokkosEvalBinary<TpetraVector, Divides>::eval(*this, Divides(), other, *this);
        }

        void e_min(const TpetraVector &other) override
        {
            KokkosEvalBinary<TpetraVector, Min>::eval(*this, Min(), other, *this);
        }

        void e_max(const TpetraVector &other) override
        {
            KokkosEvalBinary<TpetraVector, Max>::eval(*this, Max(), other, *this);
        }

        ////////////////////////////////////////////////////////////////

        void e_mul(const Scalar &other) override
        {
            scale(other);
        }

        void e_div(const Scalar &other) override
        {
           assert(false && "IMPLEMENT ME");
        }

        void e_min(const Scalar &other) override
        {
            assert(false && "IMPLEMENT ME");
        }

        void e_max(const Scalar &other) override
        {
            assert(false && "IMPLEMENT ME");
        }

        /////////////////////////////////////////////////////////////////////////////////
        //////////////////////////////// OVERRIDES for Transformable ////////////////////
        /////////////////////////////////////////////////////////////////////////////////


        void transform(const Sqrt &) override;
        void transform(const Pow2 &) override;
        void transform(const Log &)  override;
        void transform(const Exp &)  override;
        void transform(const Cos &)  override;
        void transform(const Sin &)  override;
        void transform(const Abs &)  override;
        void transform(const Minus &) override;

        void transform(const Pow &p)  override;
        void transform(const Reciprocal<TpetraScalar> &f) override;

        /////////////////////////////////////////////////////////////////////////////////
        //////////////////////////////// OVERRIDES for Comparable ////////////////////
        /////////////////////////////////////////////////////////////////////////////////

        bool equals(const TpetraVector &other, const Scalar &tol = 0.0) const override;

        //////////////////////////////////////////
        //API functions
        //////////////////////////////////////////
        inline void values(const rcp_comm_type &comm, std::size_t n_local, Tpetra::global_size_t n_global, Scalar value)
        {
            rcp_map_type map;

            if(n_local == INVALID_INDEX) {
                map = Teuchos::rcp(new map_type(n_global, 0, comm));
            } else {
                map = Teuchos::rcp(new map_type(n_global, n_local, 0, comm));
            }

            vec_.reset(new vector_type(map));
            implementation().putScalar(value);
        }

        inline void init(const rcp_map_type &map)
        {
            vec_.reset(new vector_type(map));
        }

        void ghosted(const TpetraVector::GO &local_size,
                     const TpetraVector::GO &global_size,
                     const std::vector<GO> &ghost_index
        )
        {
            ghosted(comm().get(), local_size, global_size, ghost_index);
        }

        void ghosted(const rcp_comm_type &comm,
                     const TpetraVector::GO &local_size,
                     const TpetraVector::GO &global_size,
                     const std::vector<GO> &ghost_index
        );

        inline void axpy(const Scalar &alpha, const TpetraVector &x)
        {
            implementation().update(alpha, *x.vec_, 1.);
        }

        void describe(std::ostream &os) const
        {
            auto out = Teuchos::getFancyOStream(Teuchos::rcpFromRef(os));
            implementation().describe(*out, Teuchos::EVerbosityLevel::VERB_EXTREME);
        }


        void describe() const
        {
            describe(std::cout);
        }

        inline Scalar get(const SizeType &i) const override
        {
            // assert(!read_only_data_.is_null() && "Use Read<Vector> w(v); to enable reading from this vector v!");
            // return read_only_data_[local_index(i)];

            auto local_index = view_ptr_->map.getLocalElement(i);
            return view_ptr_->view(local_index, 0);
        }

        inline Scalar operator[](const GO i) const
        {
            // assert(!read_only_data_.is_null() && "Use Read<Vector> w(v); to enable reading from this vector v!");
            return get(i);
        }

        inline LO local_index(const GO i) const
        {
            if(has_ghosts()) {
               return ghosted_vec_->getMap()->getLocalElement(i);
            } else {
                //i - implementation().getMap()->getMinGlobalIndex()
                return implementation().getMap()->getLocalElement(i);
            }
        }

        inline void set(const SizeType &i, const Scalar &value) override
        {
            assert(view_ptr_);
            auto local_index = view_ptr_->map.getLocalElement(i);

            view_ptr_->view(local_index, 0) = value;
            // if(!write_data_.is_null()) {
            //     write_data_[i - implementation().getMap()->getMinGlobalIndex()] = value;
            // } else {
            //     implementation().replaceGlobalValue(i, value);
            // }
        }

        inline void add(const SizeType &i, const Scalar &value) override
        {
            assert(view_ptr_);

            auto local_index = view_ptr_->map.getLocalElement(i);

            view_ptr_->view(local_index, 0) += value;

            // if(!ghosted_vec_.is_null()) {
            //     ghosted_vec_->sumIntoGlobalValue(i, value);
            // } else {
            //     implementation().sumIntoGlobalValue(i, value);
            // }
        }


        inline void c_set(const SizeType &i, const Scalar &value) override
        {
            assert(false && "IMPLEMENT ME");
        }


        inline void c_add(const SizeType &i, const Scalar &value) override
        {
            assert(false && "IMPLEMENT ME");
        }

        inline void set(const Scalar &value) override
        {
            implementation().putScalar(value);
        }

        template<typename Integer>
        void get(
            const std::vector<Integer> &index,
            std::vector<Scalar> &values) const
        {
            // m_utopia_warning_once(" > get does not work in parallel if it asks for ghost entries");
            //FIXME does not work in parallel
            auto data  = get_read_only_data();
            // auto offset = implementation().getMap()->getMinGlobalIndex();

            auto n = index.size();
            values.resize(n);

            for(std::size_t i = 0; i < n; ++i) {
                auto local_index = this->local_index(index[i]);// - offset;
                assert(local_index < data.size());
                values[i] = data[local_index];
            }
        }

        void set_vector(
            const std::vector<GO> &indices,
            const std::vector<Scalar> &values);

        void add_vector(
            const std::vector<GO> &indices,
            const std::vector<Scalar> &values);

        template<typename Integer>
        void generic_select(
            const std::vector<Integer> &index,
            TpetraVector &out
            ) const
        {
            //FIXME does not work in parallel
            auto data   = implementation().getData();
            auto offset = implementation().getMap()->getMinGlobalIndex();

            auto map = Teuchos::rcp(
                new map_type(
                    Teuchos::OrdinalTraits<Tpetra::global_size_t>::invalid(),
                    index.size(),
                    0,
                    communicator())
            );


            out.init(map);

            // if(communicator()->getSize() == 1) {

                std::size_t n = data.size();
                auto out_data = out.implementation().getDataNonConst();

                for(std::size_t i = 0; i < n; ++i) {
                    auto local_index = index[i] - offset;
                    assert(local_index < n);
                    out_data[i] = data[local_index];
                }
            // } else {
            //     /////////////////////////////////////////////////

            //     std::vector<GO> tpetra_index;
            //     tpetra_index.reserve(index.size());

            //     for(auto i : index) {
            //         tpetra_index.push_back(i);
            //     }

            //     const Teuchos::ArrayView<const GO>
            //        index_view(tpetra_index);

            //      auto import_map = Teuchos::rcp(new map_type(global_size, index_view, 0, comm));

            //     Tpetra::Import<
            //         LO,
            //         GO,
            //         vector_type::node_type> importer(map, import_map);

            //     implementation().doImport(out.implementation(), importer, Tpetra::INSERT);
            // }
        }

        inline Teuchos::ArrayRCP<const Scalar> get_read_only_data() const
        {
            if(!ghosted_vec_.is_null()) {
                return ghosted_vec_->getData();
            } else {
                return implementation().getData();
            }
        }


        inline void read_lock() override
        {
            // read_only_data_ = get_read_only_data();
            make_view();
        }

        inline void read_unlock() override
        {
            // read_only_data_ = Teuchos::ArrayRCP<const Scalar>();
            free_view();
        }

        inline void write_lock(WriteMode mode = utopia::AUTO) override
        {
            if(mode != GLOBAL_ADD || mode == GLOBAL_INSERT) {
                make_view();
            }
        }

        void write_unlock(WriteMode mode = utopia::AUTO) override;

        inline void read_and_write_lock(WriteMode mode = utopia::AUTO) override
        {
            // write_data_ = implementation().getDataNonConst();
            // read_only_data_ = write_data_;
            make_view();
        }

        inline void read_and_write_unlock(WriteMode mode = utopia::AUTO) override
        {
            // write_data_ = Teuchos::ArrayRCP<Scalar>();
            // read_only_data_ = Teuchos::ArrayRCP<const Scalar>();

            free_view();
        }

        inline Range range() const
        {
            return { implementation().getMap()->getMinGlobalIndex(), implementation().getMap()->getMaxGlobalIndex() + 1 };
        }

        inline SizeType size() const
        {
            if(is_null()) {
                return {0};
            }

            return implementation().getMap()->getGlobalNumElements();
        }

        inline SizeType local_size() const
        {
            if(is_null()) {
                return {0};
            }

            return implementation().getMap()->getNodeNumElements();
        }

        inline Scalar norm2() const {
           return implementation().norm2();
        }

        inline Scalar norm1() const {
           return implementation().norm1();
        }

        inline Scalar norm_infty() const {
          return implementation().normInf();
        }

 

        bool has_nan_or_inf() const;

        inline void scale(const Scalar alpha)
        {
            implementation().scale(alpha);
        }

        inline Scalar dot(const TpetraVector &other) const
        {
            return implementation().dot(other.implementation());
        }


        // inline void e_mul(const TpetraVector &right, TpetraVector &result) const
        // {
        //     //https://trilinos.org/docs/dev/packages/tpetra/doc/html/classTpetra_1_1MultiVector.html#a95fae4b1f2891d8438b7fb692a85b3bd
        //     result.values(this->communicator(), local_size().get(0), size().get(0), 0.);

        //     result.implementation().elementWiseMultiply(
        //         1.,
        //         this->implementation(),
        //         right.implementation(),
        //         0.
        //     );
        // }

        template<typename Op>
        inline void apply(const Op op)
        {
            KokkosEvalUnary<TpetraVector, Op>::eval(op, *this);
        }

        void reciprocal(TpetraVector &result) const
        {
            if(result.empty() || result.size() != this->size())
            {
                result.init(this->implementation().getMap());
            }

            result.implementation().reciprocal(this->implementation());
        }

        template<typename Op>
        inline void apply_binary(const Op op, const TpetraVector &rhs, TpetraVector &result) const
        {
            KokkosEvalBinary<TpetraVector, Op>::eval(*this, op, rhs, result);
        }

        inline vector_type &implementation()
        {
            assert(!vec_.is_null());
            return *vec_;
        }

        inline const vector_type &implementation() const
        {
            assert(!vec_.is_null());
            return *vec_;
        }

        inline rcpvector_type &raw_type()
        {
            return implementation_ptr();
        }

        inline const rcpvector_type &raw_type() const
        {
            return implementation_ptr();
        }


        inline rcpvector_type &implementation_ptr()
        {
            assert(!vec_.is_null());
            return vec_;
        }

        inline const rcpvector_type &implementation_ptr() const
        {
            assert(!vec_.is_null());
            return vec_;
        }

        inline bool is_null() const
        {
            return vec_.is_null();
        }

        bool read(const Teuchos::RCP< const Teuchos::Comm< int > > &comm, const std::string &path);
        bool write(const std::string &path) const;

        inline bool empty() const
        {
            return vec_.is_null();
        }

        void update_ghosts();
        void export_ghosts_add();

        bool has_ghosts() const
        {
            return !ghosted_vec_.is_null();
        }

        TrilinosCommunicator &comm() override
        {
            return comm_;
        }

        const TrilinosCommunicator &comm() const override
        {
            return comm_;
        }

        void clear() override;

    private:

        TrilinosCommunicator comm_;
        rcpvector_type vec_;
        rcpvector_type ghosted_vec_;
        // Teuchos::ArrayRCP<const Scalar> read_only_data_;
        // Teuchos::ArrayRCP<Scalar> write_data_;

        class View {
        public:
            using DualViewType = vector_type::dual_view_type;
            using HostViewType = DualViewType::t_host;
            using LocalMapType = vector_type::map_type::local_map_type;

            HostViewType view;
            LocalMapType map;

            View(HostViewType &&view, LocalMapType &&map)
            : view(std::move(view)), map(std::move(map)) {}
        };

        std::unique_ptr<View> view_ptr_;

        inline void make_view()
        {
            if(!view_ptr_) {

                if(has_ghosts()) {
                    view_ptr_ = utopia::make_unique<View>(
                        ghosted_vec_->getLocalView<Kokkos::HostSpace>(),
                        ghosted_vec_->getMap()->getLocalMap()
                    );

                } else {
                    view_ptr_ = utopia::make_unique<View>(
                        vec_->getLocalView<Kokkos::HostSpace>(),
                        vec_->getMap()->getLocalMap()
                    );
                }
            }
        }

        inline void free_view()
        {
            view_ptr_ = nullptr;
        }
    };
}

#endif //UTOPIA_TPETRA_VECTOR_HPP
