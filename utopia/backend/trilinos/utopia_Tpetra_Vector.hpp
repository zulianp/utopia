
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
#include "utopia_BLAS_Operands.hpp"
#include "utopia_Normed.hpp"
#include "utopia_Allocations.hpp"
#include "utopia_Select.hpp"
#include "utopia_RangeDevice.hpp"
#include "utopia_Layout.hpp"

#include "utopia_kokkos_Eval_Binary.hpp"
#include "utopia_kokkos_Eval_Unary.hpp"

#include <Tpetra_Map_decl.hpp>
#include <Tpetra_Vector_decl.hpp>
#include <Tpetra_CrsMatrix.hpp>

#include "utopia_trilinos_Traits.hpp"
#include "utopia_trilinos_Communicator.hpp"

#include <memory>

namespace utopia {

    class TpetraVector :
    public DistributedVector<TpetraScalar, TpetraSizeType>,
    public Normed<TpetraScalar>,
    public Transformable<TpetraScalar>,
    public Reducible<TpetraScalar>,
    public Constructible<TpetraScalar, TpetraSizeType, 1>,
    public ElementWiseOperand<TpetraScalar>,
    public ElementWiseOperand<TpetraVector>,
    public Comparable<TpetraVector>,
    public BLAS1Tensor<TpetraVector>,
    public Tensor<TpetraVector, 1>,
    public Selectable<TpetraVector, 1>
    {
    public:

        using SizeType      = Traits<TpetraVector>::SizeType;
        using LocalSizeType = Traits<TpetraVector>::LocalSizeType;
        using Scalar        = Traits<TpetraVector>::Scalar;
        using IndexSet      = Traits<TpetraVector>::IndexSet;
        using IndexArray    = Traits<TpetraVector>::IndexArray;
        using ScalarArray   = Traits<TpetraVector>::ScalarArray;
        using Node          = Traits<TpetraVector>::Node;

        using MapType         = Tpetra::Map<LocalSizeType, SizeType, Node>;
        using VectorType      = Tpetra::Vector<Scalar, LocalSizeType, SizeType, Node>;
        using MultiVectorType = Tpetra::MultiVector<Scalar, LocalSizeType, SizeType, Node>;
        using RCPVectorType   = Teuchos::RCP<VectorType>;
        using RCPCommType     = Teuchos::RCP<const Teuchos::Comm<int> > ;
        using RCPMapType      = Teuchos::RCP<const MapType>;
        using ExecutionSpace  = VectorType::execution_space;
        using Layout          = utopia::Layout<TrilinosCommunicator, SizeType, 1>;

        ////////////////////////////////////////////////////////////////////
        ///////////////////////// BOILERPLATE CODE FOR EDSL ////////////////
        ////////////////////////////////////////////////////////////////////

        using Super         = utopia::Tensor<TpetraVector, 1>;
        using Constructible = utopia::Constructible<Scalar, SizeType, 1>;
        using Constructible::zeros;

        using Super::Super;

         template<class Expr>
         TpetraVector(const Expression<Expr> &expr)
         {
            static_assert(!std::is_same<TpetraVector, Expr>::value, "should not come here with this derived type");
             //THIS HAS TO BE HERE IN EVERY UTOPIA TENSOR CLASS
             Super::construct_eval(expr.derived());
         }

         template<class Expr>
         inline TpetraVector &operator=(const Expression<Expr> &expr)
         {
             static_assert(!std::is_same<TpetraVector, Expr>::value, "should not come here with this derived type");
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

        RCPCommType communicator()
        {
            return implementation().getMap()->getComm();
        }

        const RCPCommType communicator() const
        {
            return implementation().getMap()->getComm();
        }

        TpetraVector &operator=(const TpetraVector &other);

        TpetraVector &operator=(TpetraVector &&other)
        {
            if(this == &other) return *this;
            comm_ = std::move(other.comm_);
            vec_ = std::move(other.vec_);
            ghosted_vec_ = std::move(other.ghosted_vec_);
            return *this;
        }

        void select(const IndexSet &index, TpetraVector &result) const override;

        void copy(const TpetraVector &other) override;

        /////////////////////////////////////////////////////////////////////////////////
        //////////////////////////////// OVERRIDES for Constructible ////////////////////
        /////////////////////////////////////////////////////////////////////////////////

        inline void values(const Layout &l, const Scalar &value)
        {
            values(l.comm(), l.local_size(), l.global_size(), value);
        }

        inline void zeros(const Layout &l)
        {
            values(l, 0.0);
        }

        void values(const SizeType &s, const Scalar &val) override;
        void local_values(const SizeType &s, const Scalar &val) override;

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

        void e_div(const TpetraVector &other) override;
        void e_min(const TpetraVector &other) override;
        void e_max(const TpetraVector &other) override;
        ////////////////////////////////////////////////////////////////

        inline void e_mul(const Scalar &other) override
        {
            scale(other);
        }

        void e_div(const Scalar &other) override;
        void e_min(const Scalar &other) override;
        void e_max(const Scalar &other) override;

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

        void zeros(
            const TrilinosCommunicator &comm,
            const SizeType &n_local,
            const SizeType &n_global)
        {
            values(comm.get(), n_local, n_global, 0.0);
        }

        void values(
            const TrilinosCommunicator &comm,
            const SizeType &n_local,
            const SizeType &n_global,
            const Scalar &value)
        {
            values(comm.get(), n_local, n_global, value);
        }

        void values(
            const RCPCommType &comm,
            const SizeType &n_local,
            const SizeType &n_global,
            const Scalar &value);

        inline void init(const RCPMapType &map)
        {
            UTOPIA_REPORT_ALLOC("TpetraVector::init");
            vec_.reset(new VectorType(map));
        }

        void ghosted(const TpetraVector::SizeType &local_size,
                     const TpetraVector::SizeType &global_size,
                     const std::vector<SizeType> &ghost_index
        )
        {
            ghosted(comm().get(), local_size, global_size, ghost_index);
        }

        void ghosted(const RCPCommType &comm,
                     const TpetraVector::SizeType &local_size,
                     const TpetraVector::SizeType &global_size,
                     const std::vector<SizeType> &ghost_index
        );

        inline void axpy(const Scalar &alpha, const TpetraVector &x) override
        {
            implementation().update(alpha, *x.vec_, 1.);
        }

        inline void swap(TpetraVector &x) override
        {
            using std::swap;

            swap(comm_, x.comm_);
            swap(vec_, x.vec_);
            swap(ghosted_vec_, x.ghosted_vec_);
        }

        inline void describe(std::ostream &os) const
        {
            auto out = Teuchos::getFancyOStream(Teuchos::rcpFromRef(os));
            implementation().describe(*out, Teuchos::EVerbosityLevel::VERB_EXTREME);
        }

        inline void describe() const override
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

        inline Scalar l_get(const SizeType &i) const //override
        {
            return view_ptr_->view(i, 0);
        }

        inline Scalar operator[](const SizeType &i) const
        {
            // assert(!read_only_data_.is_null() && "Use Read<Vector> w(v); to enable reading from this vector v!");
            return get(i);
        }

        inline LocalSizeType local_index(const SizeType &i) const
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
        }

        inline void l_set(const SizeType &i, const Scalar &value) //override
        {
            assert(view_ptr_);
            view_ptr_->view(i, 0) = value;
        }

        inline void add(const SizeType &i, const Scalar &value) override
        {
            assert(view_ptr_);

            auto local_index = view_ptr_->map.getLocalElement(i);

            view_ptr_->view(local_index, 0) += value;
        }

        void c_set(const SizeType &i, const Scalar &value) override;
        void c_add(const SizeType &i, const Scalar &value) override;

        inline void set(const Scalar &value) override
        {
            implementation().putScalar(value);
        }

        template<typename Integer>
        void get(
            const std::vector<Integer> &index,
            ScalarArray &values) const
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
            const IndexArray &indices,
            const ScalarArray &values);

        void add_vector(
            const IndexArray &indices,
            const ScalarArray &values);

        template<typename Integer>
        void generic_select(
            const std::vector<Integer> &index,
            TpetraVector &out
        ) const;

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

        inline Range range() const override
        {
            return { implementation().getMap()->getMinGlobalIndex(), implementation().getMap()->getMaxGlobalIndex() + 1 };
        }

        inline RangeDevice<TpetraVector> range_device() const
        {
            return RangeDevice<TpetraVector>(
                implementation().getMap()->getMinGlobalIndex(),
                implementation().getMap()->getMaxGlobalIndex() + 1
            );
        }

        inline SizeType size() const override
        {
            if(is_null()) {
                return 0;
            }

            return implementation().getMap()->getGlobalNumElements();
        }

        inline SizeType local_size() const override
        {
            if(is_null()) {
                return 0;
            }

            return implementation().getMap()->getNodeNumElements();
        }

        inline Scalar norm2() const override {
           return implementation().norm2();
        }

        inline Scalar norm1() const override {
           return implementation().norm1();
        }

        inline Scalar norm_infty() const override {
          return implementation().normInf();
        }

        bool has_nan_or_inf() const;

        inline void scale(const Scalar &alpha) override
        {
            implementation().scale(alpha);
        }

        inline Scalar dot(const TpetraVector &other) const override
        {
            return implementation().dot(other.implementation());
        }

        template<typename Op>
        inline void apply(const Op op)
        {
            KokkosEvalUnary<TpetraVector, Op>::eval(op, *this);
        }

        // void reciprocal(TpetraVector &result) const
        // {
        //     if(result.empty() || result.size() != this->size())
        //     {
        //         result.init(this->implementation().getMap());
        //     }

        //     result.implementation().reciprocal(this->implementation());
        // }

        template<typename Op>
        inline void apply_binary(const Op op, const TpetraVector &rhs, TpetraVector &result) const
        {
            KokkosEvalBinary<TpetraVector, Op>::eval(*this, op, rhs, result);
        }

        inline VectorType &implementation()
        {
            assert(!vec_.is_null());
            return *vec_;
        }

        inline const VectorType &implementation() const
        {
            assert(!vec_.is_null());
            return *vec_;
        }

        inline RCPVectorType &raw_type()
        {
            return implementation_ptr();
        }

        inline const RCPVectorType &raw_type() const
        {
            return implementation_ptr();
        }

        inline RCPVectorType &implementation_ptr()
        {
            assert(!vec_.is_null());
            return vec_;
        }

        inline const RCPVectorType &implementation_ptr() const
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

        inline bool empty() const  override
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


        inline std::string get_class() const override
        {
            return "TpetraVector";
        }

        inline bool is_alias(const TpetraVector &other) const
        {
            return vec_ == other.vec_;
        }

    private:
        TrilinosCommunicator comm_;
        RCPVectorType vec_;
        RCPVectorType ghosted_vec_;

        class View {
        public:
            using DualViewType = VectorType::dual_view_type;
            using HostViewType = DualViewType::t_host;
            using LocalMapType = VectorType::map_type::local_map_type;

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
                        ghosted_vec_->getLocalViewHost(),
                        ghosted_vec_->getMap()->getLocalMap()
                    );

                } else {
                    view_ptr_ = utopia::make_unique<View>(
                        vec_->getLocalViewHost(),
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
