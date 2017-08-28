#ifndef UTOPIA_BOX_ADAPTER_HPP
#define UTOPIA_BOX_ADAPTER_HPP 


#include "Box.hpp"
#include "moonolith_serializable.hpp"
#include "moonolith_describable.hpp"
#include "moonolith_input_stream.hpp"
#include "moonolith_output_stream.hpp"
#include "moonolith_bounding_volume_with_span.hpp"


namespace utopia {

	    // class BoxAdapter : public moonolith::Serializable, public moonolith::Describable, public Box {
	    // public:
	    //     void read(moonolith::InputStream &is) override;
	    //     void write(moonolith::OutputStream &os) const override;
	        
	    //     void describe(std::ostream &os) const override
	    //     {
	    //         print(os);
	    //     }
	        
	    //     inline bool isEmpty() const
	    //     {
	    //         return empty();
	    //     }
	        
	    //     inline double getMinAt(const int coord) const
	    //     {
	    //         return get_min(coord);
	    //     }
	        
	    //     inline double getMaxAt(const int coord) const
	    //     {
	    //         return get_max(coord);
	    //     }
	        
	    //     inline void setMinAt(const int coord, const double value)
	    //     {
	    //         get_min().el(coord) = value;
	    //     }
	        
	    //     inline void setMaxAt(const int coord, const double value)
	    //     {
	    //         get_max().el(coord) = value;
	    //     }
	        
	    //     inline void clear()
	    //     {
	    //         reset();
	    //     }
	        
	    //     inline int nDims() const {
	            
	    //         return get_dims();
	    //     }
	    // };
	    
	    // template<int Dimension>
	    // class BoxBoxAdapter : public moonolith::Describable, public moonolith::Serializable {
	    // public:
	    //     typedef utopia::BoxAdapter StaticBound;
	        
	    //     void read(moonolith::InputStream &is)
	    //     {
	    //         is >> static_;
	    //         bool is_empty;
	    //         is >> is_empty;
	    //         if(!is_empty) { is >> dynamic_; };
	    //     }
	        
	        
	    //     void write(moonolith::OutputStream &os) const
	    //     {
	    //         os << static_;
	    //         bool is_empty = dynamic_.isEmpty();
	    //         os << is_empty;
	    //         if(!is_empty) { os << dynamic_; }
	    //     }
	        
	    //     bool intersects(const BoxBoxAdapter &bound) const
	    //     {
	    //         return static_.intersects(bound.static_) && dynamic_.intersects(bound.dynamic_);
	    //     }
	        
	    //     bool intersects(const BoxBoxAdapter &bound, const double tol) const
	    //     {
	    //         return static_.intersects(bound.static_, tol) && dynamic_.intersects(bound.dynamic_, tol);
	    //     }
	        
	    //     bool intersects(const BoxAdapter &bound) const
	    //     {
	    //         return static_.intersects(bound);
	    //     }
	        
	    //     inline double getMinAt(const int coord) const
	    //     {
	    //         return static_.getMinAt(coord);
	    //     }
	        
	    //     inline double getMaxAt(const int coord) const
	    //     {
	    //         return static_.getMaxAt(coord);
	    //     }
	        
	    //     inline void setMinAt(const int coord, const double value)
	    //     {
	    //         static_.setMinAt(coord, value);
	    //     }
	        
	    //     inline void setMaxAt(const int coord, const double value)
	    //     {
	    //         static_.setMaxAt(coord, value);
	    //     }
	        
		   //     //expands to contain the union of this and CompositeBound
	    //     BoxBoxAdapter &operator +=(const BoxBoxAdapter &bound)
	    //     {
	    //         static_ += bound.static_;
	    //         if(dynamic_.isEmpty()) {
	    //             dynamic_ = bound.dynamic_;
	    //         } else if(!bound.dynamic_.isEmpty()) {
	    //             dynamic_ += bound.dynamic_;
	    //         }
	    //         return *this;
	    //     }
	        
	    //     bool isEmpty() const
	    //     {
	    //         return static_.isEmpty();
	    //     }
	        
	    //     void clear()
	    //     {
	    //         static_.reset(Dimension);
	    //         dynamic_.reset(Dimension);
	    //     }
	        
	    //     BoxBoxAdapter()
	    //     {
	    //         clear();
	    //     }
	        
	    //     void describe(std::ostream &os) const
	    //     {
	    //         os << "Static bound:\n"  << static_  << "\n";
	    //         os << "Dynamic bound:\n";
	    //         dynamic_.describe(os);
	    //         os << "\n";
	    //     }
	        
	    //     inline BoxAdapter &staticBound() { return static_; }
	    //     inline const BoxAdapter &staticBound() const { return static_; }
	        
	    //     inline BoxAdapter &dynamicBound() { return dynamic_; }
	    //     inline const BoxAdapter &dynamicBound() const { return dynamic_; }
	        
	    // private:
	    //     BoxAdapter static_;
	    //     BoxAdapter dynamic_;
	    // };

	template<int Dimension>
	using BoxBoxAdapter = moonolith::AABBWithSpan<Dimension, double>;
}

#endif //UTOPIA_BOX_ADAPTER_HPP

