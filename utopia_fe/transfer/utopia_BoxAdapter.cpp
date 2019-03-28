#include "utopia_BoxAdapter.hpp"

// namespace utopia {


// 	void BoxAdapter::read(moonolith::InputStream &is)
// 	{
// 	    auto &min = get_min();
// 	    auto &max = get_max();

// 	    int n;
// 	    is >> n;

// 	    reset(n);

// 	    for (int i = 0; i < n; ++i) {
// 	        min.el(i) << is;
// 	        max.el(i) << is;
// 	    }
// 	}

// 	void BoxAdapter::write(moonolith::OutputStream &os) const
// 	{
// 	    const int n = get_dims();
// 	    auto &min = get_min();
// 	    auto &max = get_max();

// 	    os << n;

// 	    for (int i = 0; i < n; ++i) {
// 	        os << min(i);
// 	        os << max(i);
// 	    }
// 	}
// }
