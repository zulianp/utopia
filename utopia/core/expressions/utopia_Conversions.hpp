//
// Created by Patrick Zulian on 06/07/15.
//

#ifndef UTOPIA_UTOPIA_CONVERSIONS_HPP
#define UTOPIA_UTOPIA_CONVERSIONS_HPP

#include "utopia_Range.hpp"
#include "utopia_Factory.hpp"
#include "utopia_Base.hpp"

namespace utopia {
//    template<class Vector, class Matrix>
//    void vec2mat(const Vector &v, Matrix &m, const bool transpose = false)
//    {
//        if(transpose) {
//            m = values(1, v.size().get(0), 0);
//            const Read<Vector> read(v);
//            const Write<Matrix> write(m);
//            const Range r = range(v);
//
//            for(int i = r.begin(); i < r.end(); ++i) {
//                m.set(0, i, v.get(i));
//            }
//
//        } else {
//            m = values(v.size().get(0), 1, 0);
//            const Read<Vector> read(v);
//            const Write<Matrix> write(m);
//            const Range r = range(v);
//
//            for(int i = r.begin(); i < r.end(); ++i) {
//                m.set(i, 0, v.get(i));
//            }
//        }
//    }
}

#endif //UTOPIA_UTOPIA_CONVERSIONS_HPP


