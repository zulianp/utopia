#ifndef TRIANGULATE_H
#define TRIANGULATE_H

#include <vector>
#include "utopia_Polygon2.hpp"

void triangulate_polygon(const int n_vertices, const double *in_polygon, std::vector<int> &result);
void triangulate_polygon(const utopia::Polygon2 &in_polygon, std::vector<int> &result);

#endif  // TRIANGULATE_H