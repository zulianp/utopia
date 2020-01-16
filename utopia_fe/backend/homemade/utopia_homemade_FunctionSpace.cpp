// #include "utopia_homemade_FunctionSpace.hpp"
// #include <map>

// #include "utopia_intersector.hpp"


// namespace utopia {

//     void HMFESpace::dof_indices(const int element_index, std::vector<int> &indices)
//     {
//         indices.clear();
//         indices.insert(indices.end(),
//                        dof_index.begin() + dof_ptr[element_index],
//                        dof_index.begin() + dof_ptr[element_index + 1]);
//     }

//     void HMFESpace::make_dof_map()
//     {
//         std::vector<int> node_indices;

//         if(order_ == 2) {
//             std::map<std::pair<int, int>, int> edges;
//             dof_index.reserve(mesh().n_nodes() * 2);
//             dof_ptr.resize(mesh().n_elements() + 1, 0);

//             int edge_index = 0;
//             for(std::size_t i = 0; i < mesh().n_elements(); ++i) {
//                 mesh().node_indices(i, node_indices);
//                 dof_index.insert(dof_index.end(), node_indices.begin(), node_indices.end());

//                 switch(mesh().element_type(i)) {
//                     case Intersector::ELEMENT_TYPE_TRIANGLE_ORDER_2:
//                     {
//                         dof_ptr[i + 1] = dof_ptr[i] + node_indices.size() + 3;

//                         for(int k = 0; k < 3; ++k) {
//                             const int kp1 = k == 2? 0 : (k + 1);
//                             if(node_indices[k] > node_indices[kp1]) {
//                                 auto e  = std::pair<int, int>(node_indices[kp1], node_indices[k]);
//                                 auto ret = edges.insert(std::make_pair(e, edge_index));

//                                 dof_index.push_back(ret.first->second);

//                                 if(ret.second) {
//                                     ++edge_index;
//                                 }
//                             }
//                         }

//                         break;
//                     }

//                     default:
//                     {
//                         assert(false);
//                         break;
//                     }
//                 }


//             }
//         } else {
//             dof_index = mesh().el_index;
//             dof_ptr = mesh().el_ptr;
//         }
//     }

// }