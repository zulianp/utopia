// #ifndef UTOPIA_HOMEMADE_MESH_HPP
// #define UTOPIA_HOMEMADE_MESH_HPP

// #include "utopia_Base.hpp"
// #include "utopia_homemade_FEForwardDeclarations.hpp"

// #include <vector>
// #include <memory>

// namespace utopia {

//     class Mesh {
//     public:
//         class Impl;

//         void make_triangle(const int order = 1);
//         int element_order(const int element_index) const;
//         int element_type(const int element_index) const;
//         void node_indices(const int elem, std::vector<int> &index);
//         int n_dims() const;

//         std::size_t n_nodes() const;
//         std::size_t n_elements() const;

//         Mesh();
//         ~Mesh();

//         void *mesh_impl_ptr() const;

//         //memory
//         std::vector<int> el_ptr;
//         std::vector<int> el_index;
//         std::vector<int> el_type;
//         std::vector<int> meta;
//         std::vector<double> points;

//         inline friend auto elements_begin(const Mesh &m) -> int
//         {
//             return 0;
//         }

//         inline friend auto elements_end(const Mesh &m) -> int
//         {
//             return m.n_elements();
//         }

//     private:
//         std::unique_ptr<Impl> impl_ptr;
//     };

// }

// #endif //UTOPIA_HOMEMADE_MESH_HPP
