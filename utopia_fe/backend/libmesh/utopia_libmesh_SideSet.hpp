#ifndef UTOPIA_LIBMESH_SIDESET_HPP
#define UTOPIA_LIBMESH_SIDESET_HPP

namespace utopia {
    namespace libmesh {
        class SideSet {
        public:
            using BoundaryIdType = int;
            inline static constexpr BoundaryIdType left_id() { return 1; }
            inline static constexpr BoundaryIdType bottom_id() { return 2; }
            inline static constexpr BoundaryIdType back_id() { return 3; }

            inline static constexpr BoundaryIdType right_id() { return 4; }
            inline static constexpr BoundaryIdType top_id() { return 5; }
            inline static constexpr BoundaryIdType front_id() { return 6; }
            inline static constexpr BoundaryIdType invalid_id() { return 0; }

            inline static SideSet left() { return SideSet(left_id()); }
            inline static SideSet bottom() { return SideSet(bottom_id()); }
            inline static SideSet back() { return SideSet(back_id()); }

            inline static SideSet right() { return SideSet(right_id()); }
            inline static SideSet top() { return SideSet(top_id()); }
            inline static SideSet front() { return SideSet(front_id()); }
            inline static SideSet invalid() { return SideSet(invalid_id()); }

            inline SideSet(const std::string &name) : id_(id_from_user_space_name(name)) {}
            inline SideSet(const int id) : id_(id) {}

            inline operator int() const { return id_; }

            inline static BoundaryIdType id_from_user_space_name(const std::string &str_in) {
                if (str_in.empty()) {
                    return invalid_id();
                }

                std::string str;
                str.resize(str_in.size(), '0');

                std::transform(
                    str_in.begin(), str_in.end(), str.begin(), [](unsigned char c) { return std::tolower(c); });

                if (str == "left") {
                    return left_id();
                }

                if (str == "right") {
                    return right_id();
                }

                if (str == "bottom") {
                    return bottom_id();
                }

                if (str == "top") {
                    return top_id();
                }

                if (str == "front") {
                    return front_id();
                }

                if (str == "back") {
                    return back_id();
                }

                return invalid_id();
            }

        private:
            int id_;
        };
    }  // namespace libmesh
}  // namespace utopia

#endif  // UTOPIA_LIBMESH_SIDESET_HPP
