#ifndef UTOPIA_LIBMESH_SIDESET_HPP
#define UTOPIA_LIBMESH_SIDESET_HPP

namespace utopia {
    namespace libmesh {
        class SideSet {
        public:
            using BoundaryIdType = int;

            inline static constexpr BoundaryIdType invalid_id() { return -1; }

            class Cube {
            public:
                // back: 0
                // bottom: 1
                // right: 2
                // top: 3
                // left: 4
                // front: 5

                inline static constexpr BoundaryIdType back_id() { return 0; }
                inline static constexpr BoundaryIdType bottom_id() { return 1; }
                inline static constexpr BoundaryIdType right_id() { return 2; }

                inline static constexpr BoundaryIdType top_id() { return 3; }
                inline static constexpr BoundaryIdType left_id() { return 4; }
                inline static constexpr BoundaryIdType front_id() { return 5; }

                inline static SideSet left() { return SideSet(left_id(), "left"); }
                inline static SideSet bottom() { return SideSet(bottom_id(), "bottom"); }
                inline static SideSet back() { return SideSet(back_id(), "back"); }

                inline static SideSet right() { return SideSet(right_id(), "right"); }
                inline static SideSet top() { return SideSet(top_id(), "top"); }
                inline static SideSet front() { return SideSet(front_id(), "front"); }
                inline static SideSet invalid() { return SideSet(invalid_id(), ""); }

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
            };

            inline SideSet(const std::string &name) : name_(name) {}
            inline SideSet(const int id) : id_(id) {}
            inline SideSet(const int id, const std::string &name) : id_(id), name_(name) {}

            inline operator int() const { return id_; }
            // inline operator const std::string &() const { return name_; }

        private:
            int id_{invalid_id()};
            std::string name_;
        };
    }  // namespace libmesh
}  // namespace utopia

#endif  // UTOPIA_LIBMESH_SIDESET_HPP
