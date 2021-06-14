#ifndef UTOPIA_STK_SIDESET_HPP
#define UTOPIA_STK_SIDESET_HPP

namespace utopia {
    namespace stk {
        class SideSet {
        public:
            using BoundaryIdType = int;

            class Cube {
            public:
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

                inline static std::string name_from_id(const BoundaryIdType id) {
                    return "surface_" + std::to_string(id);
                }

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

            inline static SideSet invalid() { return SideSet(""); }

            inline SideSet(int num) : name_("surface_" + std::to_string(num)) {}
            inline SideSet(std::string name) : name_(std::move(name)) {}

            inline operator const std::string &() const { return name_; }

        private:
            std::string name_;
        };
    }  // namespace stk
}  // namespace utopia

#endif  // UTOPIA_STK_SIDESET_HPP
