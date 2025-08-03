#include <iostream>
#include <ranges>
#include <vector>
#include <string>
#include <string_view>
#include <type_traits>

// ─── 概念：可遍历、元素可用 << 输出，排除 string / C-string 等 ───────────
template<class R>
concept printable_range =
       std::ranges::input_range<R> &&
       requires(std::ostream& os, std::ranges::range_reference_t<R> ref) {
           os << ref;
       } &&
       !std::same_as<std::remove_cvref_t<R>, std::string> &&
       !std::same_as<std::remove_cvref_t<R>, std::string_view> &&
       !(std::is_pointer_v<std::remove_cvref_t<R>> &&
         std::is_same_v<std::remove_cvref_t<std::remove_pointer_t<R>>, char>) &&
       !(std::is_array_v<std::remove_cvref_t<R>> &&
         std::is_same_v<std::remove_all_extents_t<std::remove_cvref_t<R>>, char>);

// ─── 包装器：让 operator<< 可用在任意 printable_range 上 ───────────────
template<printable_range R>
struct printable_wrapper {
    R r;
};

// 让 CTAD 推导出 printable_wrapper 类型（关键：解决 no deduction guide 错误）
template<printable_range R>
printable_wrapper(R&&) -> printable_wrapper<std::remove_cvref_t<R>>;

// 实际的 operator<< 重载（作用于包装后的 range）
template<printable_range R>
std::ostream& operator<<(std::ostream& os, printable_wrapper<R> const& w) {
    using std::ranges::begin;
    using std::ranges::end;
    using std::ranges::size;
    using value_t = std::remove_cvref_t<std::ranges::range_value_t<R>>;

    constexpr bool elem_is_int = std::integral<value_t>;

    std::size_t n;
    if constexpr (std::ranges::sized_range<R>) {
        n = size(w.r);
    } else {
        n = std::ranges::distance(begin(w.r), end(w.r));
    }

    bool isClique = elem_is_int && n > 1;

    if (isClique) {
        os << "[";
    } else {
        os << "vec size: " << n << " [";
    }

    bool first = true;
    for (auto&& elem : w.r) {
        if (!first) os << ", ";
        first = false;
        os << elem;
    }

    os << "]";
    return os;
}

// ─── 示例用法 ───────────────────────────────────────────────────────
int main() {
    std::vector<unsigned long long> v{1, 2, 3, 4};
    auto transformed = v | std::views::transform([](auto x){ return x * 10; });

    std::cout << "原始 vector: " << printable_wrapper(v) << "\n";
    std::cout << "transform_view: " << printable_wrapper(transformed) << "\n";

    // 结合你的打印语境（类似 conflictSets / maxRClique）
    unsigned long long maxRClique = 42;
    std::cout << "Conflict sets: " << printable_wrapper(transformed)
              << ", maxRClique: " << maxRClique << std::endl;

    // 非整数元素的例子（用于验证 clique 判断分支）
    std::vector<std::string> names{"alice", "bob"};
    std::cout << "Names: " << printable_wrapper(names) << "\n";
}