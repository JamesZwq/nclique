#include <iostream>
#include <map>

int main() {
    std::map<double, double, std::greater<std::size_t>> myMap;
    myMap[3] = 123;
    myMap[1] = 12341234;
    myMap[2000000000000000000] = 12398471;

    // 按 key 降序遍历
    std::cout << "Descending order:" << std::endl;
    for (const auto& [key, value] : myMap) {
        std::cout << key << ": " << value << std::endl;
    }
    std::cout << ++myMap[10] << std::endl;


    for (const auto& [key, value] : myMap) {
        std::cout << key << ": " << value << std::endl;
    }

    return 0;
}