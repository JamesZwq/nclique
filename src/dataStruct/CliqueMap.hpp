//
// Created by 张文谦 on 25-5-8.
//

#ifndef CLIQUEMAP_HPP
#define CLIQUEMAP_HPP

#endif //CLIQUEMAP_HPP

#include <Global/Global.h>

class CliqueMap {
    daf::Size n;
    daf::StaticVector<std::unordered_map<daf::Size>> cliqueMap;
    daf::StaticVector<daf::StaticVector<daf::Size>> cliques;
};