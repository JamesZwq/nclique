//
// Created by 张文谦 on 25-4-16.
//

#ifndef EDGELIST_H
#define EDGELIST_H
#include "Graph.h"

class EdgeList {

    public:
        explicit EdgeList(const Graph &graph);


        const Graph &graph;
        unsigned int *edgeList;
        unsigned int *edgeListOffsets;
};

#endif //EDGELIST_H
