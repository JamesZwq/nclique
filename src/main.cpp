//
// Created by 张文谦 on 25-3-18.
//
#include<assert.h>
#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<time.h>
#include<limits.h>
#include<unistd.h>
#include<libgen.h>
#include <tree/NCliqueCoreDecomposition.h>

#include"misc.h"
#include"LinkedList.h"
#include"MemoryManager.h"


int main(int argc, char **argv) {
    // usage: ./main treeFile r s
    if (argc < 4) {
        printf("Incorrect number of arguments.\n");
        printf("./main <treeFile> <r> <s>\n");
        printf("treeFile: path to file\n");
        printf("r: r\n");
        printf("s: s\n");
        return 0;
    }
    const char *fpath = argv[1];
    const daf::CliqueSize r = strtol(argv[2], nullptr, 10);
    const daf::CliqueSize s = strtol(argv[3], nullptr, 10);
    std::cout << fpath << " " << r << " " << s << std::endl;
    populate_nCr();

    auto tree = MultiBranchTree::deserialize(fpath);
    // tree->printTree();
    tree->cliqueCount();
    baseNucleusCoreDecompositionPar(*tree,3);

    // baseNucleusCoreDecomposition(*tree,3);
}