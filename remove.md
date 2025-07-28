0 1 2 3 4 5 {0 1 2}

CONF {0 1 2} 3
0 1 3 4 5     -- 跳过2因为size满了 1
{2} 0 3 4 5   -- 上轮跳过了2, 所以2是keep节点 2
{1 2} 3 4 5   -- 上轮跳过了1, 所以1是keep节点 3



[1,2,3,4,5,6,7,8]
[1,2,3,4,5,6,7]

[1,2,3,4,5,6]
[1,2,3,4,5] [7]
[1,2,3,4] [6,7]
[1,2,3] [5,6,7]
[1,2] [4,5,6,7]
[1] [3,4,5,6,7]
[] [2,3,4,5,6,7]






















1 2 3 4 5     -- 跳过2因为size满了
{0} 2 3 4 5   -- 上轮跳过了2, 所以2是keep节点
{0 1} 3 4 5   -- 上轮跳过了1, 所以1是keep节点


CONF {0 1 2} {1 2 3}
0 1 2 3 4 5

0 1 3 4 5     ok
{2} 0 3 4 5   ok
{2 1} 3 4 5   not ok, not correct with {1 2 3}
    {2 1} 4 5     ok

CONF {0 1 2} {1 3 4}
0 1 2 3 4 5

0 1 3 4 5     not ok, not correct with {1 3 4}
    0 1 3 5      ok
    {4} 0 1 5    ok
    {4 3} 0 5    ok
{2} 0 3 4 5   ok
{2 1} 3 4 5   not ok
    {2 1} 3 5     ok
    {2 1 4} 5     ok

先对于CONF[1] 进行拆分:
0 1 3 4 5
    发现这个有冲突, 和CONF[2], 基于CONF[2]的拆分:
    0 1 3 5      ok
    {4} 0 1 5    ok
    {4 3} 0 5    ok
{2} 0 3 4 5
    没有冲突
{2 1} 3 4 5












0 1 3 4
0 1 3 {5}

{2} 0 3 4
{2} 0 3 {5}

{2 1} 3 4
{2 1} 3 {5}










0 1 3 5       -- 跳过2因为size 1满了 跳过4因为size 2满了
    {2} 0 3 5     -- 上轮跳过了2, 所以2是keep; 
                        在这里跳过1, 因为size 1满了;  
                        在这里跳过4, 因为size 2满了;
        {2 1} 3 5     -- 上轮跳过了1, 所以1是keep; 
                            在这里跳过0, 因为size 1满了;  
                            在这里跳过4, 因为size 2满了;
            {2 1 4} 5     -- 上轮跳过了4, 所以4是keep; 
                                在这里跳过0, 因为size 1满了;  
                                在这里跳过3, 因为size 2满了;
            -- 0 3都无法在加入到hold, 本次递归结束.
        {2 4} 0 5     -- 上轮跳过了4, 所以4是keep; 
                            在这里跳过1, 因为size 1满了;  
                            在这里跳过3, 因为size 2满了;
            {2 4 1} 5
    {4} 0 1 5 -- 上轮跳过了4, 所以4是keep; 
                        在这里跳过2, 因为size 1满了;
                        在这里跳过3, 因为size 2满了;
        {4 2} 0 5 -- 上轮跳过了2, 所以2是keep; 
                        在这里跳过1, 因为size 1满了;
                        在这里跳过3, 因为size 2满了; 
            {4 2 1} 5 -- 上轮跳过了1, 所以1是keep; 
                        在这里跳过0, 因为size 1满了;
                        在这里跳过3, 因为size 2满了;




CONF {0 1 2} {2 3}
0 1 3 4 5
{2} 1 4 5
{2 3} 1 4 5
{0 2} 4 5
{0 2 3} 4 5

CONF {0 1 2 3} 
0 1 2 4
｛3｝ 1 2 4
｛0 3｝ 1 4
｛0 2 3｝ 4

CONF {0 1 2 3} 
0 1 2 4
｛3｝ 1 2 4
｛0 3｝ 1 4
｛0 2 3｝ 4


CONF {0 1 2 3} 
1 2 3 4
｛0｝ 2 3 4
｛0 1｝ 3 4
｛0 1 2｝ 4

0 1 2 3 4 5
CONF {0 1 2 3}  {0,1,2,4}
0 1 2 5
{3} 0 1 5
{4} 0 1 5
{3 4} 0 1 5



std::vector<TreeGraphNode> vList = {
    {2, false},
    {4, true},
    {5, true},[BronKerboschRmRClique.hpp](src/BK/BronKerboschRmRClique.hpp)
    {6, true},
    {9, true},[BronKerboschRmRClique.hpp](src/BK/BronKerboschRmRClique.hpp)
};
std::vector<std::vector<daf::Size>> conflictSets;
conflictSets.emplace_back(std::vector<daf::Size>{5, 6, 9});
conflictSets.emplace_back(std::vector<daf::Size>{4, 6, 9});
conflictSets.emplace_back(std::vector<daf::Size>{2, 4, 6});
conflictSets.emplace_back(std::vector<daf::Size>{2, 5, 6});
conflictSets.emplace_back(std::vector<daf::Size>{2, 6, 9});
conflictSets.emplace_back(std::vector<daf::Size>{4, 5, 6});


{2} 4 5 6 9

{2} 4 5 6 // 2 4 6
    {2} 4 5
    {2 6} 5 // 2 5 6
        {2 6}
{2 9} 4 5


[//]: # ({2 6 9} 4 // 4 6 9)

[//]: # (    {2 6 9})


changed leafId: 2 leaf index: 0 leaf: vec size: 6 [(3, Keep), (4, Drop), (6, Drop), (7, Drop), (8, Drop), (9, Drop)]
removedRCliqueIdForLeaf: 28 ([3, 4, 9]) 31 ([3, 6, 9]) 44 ([7, 8, 9]) 43 ([6, 8, 9]) 33 ([3, 7, 9]) 34 ([3, 8, 9]) 37 ([4, 6, 9]) 40 ([4, 8, 9]) 39 ([4, 7, 9]) 42 ([6, 7, 9]) 

{3} 4 6 7 9

[//]: #  4 6 9
{3} 6 7 9
{3 4} 6 7
}



std::vector<TreeGraphNode> vList = {
    {1, false},  // 原来是 3
    {2, true},   // 原来是 4
    {3, true},   // 原来是 6
    {4, true},   // 原来是 7
    {5, true},   // 原来是 8
    {6, true},   // 原来是 9
};

std::vector<std::vector<daf::Size>> conflictSets;
conflictSets.emplace_back(std::vector<daf::Size>{1, 2, 6});  // {3,4,9}
conflictSets.emplace_back(std::vector<daf::Size>{1, 3, 6});  // {3,6,9}
conflictSets.emplace_back(std::vector<daf::Size>{4, 5, 6});  // {7,8,9}
conflictSets.emplace_back(std::vector<daf::Size>{3, 5, 6});  // {6,8,9}
conflictSets.emplace_back(std::vector<daf::Size>{1, 4, 6});  // {3,7,9}[NucleusCoreDecompositionRemoveSclique.cpp](src/NucleusDecomposition/NucleusCoreDecompositionRemoveSclique.cpp)
conflictSets.emplace_back(std::vector<daf::Size>{1, 5, 6});  // {3,8,9}
conflictSets.emplace_back(std::vector<daf::Size>{2, 3, 6});  // {4,6,9}
conflictSets.emplace_back(std::vector<daf::Size>{2, 5, 6});  // {4,8,9}
conflictSets.emplace_back(std::vector<daf::Size>{2, 4, 6});  // {4,7,9}
conflictSets.emplace_back(std::vector<daf::Size>{3, 4, 6});  // {6,7,9}

{1} 2 3 4 5 6
{1} 3 4 5 6
{1 2} 3 4 5 