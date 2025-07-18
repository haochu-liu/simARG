

ARG <- testClonalOrigin_ARG_based(5L, 1, 10L, 1, optimise_recomb=TRUE)

for (i in 1:100) {
  set.seed(i)
  ARG <- testClonalOrigin_ARG_based(6L, 0.5, 10L, 3, optimise_recomb=TRUE)
  tree1 <- local_tree(ARG, 1L)
  tree3 <- local_tree(ARG, 3L)
  tree6 <- local_tree(ARG, 6L)

  if (tree1$sum_time != tree3$sum_time |
      tree3$sum_time != tree6$sum_time |
      tree6$sum_time != tree1$sum_time) {
    print(i)
  }
}

set.seed(4)
ARG <- testClonalOrigin_ARG_based(6L, 0.5, 10L, 3, optimise_recomb=TRUE)


ARG <- testClonalOrigin_tree_based(5L, 1, 10L, 3)

for (i in 1:100) {
  set.seed(i)
  ARG <- testClonalOrigin_tree_based(6L, 1, 10L, 3)
  tree1 <- local_height_ClonalOrigin(ARG, 1)
  tree3 <- local_height_ClonalOrigin(ARG, 3)
  tree6 <- local_height_ClonalOrigin(ARG, 6)

  if (tree1 != tree3 |
      tree3 != tree6 |
      tree6 != tree1) {
    print(i)
  }
}

set.seed(2)
ARG <- testClonalOrigin_tree_based(6L, 1, 10L, 3)
local_height_ClonalOrigin(ARG, 1)
local_height_ClonalOrigin(ARG, 3)
local_height_ClonalOrigin(ARG, 6)


