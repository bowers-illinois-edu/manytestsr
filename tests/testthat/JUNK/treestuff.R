
### Trying to figure out how to extract the subtrees:
### We want, for example, for the 5th split, the nodes labeled "0","00","000","0000","00000","00001"
### Or, for the 4th split, "0","01","010","0100","0101" (all the 010 plus 01 and 0)

patts_lst <- list()
patts_lst[[1]] <- "0"
for (i in 2:5) {
  patts_lst[[i]] <- levels(interaction(patts_lst[[i - 1]], c("0", "1")))
}

pattsC <- gsub("\\.", "", unlist(patts_lst))
vertLabs <- 1:31
# FYI children are "index*2" and "index*2+1
## For any node x -> i It's children will be 2*i and 2*i + 1 And it's parent will be floor(i/2)

i <- 4 # level of tree, 5 levels here
patts <- pattsC[nchar(pattsC) == i]
j <- 1 # length(patts)
pattsC[pattsC == "0" | stri_sub(pattsC, 1, i) == patts[j] | stri_sub(pattsC, 1, i - 1) == patts[j] |
  stri_sub(pattsC, 1, i - 2) == stri_sub(patts[j], 1, i - 2)][1:i]
j <- 2 # length(patts)
pattsC[pattsC == "0" | stri_sub(pattsC, 1, i) == patts[j] | stri_sub(pattsC, 1, i - 1) == patts[j]][1:i]


## https://stackoverflow.com/questions/6260606/adjacency-matrix-of-binary-tree-of-depth-4-in-c
## FYI children are "index*2" and "index*2+1
## For any node x -> i It's children will be 2*i and 2*i + 1 And it's parent will be floor(i/2)
library(igraph)
depth <- 5
total_nodes <- 2^depth - 1
t1 <- graph.tree(total_nodes, children = 2, mode = "out")

plot(t1, layout = layout_as_tree(t1, root = 1))
mat <- as_adjacency_matrix(t1, sparse = FALSE)
t1_dat <- data.table(nodenum = 1:total_nodes)

hrg1 <- hrg(t1, prob = rep(1, 31))

all_simple_paths(t1, from = c("16"), to = "1", mode = "all")
all_simple_paths(t1, from = c("17"), to = "1", mode = "all")

devtools::install_github("ben519/btree")
library(btree)
t2 <- make_perfect_btree(5)

# https://stackoverflow.com/questions/33768841/r-tree-with-n-branches
# https://stackoverflow.com/questions/33772031/generating-a-k-nary-tree-in-r-recursively-defined-by-a-node-wise-function
#

#### OLD STUFF BELOW

# ## The question is whether one or more of the FDR control algorithms
# ## can work at levels 2 and 3 here to control FDR
#
# bdatps1 <- melt(bdat4grps, id.vars = c("bF","gf0","gf1","gf2","gf3"),
#      measure.vars = list(c("gf0","gf1","gf2","gf3"),c("p0","p1","p2","p3")),
#      value.name = c("g","p"))
#
# pdat <- bdatps1[,.(pval=unique(p)),by=g]
# pdat[,date:=nchar(g)] ## using data for batch or level of testing
# pdat[g=="Z",date:=0]
# setnames(pdat,'g','id')
#
# p_ex <- data.frame(pval=c(.005,.01,.03,.02,.1,.06,.04),
# 		   group=c('0','0.1','0.0','0.1.0','0.1.1',
# 			   '0.0.0','0.0.1'),
#                    date=c(1,2,2,3,3,3,3))
#
# ## If indiv blocks are level 3:
#
# ### Simple FWER
# library(hommel)
#
# p_ex_b <- p_ex[4:7,'pval']
# p.adjust(p_ex_b,method="holm")
# p.adjust(p_ex_b,method="BH")
# hom.obj <- hommel(p_ex_b,simes=TRUE)
# p.adjust(hom.obj) ## hommel's adjustment
#
#
#
#
# SAFFRONstar(pdat$pval[1:3],version='batch',batch.sizes=c(1,2))
#
#
# ps <- unique(bdatps1$p)
#
# ## Perhaps each set of nested p-values can be thought of as arriving in sequence. For example here, one set would be:
# ## p0|gf0 -> {
# ## { p1|gf1==level(gf1)[1] -> { p2|gf2==level(gf2)[1] }, p1|gf1==level(gf1)[2] }
#
# ## library(DiagrammeR)
# #library(igraph)
# #tree1 <- make_tree(15,2)
# #plot(tree1)
#
#
# ## Simes, R. J. (1986). An improved Bonferroni procedure for multiple tests of significance. Biometrika, 73(3):751-754.
#
