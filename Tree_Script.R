library(ape)

setwd("~/Desktop/frogs")

#tree file
tree <- read.tree("Tree/rooted_tree.nh")

node <- c(
  getMRCA(tree, tip = c("DQ283452.1", "KP295607.1")), #Neobatrachia node
  getMRCA(tree, tip = c("KP295607.1", "KP295606.1")), #Ceratophrys node
  getMRCA(tree, tip = c("AY843704.1", "KF214101.1")), #Odontophrynus node
  getMRCA(tree, tip = c("AY843733.1", "JQ937192.1")), #Pleurodema node
  getMRCA(tree, tip = c("GQ366253.1", "AY326046.1"))  #Phyllomedusa node
  )

#age intervals from TimeTree
age.min <- c(114.1,
             11.3,
             13.3,
             51.1,
             20.4)

age.max <- c(191,
             18.5,
             13.3,
             51.1,
             29.8)

soft.bounds <- c(FALSE,
                 FALSE,
                 FALSE,
                 FALSE,
                 FALSE)

mycalibration <- data.frame(node, age.min, age.max, soft.bounds)

ultrametric_tree <- chronos(tree, calibration = mycalibration)

key <- read.delim("Tree/key.txt", stringsAsFactors = F)

ultrametric_tree[["tip.label"]] <- key[match(ultrametric_tree[["tip.label"]], key[['Acc']] ) , 'Species']

plot(ultrametric_tree, show.tip.label = T)

