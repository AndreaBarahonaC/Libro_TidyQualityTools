setClass(Class = "facDesign", representation = representation(name = "character", factors = "list",
                                                              cube = "data.frame", star = "data.frame", centerCube = "data.frame", centerStar = "data.frame",
                                                              generator = "ANY", response = "data.frame", block = "data.frame", blockGen = "data.frame", runOrder = "data.frame",
                                                              standardOrder = "data.frame", desireVal = "list", desirability = "list", fits = "list"))



function (k = 3, p = 0, replicates = 1, blocks = 1, centerCube = 0)
{
  frameOut = fracDesign(k = k, p = p, gen = NULL, replicates = replicates,
                        blocks = blocks, centerCube = centerCube)
  return(frameOut)
}
