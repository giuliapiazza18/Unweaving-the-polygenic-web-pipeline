# Function to unfade edges between PGS and items and modify legend for ALSPAC/TEDS project. 
# IMPORTANT: this function assumes the PGS is the last item in the network.

# network = network graph: usually found in the results of estimatenetwork $ graph. Matrix of edge weights. 
# PGSnode = character indicating the abbreviation for the polygenic score, e.g. "ANX" or "EA"
# PGSdesc = character description of the polygenic score to use in the groups to decide colours
# PGScolor = character description of the color of polygenic score node to use in the graph 
# ("grey" (anx) "#9CB469"(ea) "#56BD96"(dep) "#46BAC8"(adhd) "#99A9E2" (bmi))
# avglay = layour input.

# Giulia Oct 2022

  # NETWORK ====
  # create regular network
  unfade.nolegend <- function (network, PGSnode, PGSdesc, PGScolor, avglay, curve.matrix) {
    graph1 <- qgraph(
      network,
      layout = "spring",
      theme = "colorblind",
      fade = FALSE,
      nodeNames = c(quest.sub$abbreviation, PGSnode),
      groups = c(quest.sub$subscale, PGSnode),
      labels = c(quest.sub$subscale.node, PGSnode),
      DoNotPlot = TRUE
    )
    
    # create temporary network matrix where edges are unfaded (weights below abs(0.1) are rounded up)
    temp <- as.data.frame(network)
    
    i <- nrow(temp)
    
    index_pos <- which(temp[i] >0 & temp[i] < 0.1)
    index_neg <- which(temp[i] < 0 & temp[i] > -0.1)
    
    temp[i,index_pos] <- temp[index_pos,i] <- 0.11
    temp[i,index_neg] <- temp[index_neg,i] <- -0.11
    
    # create a matrix of colour grey
    edgecolor = ifelse(abs(temp) > 0, "grey90", NA)
    edgecolor[i,] <- NA ; edgecolor[,i] <- NA
    
    # create line type vector (dotted for negative, solid for positive, solid for all PGS connections)
    edgeline <- graph1$graphAttributes$Edges$color
    index_line_neg <- which(edgeline == "red" | edgeline == "#BF0000FF")
    edgeline[index_line_neg] <- "dotted"
    index_line_pos <- which(edgeline == "darkblue" | edgeline == "#0000D5FF")
    edgeline[index_line_pos] <- "solid"
    
    edgeline[(length(edgeline) - (length(index_neg) + length(index_pos))):length(edgeline)] <- "solid"
    n_edges = sum(abs(network[,i])>0)
    set.seed(1)
    edgeSort = c(sample(1:(length(edgeline)-n_edges)), 
                 (length(edgeline)- n_edges):length(edgeline))
    
    
    # create unfaded graph with same layout as graph1
    unfaded.graph <- qgraph.alt(temp, 
                                layout = avglay, 
                                labels = c(quest.sub$subscale.node, PGSnode),
                                legend.cex = 0.3, 
                                theme = "colorblind",
                                nodeNames = c(quest.sub$abbreviation, "Polygenic score for ", PGSnode),
                                groups = c(quest.sub$subscale, PGSdesc),
                                vsize = c(rep(5,36), 7),
                                border.width = c(rep(1,36), 2),
                                label.cex = 1.2,
                                minimum = 0.1,
                                lty = edgeline,
                                edge.color = edgecolor,
                                curveAll = T,
                                curve = curve.matrix,
                                mar = c(3,2,3,2),
                                GLratio =2.5,
                                legend = FALSE,
                                edgeSort = edgeSort,
                                color = c("#E69F00", "#56B4E9","#009E73", "#F0E442", "#CC79A7", PGScolor,  "#0072B2"),
                                normalize = TRUE)
    
  }
  
  
  