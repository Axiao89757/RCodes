setwd("D:\\Project\\RData\\GSE11051plus")

getGse <- function (number) {
  file_name <- paste0("GSE", number, "_family.soft.gz")
  gse <- getGEO(filename = file_name)
}
getEset <- function(gse) {
  gsmlist = GSMList(gse)
  # get the probeset ordering
  probesets <- Table(GPLList(gse)[[1]])$ID
  
  # make the data matrix from the VALUE columns from each GSM
  # being careful to match the order of the probesets in the platform
  # with those in the GSMs
  data.matrix <- do.call('cbind',lapply(gsmlist,function(x){
    tab <- Table(x)
    mymatch <- match(probesets,tab$ID_REF)
    return(tab$VALUE[mymatch])}))
  data.matrix <- apply(data.matrix,2,function(x) {as.numeric(as.character(x))})
  #data.matrix <- log2(data.matrix)
  
  require(Biobase)
  # go through the necessary steps to make a compliant ExpressionSet
  rownames(data.matrix) <- probesets
  colnames(data.matrix) <- names(gsmlist)
  
  # make a data.frame of pdata
  pdata <- do.call("rbind", lapply(gsmlist, function(x) {
    tab <- Meta(x)
    chs <- as.vector(tab[["characteristics_ch1"]])
    tmp <- strsplit(chs, ': ')
    name <- c("Tilte", "Sorce name", sapply(tmp, function(x)x[1]))
    value <- c(tab[["title"]], tab[["source_name_ch1"]], sapply(tmp, function(x)x[2]))
    names(value) <- name
    
    return(value)
  }))
  pdata <- data.frame(pdata)
  pheno <- as(pdata,"AnnotatedDataFrame")
  eset <- new('ExpressionSet', exprs=data.matrix, phenoData=pheno)
  
  return(eset)
}
getEx <- function(number) {
  file_name <- paste0("GSE", number, "_series_matrix.txt.gz")
  ex <- read.table(file_name, 
                   sep = "\t",
                   quote = "", 
                   fill = T, 
                   comment.char = "!", 
                   header = T)
  colnames(ex)[1] <- "ID"
  colnames(ex)[2:ncol(ex)] <- substr(colnames(ex)[2:ncol(ex)], 3, 12)
  ex[, 1] <- sapply(ex[, 1], function(x) gsub('["]', '', x))
  return(ex)
}
getGrpMat <- function(exprMatrix, phe, gist) {
  group_list <- phe[rownames(phe) %in% colnames(exprMatrix), gist]
  group.matrix <- model.matrix(~0+factor(group_list))
  rownames(group.matrix) <- colnames(exprMatrix)
  return(group.matrix)
}
getCtrMat <- function(group.matrix, formatStr) {
  contrast.matrix <- makeContrasts(formatStr, levels = group.matrix)
  return(contrast.matrix)
}
difanl = function(exprMatrix, group.matrix, contrast.matrix){
  ##step1
  fit <- lmFit(exprMatrix, group.matrix)
  ##step2
  fit2 <- contrasts.fit(fit, contrast.matrix) 
  fit2 <- eBayes(fit2)
  ##step3
  tempOutput = topTable(fit2, coef=1, n=Inf)
  nrDEG = na.omit(tempOutput) 
  
  return(nrDEG)
}
idtrf <- function (exprMatrix, ids) {
  exprMatrix <- as.data.frame(exprMatrix)
  colnames(ids)[1] <- "ID"
  exprMatrix <- merge(exprMatrix, ids, by = "ID")
  
  exprMatrix <- aggregate(x = exprMatrix[, 2:ncol(exprMatrix)-1],by = list(exprMatrix$symbol),FUN = mean)
  
  exprMatrix <- exprMatrix[, -2]
  rownames(exprMatrix) <- exprMatrix[ ,1]
  exprMatrix <-  exprMatrix[ , -1]
  
  return(exprMatrix)
}
idtrf0 <- function (exprMatrix, ids) {
  exprMatrix <- as.data.frame(exprMatrix)
  colnames(ids)[1] <- "ID"
  #exprMatrix <- merge(exprMatrix, ids, by = "ID")
  rownames(exprMatrix) <- exprMatrix[ ,"ID"]
  exprMatrix <-  exprMatrix[ , -1]
  
  exprMatrix <- exprMatrix[rownames(exprMatrix) %in% ids$ID, ]
  tmp <- by(exprMatrix,
            ids$symbol,
            function(x) rownames(x)[which.max(abs(rowMeans(x)))])
  probes <- as.character(tmp)
  exprMatrix <- exprMatrix[rownames(exprMatrix) %in% probes, ]
  
  exprMatrix <- mutate(exprMatrix, ID = rownames(exprMatrix))#增加一列，名为probe_id，值为行名
  exprMatrix <- merge(exprMatrix, ids, by = "ID") # 合并数据
  rownames(exprMatrix) <- exprMatrix$symbol # 行名变为基因名
  exprMatrix <- exprMatrix[ , !colnames(exprMatrix) %in% c("ID", "symbol")] # 删除ID列和Symbol列
  
  return(exprMatrix)
}
idtrf1 <- function (exprMatrix, ids) {
  exprMatrix <- as.data.frame(exprMatrix)
  exprMatrix$ID = rownames(exprMatrix)
  
  exprMatrix <- merge(exprMatrix, ids, by.x = "ID", by.y = "probe_id")
  
  exprMatrix <- aggregate(x = exprMatrix[, 2:ncol(exprMatrix)-1],by = list(exprMatrix$symbol),FUN = mean)
  
  exprMatrix <- exprMatrix[, -2]
  rownames(exprMatrix) <- exprMatrix[ ,1]
  exprMatrix <-  exprMatrix[ , -1]
  
  return(exprMatrix)
}
drawVlcnMap <- function(nrDEG, VPT_PValue, VPT_logFC, VPT_PValue_text, VPT_logFC_text) {
  nrDEG$threshold = as.factor(ifelse(nrDEG$P.Value < VPT_PValue & abs(nrDEG$logFC) >= VPT_logFC, 
                                     ifelse(nrDEG$logFC> VPT_logFC ,'Up','Down'),'NoSignificance'))
  
  subNrDeg <- subset(nrDEG, nrDEG$P.Value < VPT_PValue_text & abs(nrDEG$logFC) >= VPT_logFC_text)
  
  rownames(subNrDeg) <- lapply(strsplit(rownames(subNrDeg), " /// "), function(x) {
    x[1]
  })
  
  ggplot(data = nrDEG, aes(x = logFC, y = -log10(P.Value), color=threshold,label = rownames(nrDEG))) +
    geom_point(alpha=0.8, size=3) +
    scale_color_manual(values=c("blue", "grey","red")) +
    geom_vline(xintercept=c(-VPT_logFC,VPT_logFC),lty=4,col="black",lwd=0.6) +
    geom_hline(yintercept = -log10(VPT_PValue),lty=4,col="black",lwd=0.6) +
    labs(x="log2(fold change)",y="-log10 (p-value)",title="Volcano map") +
    theme(plot.title = element_text(hjust = 0.5), legend.position="right", legend.title = element_blank()) +  
    ggrepel::geom_label_repel(
      data = subNrDeg,
      aes(label = rownames(subNrDeg)),
      size = 3,
      box.padding = unit(0.5, "lines"),
      point.padding = unit(0.8, "lines"), segment.color = "black", show.legend = FALSE )
}
drawHeatMap <- function(deg, exprMatrix, foreNum) {
  choose_gene <- head(rownames(deg), foreNum)
  choose_matrix <- exprMatrix[choose_gene, ]
  rownames(choose_matrix) <- lapply(strsplit(rownames(choose_matrix), " /// "), function(x) {
    x[1]
  })
  choose_matrix <- t(scale(t(choose_matrix)))
  pheatmap(choose_matrix, 
           color = c(colorRampPalette(colors = c("green","black"))(length(choose_matrix)/2),
                     colorRampPalette(colors = c("black","red"))(length(choose_matrix)/2)
           )
  )
} 
drawHeatMapSort <- function(deg, exprMatrix, foreNum) {
  choose_gene <- c(head(rownames(deg[order(deg$logFC,decreasing = F), ] ), foreNum),
                  head(rownames(deg[order(deg$logFC,decreasing = T), ] ), foreNum))
  choose_matrix <- exprMatrix[choose_gene, ]
  rownames(choose_matrix) <- lapply(strsplit(rownames(choose_matrix), " /// "), function(x) {
    x[1]
  })
  choose_matrix <- t(scale(t(choose_matrix)))
  pheatmap(choose_matrix, 
           color = c(colorRampPalette(colors = c("green","black"))(length(choose_matrix)/2),
                     colorRampPalette(colors = c("black","red"))(length(choose_matrix)/2)
           )
  )
} 
drawHeatMapAnot <- function(deg, exprMatrix, foreNum, anot) {
  choose_gene <- head(rownames(deg), foreNum)
  choose_matrix <- exprMatrix[choose_gene, ]
  rownames(choose_matrix) <- lapply(strsplit(rownames(choose_matrix), " /// "), function(x) {
    x[1]
  })
  choose_matrix <- t(scale(t(choose_matrix)))
  pheatmap(choose_matrix, 
           col = c(colorRampPalette(colors = c("green","black"))(length(choose_matrix)/2),
                     colorRampPalette(colors = c("black","red"))(length(choose_matrix)/2)
           ),
          annotation_col = anot
  )
} 
wipeOffPart <- function(exprMatrix, names) {
  exprMatrix <- exprMatrix[, -which(colnames(exprMatrix) %in% names)]
  return(exprMatrix)
}
wipeOffPart <- function(matrix, names, margin) {
  if (margin == 1) {
    matrix <- matrix[-which(rownames(matrix) %in% names), ]
  }
  if (margin == 2) {
    matrix <- matrix[, -which(colnames(matrix) %in% names)]
  }
  
  return(matrix)
}
choosePart <- function(matrix, names, margin) {
  if(margin == 1) {
    matrix <- matrix[which(rownames(matrix) %in% names), ]
  }
  
  if(margin == 2) {
    matrix <- matrix[, which(colnames(matrix) %in% names)]
  }
  return(matrix)
}
getEnrResTab <- function(enrResTab) {
  symbolList <- lapply(strsplit(enrResTab$geneID, '/'), function (x) {
    bitr(x, fromType = "ENTREZID", toType = "SYMBOL", OrgDb = org.Hs.eg.db)
  })
  
  enrResTab$geneSymbol <- unlist(lapply(symbolList, function(x) {
    paste0(x$SYMBOL, collapse = '/')
  }))
  
  return(enrResTab)
}
ecGetf <- function(ec) {
  symbolList <- lapply(strsplit(ec$Genes, ','), function (x) {
    bitr(x, fromType = "ENTREZID", toType = "SYMBOL", OrgDb = org.Hs.eg.db)
  })
  
  ec$Genes <- unlist(lapply(symbolList, function(x) {
    paste0(x$SYMBOL, collapse = ',')
  }))
  
  return(ec)
}

library(GEOquery)
library(hgu133plus2.db)
library(limma)
library(dplyr)
library(ggplot2)
library(pheatmap)
library(ComplexHeatmap)
library(clusterProfiler)
library(ReactomePA)
library(org.Hs.eg.db)

# 流程！！！
if (F) {
  rm(list = ls())
  gse <- getGse("110551")
  save(gse, file = "gse.Rdata")
  eset <- getEset(gse)
  save(eset, file = "eset.Rdata")
  ex.rma <- getEx("110551")
  save(ex.rma, file = "ex.Rdata")
  ids <- toTable(hgu133plus2SYMBOL)
  save(ids, file = "ids.Rdata")
  phe <- pData(eset)
  save(phe, file = "phe.Rdata")
  ex.rma <- idtrf0(ex.rma, ids)
  save(ex.rma, file = "ex_tf0ed.Rdata")
  ex.rma <- idtrf(ex.rma, ids)
  save(ex.rma, file = "ex_tfed.Rdata")
  ex.rma <- idtrf1(ex.rma, ids)
  save(ex.rma, file = "ex_tf1ed.Rdata")
}
if (T) {
  load(file = "gse.Rdata")
  load(file = "eset.Rdata")
  load(file = "ex.Rdata")
  load(file = "ids.Rdata")
  load(file = "phe.Rdata")
  load(file = "ex_tfed.Rdata")
  groupMat <- getGrpMat(ex.rma, phe, "asthma.status")
  colnames(groupMat) <- c("Asthma", "Normal")
  ctrMat <- makeContrasts("Asthma - Normal", levels = groupMat)
  nrDeg <- difanl(ex.rma, groupMat, ctrMat)
}

# 分组1：哮喘中 肥胖组 - 非肥胖组
if (T) {
  ex1 <- ex.rma[, colnames(ex.rma) %in% rownames(phe[phe$asthma.status %in% "Asthma", ])]
  groupMat1 <- getGrpMat(ex1, phe, "obesity.status")
  colnames(groupMat1) <- c("Normal", "Obese")
  ctrMat1 <- makeContrasts("Obese - Normal", levels = groupMat1)
  nrDeg1 <- difanl(ex1, groupMat1, ctrMat1)
  drawVlcnMap(nrDeg1, 0.05, 0.433, 0.005, 0.6)
  drawHeatMapAnot(nrDeg1, scale(ex1), 30, anot)
  # 筛选基因
  if (T) {
    VPT_PValue <- 0.05
    VPT_logFC <- 0.433
    deg <- nrDeg1[nrDeg1$P.Value <= VPT_PValue & abs(nrDeg1$logFC) >= VPT_logFC, ]
    deg_up <- nrDeg1[nrDeg1$P.Value <= VPT_PValue & nrDeg1$logFC >= VPT_logFC, ]
    deg_down <- nrDeg1[nrDeg1$P.Value <= VPT_PValue & nrDeg1$logFC <= -VPT_logFC, ]
    deg$type <- ifelse(deg$logFC >= VPT_logFC, "up", "down")
    deg$absLogFC <- abs(deg$logFC)
  }
  
  # 富集分析
  if (T) {
    degEtrzId <- bitr(rownames(deg), "SYMBOL", "ENTREZID", org.Hs.eg.db, drop = TRUE)
    
    BP <- enrichGO(degEtrzId$ENTREZID, org.Hs.eg.db, "ENTREZID", "BP", pvalueCutoff = 0.05, qvalueCutoff = 0.05)
    CC <- enrichGO(degEtrzId$ENTREZID, org.Hs.eg.db, "ENTREZID", "CC", pvalueCutoff = 0.05, qvalueCutoff = 0.05)
    MF <- enrichGO(degEtrzId$ENTREZID, org.Hs.eg.db, "ENTREZID", "MF", pvalueCutoff = 0.5, qvalueCutoff = 0.5)
    
    KEGG <- enrichKEGG(degEtrzId$SYMBOL, "hsa", "kegg", pvalueCutoff = 0.05)
    Rtom <- enrichPathway(degEtrzId$ENTREZID, organism = "human", pvalueCutoff = 0.5)
    dotplot(BP)
    pdf(file="./result/BP_tree.pdf",width = 15,height = 10)
    plotGOgraph(BP)
    dev.off()
    
    barplot(CC)
    barplot(MF)
    dotplot(Rtom)
    dotplot(KEGG)
    
    ## 输出富集结果表格
    if (T) {
      write.csv(BP, "./result/BP.csv")
      bpResTab <- read.csv("./result/BP.csv", header = T)
      bpResTab <- getEnrResTab(bpResTab)
      bpResTab <- bpResTab[, -1]
      write.csv(bpResTab, "./result/BP.csv")
      
      write.csv(CC, "./result/CC.csv")
      ccResTab <- read.csv("./result/CC.csv")
      ccResTab <- getEnrResTab(ccResTab)
      ccResTab <- ccResTab[, -1]
      write.csv(ccResTab, "./result/CC.csv")
      
      write.csv(MF, "./result/MF.csv")
      mfResTab <- read.csv("./result/MF.csv")
      mfResTab <- getEnrResTab(mfResTab)
      mfResTab <- mfResTab[, -1]
      write.csv(mfResTab, "./result/MF.csv")
      
      write.csv(reactomeRes, "reactomeRes.csv")
      pathwayResTab <- read.csv("reactomeRes.csv")
      pathwayResTab <- getEnrResTab(pathwayResTab)
      pathwayResTab <- pathwayResTab[, -1]
      write.csv(pathwayResTab, "pathwayResTab.csv")
    }
  }
  
  if (T) {
    writeClipboard(rownames(deg))
    write.csv(deg, "./result/deg.csv")
    keyGene <- readClipboard()
    kgBP <- enrichGO(keyGene, org.Hs.eg.db, "SYMBOL", "BP", pvalueCutoff = 0.05, qvalueCutoff = 0.05)
    kgCC <- enrichGO(keyGene, org.Hs.eg.db, "SYMBOL", "CC", pvalueCutoff = 0.05, qvalueCutoff = 0.05)
    kgMF <- enrichGO(keyGene, org.Hs.eg.db, "SYMBOL", "MF", pvalueCutoff = 0.5, qvalueCutoff = 0.5)
    dotplot(kgBP)
    barplot(kgCC)
    barplot(kgMF)
  }
}

# 热图注解
anot <- data.frame(obesity.status = ifelse(groupMat1[, 1] == 1, colnames(groupMat1)[1], colnames(groupMat1)[2]))

ex2 <- wipeOffPart(ex1, c(paste0("GSM299", c("7714", 7724, 7743, 7752, 7755, 7775, 7794, 7798, 7801, 7810, 7834))))

# 分组2: 删除离群数据
if (T) {
  groupMat2 <- getGrpMat(ex2, phe, "obesity.status")
  colnames(groupMat2) <- c("Normal", "Obese")
  ctrMat2 <- makeContrasts("Obese - Normal", levels = groupMat2)
  nrDeg2 <- difanl(ex2, groupMat2, ctrMat2)
  drawVlcnMap(nrDeg2, 0.05, 0.433, 0.005, 0.6)
  drawHeatMap(nrDeg2, scale(ex2), 100)
  # 筛选基因
  if (T) {
    VPT_PValue <- 0.05
    VPT_logFC <- 0.263
    deg2 <- nrDeg2[nrDeg2$P.Value <= VPT_PValue & abs(nrDeg2$logFC) >= VPT_logFC, ]
    deg_up2 <- nrDeg2[nrDeg2$P.Value <= VPT_PValue & nrDeg2$logFC >= VPT_logFC, ]
    deg_down2 <- nrDeg2[nrDeg2$P.Value <= VPT_PValue & nrDeg2$logFC <= -VPT_logFC, ]
    deg2$type <- ifelse(deg2$logFC >= VPT_logFC, "up", "down")
    deg2$absLogFC <- abs(deg2$logFC)
  }
  
  # 富集分析
  if (T) {
    degEtrzId2 <- bitr(rownames(deg2), "SYMBOL", "ENTREZID", org.Hs.eg.db, drop = TRUE)
    
    BP2 <- enrichGO(degEtrzId2$SYMBOL, org.Hs.eg.db, "SYMBOL", "BP", pvalueCutoff = 0.05, qvalueCutoff = 0.05)
    CC2 <- enrichGO(degEtrzId2$SYMBOL, org.Hs.eg.db, "SYMBOL", "CC", pvalueCutoff = 0.5, qvalueCutoff = 0.5)
    MF2 <- enrichGO(degEtrzId2$ENTREZID, org.Hs.eg.db, "ENTREZID", "MF", pvalueCutoff = 0.5, qvalueCutoff = 0.5)
    
    KEGG2 <- enrichKEGG(degEtrzId2$ENTREZID, "hsa", "kegg", pvalueCutoff = 0.05)
    Rtom2 <- enrichPathway(degEtrzId2$ENTREZID, organism = "human", pvalueCutoff = 0.5)
    dotplot(BP2)
    pdf(file="./result/BP2_tree.pdf",width = 15,height = 10)
    plotGOgraph(BP2)
    dev.off()
    
    dotplot(CC2)
    dotplot(MF2)
    dotplot(Rtom2)
    dotplot(KEGG2)
  }
  
  if (T) {
    writeClipboard(rownames(deg))
    write.csv(deg, "./result/deg.csv")
  }
}

getEC <- function(enrichResult_BP, enrichResult_CC, enrichResult_MF) {
  library(stringr)
  
  EC_BP <- enrichResult_BP[, c(1, 2, 8, 6)]
  EC_BP$geneID <- str_replace_all(EC_BP$geneID, "/", ",")
  names(EC_BP) <- c("ID", "Term", "Genes","adj_pval")
  EC_BP$Category <- "BP"
  
  EC_CC <- enrichResult_CC[, c(1, 2, 8, 6)]
  EC_CC$geneID <- str_replace_all(EC_CC$geneID, "/", ",")
  names(EC_CC) <- c("ID", "Term", "Genes","adj_pval")
  EC_CC$Category <- "CC"
  
  EC_MF <- enrichResult_MF[, c(1, 2, 8, 6)]
  EC_MF$geneID <- str_replace_all(EC_MF$geneID, "/", ",")
  names(EC_MF) <- c("ID", "Term", "Genes","adj_pval")
  EC_MF$Category <- "MF"
  
  EC <- rbind(EC_BP, EC_CC, EC_MF)
  
  return(EC)
}
getECSingle <- function(enrichResult, Category) {
  library(stringr)
  
  EC <- enrichResult[, c(1, 2, 8, 6)]
  EC$geneID <- str_replace_all(EC$geneID, "/", ",")
  names(EC) <- c("ID", "Term", "Genes","adj_pval")
  EC$Category <- Category
  return(EC)
}
getGe <- function(deg, geneName) {
  gene <- subset.data.frame(deg, rownames(deg) %in% geneName, select = "logFC")
  gene$ID <- rownames(gene)
  
  return(gene)
}
# 可视化GO结果
if(T) {
  library(GOplot)
  EC1 <- getEC(BP, CC, MF)
  EC1 <- ecGetf(EC1)
  gene <- getGe(deg, rownames(deg))
  circ <- circle_dat(EC1, gene)
  
  EC2 <- getEC(kgBP, kgCC, kgMF)
  gene2 <- getGe(deg, keyGene)
  circ2 <- circle_dat(EC2, gene2)
  
  #### GOBar
  GOBar(circ2)
  GOBar(circ2, display = 'multiple')
  GOBar(circ2, display = 'multiple', title = 'Z-score coloured barplot', zsc.col = c('yellow', 'black', 'cyan'))
  #### GOBubble
  GOBubble(reduced_circ1, labels = 2)
  # Add a title, change the colour of the circles, facet the plot according to the categories and change the label threshold
  GOBubble(reduced_circ1, title = 'Bubble plot', colour = c('orange', 'darkred', 'gold'), display = 'multiple', labels = 2)
  # Colour the background according to the category
  GOBubble(reduced_circ1, title = 'Bubble plot with background colour', display = 'multiple', bg.col = T, labels = 3)
  # Reduce redundant terms with a gene overlap >= 0.75...
  reduced_circ <- reduce_overlap(circ, overlap = 0.95)
  # ...and plot it
  GOBubble(reduced_circ, labels = 1.3)
  #### GOCircle
  GOCircle(reduced_circ1)
  #### GOChord
  chord <- chord_dat(circ2,process = unique(circ2$term)[1:8])
  # Generate the matrix with a list of selected genes
  chord <- chord_dat(data = reduced_circ1, genes = gene)
  # Generate the matrix with selected processes
  chord <- chord_dat(data = circ, process = EC$process)
  # Create the plot
  GOChord(chord, space = 0.02, gene.order = 'logFC', gene.space = 0.25, gene.size = 4)
  
  ################
  circBP <- reduce_overlap(subset(circ, circ$category == "BP"), overlap = 1)
  circCC <- reduce_overlap(subset(circ, circ$category == "CC"), overlap = 1)
  circMF <- reduce_overlap(subset(circ, circ$category == "MF"), overlap = 1)
  reduced_circ1 <- rbind(circBP, circCC, circMF)
  reduced_circ1 <- reduced_circ1[!(reduced_circ1$ID %in% c("GO:0031838", "GO:0034774")), ]
}




