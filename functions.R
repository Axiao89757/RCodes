# 有用函数

#####################################################
####################### 获取 ########################
#####################################################
# 1. 获取 GSE 对象
getGse <- function (number) {
	file_name <- paste0("GSE", number, "_family.soft.gz")
	gse <- getGEO(filename = file_name)
}
# 使用方法：
##	去ncbi官网下载soft文件到工作目录，用编号字符串作为输入参数,
##	即可得到GSE对象

# 2. 获取 ids
getIds <- function (gse) {
	# 获取GPL
	gpl <- GPLList(gse)[[1]]
	# 获取ids
	ids <- Table(gpl)[, c("ID", "Gene Symbol")]
	# 去除空值
	ids <- ids[!ids$`Gene Symbol` == '', ]

	return(ids)
}
# 说明：
##	输入：GSE对象
##  输出：ids data.frame
##	有可能列名不兼容，需要手动修改

# 3. 获取表达矩阵对象
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
# 说明：
##	输入：GSE对象
##	输出：表达矩阵对象
##	可通过 exprs(eset) 和 pData(eset)函数来分别获取表达矩阵和样本说明矩阵

# 4. 获取表达矩阵
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
# 说明：
##	输入：编号字符串
##	输出：表达矩阵
##	去ncbi官网下载表达矩阵文件到工作目录，用编号字符串作为输入参数,即可得到表达矩阵

# 5. 获取分组矩阵
getGrpMat <- function(exprMatrix, phe, gist) {
	group_list <- phe[rownames(phe) %in% colnames(choose_name), gist]
	group.matrix <- model.matrix(~0+factor(group_list))
	colnames(group.matrix) <- c("normal", "obese")
	rownames(group.matrix) <- colnames(exprMatrix)
	return(group.matrix)
}
# 说明：
##	输入：
##		exprMatrix：data.frame 表达矩阵，列名是GSE编号，行名是基因名
##		phe：data.frame 样本的说明矩阵
##		gist：character 分组依据，phe的一个column
##	输出：data.frame 分组矩阵

# 6. 获取比较矩阵
getCtrMat <- function(group.matrix, formatStr) {
	contrast.matrix <- makeContrasts(formatStr, levels = group.matrix)
	return(contrast.matrix)
}
# 说明：
##	输入：
##		group.matrix：data.frame 分组矩阵
##		formatStr：character 比较的格式，比如 "obese - normal"
##	输出：matrix array 比较矩阵


######################################################
###################### 分析 ##########################
######################################################
# 1. 差异分析
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
# 说明：
##	输入：
##		exprMatrix：data.frame 表达矩阵，列名是GSE编号，行名是基因名
##		group.matrix：分组矩阵
##		contrast.matrix：比较矩阵，注意 谁比谁，对照组在后面
##	输出：nrDEG

# 2. 富集分析
# 载入包
erianl = function(deg) {
	library(clusterProfiler)
	library(org.Hs.eg.db)
	library(ReactomePA)
	
	# id转换
	DEG.entrez_id = mapIds(x = org.Hs.eg.db,
						   keys = rownames(deg),
						   keytype = "SYMBOL",
						   column = "ENTREZID")
	DEG.entrez_id <- na.omit(DEG.entrez_id)
	
	# BP
	erich.go.BP = enrichGO(gene = DEG.entrez_id,
                       OrgDb = org.Hs.eg.db,
                       keyType = "ENTREZID",
                       ont = "BP",
                       pvalueCutoff = 0.05,
                       qvalueCutoff = 0.05)
	##分析完成后，作图
	dotplot(erich.go.BP)
	
	# CC 
	erich.go.CC = enrichGO(gene = DEG.entrez_id,
                       OrgDb = org.Hs.eg.db,
                       keyType = "ENTREZID",
                       ont = "CC",
                       pvalueCutoff = 0.05,
                       qvalueCutoff = 0.05)
	## 画图
	barplot(erich.go.CC)
	
	# MF
	erich.go.MF = enrichGO(gene = DEG.entrez_id,
                       OrgDb = org.Hs.eg.db,
                       keyType = "ENTREZID",
                       ont = "MF",
                       pvalueCutoff = 0.5,
                       qvalueCutoff = 0.5)
	## 画图
	barplot(erich.go.MF)
	
	# KEGG
	erich.KEGG = enrichKEGG(gene = DEG.entrez_id,
                        organism  = "hsa",
                        pvalueCutoff  = 5)
	dotplot(erich.KEGG)	# 画气泡图
	
	# Reactome
	reactomeRes <- enrichPathway(gene = DEG.entrez_id, pvalueCutoff = 0.05, qvalueCutoff = 0.05, readable = FALSE)
	
	
	getEnrResTab <- function(enrResTab) {
		symbolList <- lapply(strsplit(enrResTab$geneID, '/'), function (x) {
			bitr(x, fromType = "ENTREZID", toType = "SYMBOL", OrgDb = org.Hs.eg.db)
		})

		enrResTab$geneSymbol <- unlist(lapply(symbolList, function(x) {
			paste0(x$SYMBOL, collapse = '/')
		}))

		return(enrResTab)
	}
	write.csv(erich.go.BP, "erich.go.BP.csv")
	bpResTab <- read.csv("erich.go.BP.csv", header = T)
	bpResTab <- getEnrResTab(bpResTab)
	bpResTab <- bpResTab[, -1]
	write.csv(bpResTab, "bpResTab.csv")

	write.csv(erich.go.CC, "erich.go.CC.csv")
	ccResTab <- read.csv("erich.go.CC.csv")
	ccResTab <- getEnrResTab(ccResTab)
	ccResTab <- ccResTab[, -1]
	write.csv(ccResTab, "ccResTab.csv")

	write.csv(erich.go.MF, "erich.go.MF.csv")
	mfResTab <- read.csv("erich.go.MF.csv")
	mfResTab <- getEnrResTab(mfResTab)
	mfResTab <- mfResTab[, -1]
	write.csv(mfResTab, "mfResTab.csv")

	write.csv(reactomeRes, "reactomeRes.csv")
	pathwayResTab <- read.csv("reactomeRes.csv")
	pathwayResTab <- getEnrResTab(pathwayResTab)
	pathwayResTab <- pathwayResTab[, -1]
	write.csv(pathwayResTab, "pathwayResTab.csv")
	
}





##########################################################
###################### 数据转换 ##########################
##########################################################
# id转换
idtrf <- function (exprMatrix, ids) {
	exprMatrix <- as.data.frame(exprMatrix)
	exprMatrix$ID <- rownames(exprMatrix)
	colnames(ids)[1] <- "ID"
	exprMatrix <- merge(exprMatrix, ids, by = "ID")

	exprMatrix <- aggregate(x = exprMatrix[, 2:ncol(exprMatrix)-1],by = list(ex$symbol),FUN = mean)

	exprMatrix <- exprMatrix[, -2]
	rownames(exprMatrix) <- exprMatrix[ ,1]
	exprMatrix <-  exprMatrix[ , -1]
	
	return(exprMatrix)
}
# 说明：
##	输入：
##		exprMatrix：表达矩阵，列名是GSE编号，行名是探针
##		ids：基因id和symbol对应表
##	输出：id转换后的表达矩阵，列名是GSE编号，行名是基因名


##########################################################
###################### 绘图 ##############################
##########################################################
# 1. 画火山图
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
# 说明：
##	输入：
##		nrDEG：data.frame 差异分析结果
##		VPT_PValue：numeric p值阈值
##		VPT_logFC：numeric logFC值阈值
##		VPT_PValue_text：numeric p值标注点阈值
##		VPT_logFC_text：numeric logFC值标注点阈值

# 2. 画热图
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
# 说明：
##	输入：
##		deg：data.frame 差异分析结果，或者差异基因
##		exprMatrix：data.frame 表达矩阵，列名是GSE编号，行名是探针
##		foreNum: numeric 选择deg的前foreNum个基因















