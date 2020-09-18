# 有用代码块

# 改名
if (T) {
	phe$sample <- rownames(phe)
	temp <- phe[phe$asthma.status %in% "Asthma", c(4,11)]
	temp$sample_name <- ifelse(temp$obesity.status == "Non-Obese", paste("Normal", 1:39), paste("Obese", 1:39))
	colnames(choose_ex) <- temp$sample_name  
}
# 取子集分组
if (T) {
	choose_name <- rownames(phe[phe$asthma.status %in% "Asthma", ])
	choose_ex <- ex[, colnames(ex) %in% choose_name]
}

# 筛选基因
if (T) {
	VPT_PValue <- 0.05
	VPT_logFC <- 0.433
	deg <- nrDEG[nrDEG$P.Value <= VPT_PValue & abs(nrDEG$logFC) >= VPT_logFC, ]
	deg_up <- nrDEG[nrDEG$P.Value <= VPT_PValue & nrDEG$logFC >= VPT_logFC, ]
	deg_down <- nrDEG[nrDEG$P.Value <= VPT_PValue & nrDEG$logFC <= -VPT_logFC, ]
	deg$type <- ifelse(deg$logFC >= VPT_logFC, "up", "down")
	deg$absLogFC <- abs(deg$logFC)
}
write.csv(nrDEG, file = "nrDEG.csv")
write.csv(deg, file = "deg0.433.csv")
writeClipboard(rownames(deg))

# 流程！！！
if (T) {
	gse <- getgse("110551")
	ids <- getids(gse)
	#ids=toTable(hgu133plus2SYMBOL)
	eset <- geteset(gse)
	ex <- exprs(eset)
	phe <- pData(eset)
	group_list <- phe$obesity.status
}