library(dplyr)
library(tidyr)
library(tibble)
library(reshape2)
library(ggplot2)
library(tidybayes)
library(car)

ime.result <- read.csv(file.choose())
ddl <- ime.result %>% 
  as.data.frame() %>% 
  pivot_longer(cols=2:23,
               names_to= "celltype",
               values_to = "Proportion")

if (T) {
  library(ggplot2)
  colTab <- rainbow(length(colTab))
  names(colTab) <- as.character(levels(as.factor(ddl$celltype)))
  ggplot(ddl, aes(Mixture, Proportion)) + 
    geom_bar(position = "stack",stat = "identity", aes(fill = celltype, color = celltype))+
    theme_bw()+
    guides(fill=guide_legend(ncol=1)) +
    coord_flip() +
    scale_fill_manual(values = colTab)
}

###################### 鎻愮惔鍥? ############################
# 婊ら櫎澶ч儴鍒嗕负闆剁殑鏁版嵁
if (T) {
  choose.names <- character(0)
  ime.sum <- apply(ime.result[, 2:23], 2, sum)
  for (i in 1:length(ime.sum)) {
    if (ime.sum[i] > 1) {
      print(names(ime.sum[i]))
      choose.names <- c(choose.names, names(ime.sum[i]))
    }
  }
  choose.ime <- subset(ime.result, select = c("Mixture", choose.names))
}

# 鏁版嵁鍒嗕负涓ょ粍
if (T){
  grp.obese <- subset(phe, asthma.status == "Asthma" & obesity.status == "Obese",select = c(4, 5)) %>%
    rownames()
  grp.normal <- subset(phe, asthma.status == "Asthma" & obesity.status == "Non-Obese",select = c(4, 5)) %>%
    rownames()
  choose.ime$Group <- ifelse (choose.ime$Mixture %in% grp.obese, "Obese", "Non-Obese")
  
  vl.data <- choose.ime[, 2:length(choose.ime)] %>%
    melt(variable.name = "celltype")
  
}

# 鐢绘彁鐞村浘
if (T) {
  ggplot(vl.data, aes(x = celltype, y = value, fill = Group)) +
    stat_halfeye(data = vl.data[vl.data$Group == "Non-Obese",], show_interval = FALSE, trim=FALSE,side = "left", alpha = .5) +
    stat_halfeye(data = vl.data[vl.data$Group == "Obese",], show_interval = FALSE, trim=FALSE,side = "right", alpha = .5) +
    geom_boxplot( width =.25, outlier.alpha = 0)+
    coord_flip() +
    theme_bw()  
}

# t检验分组差异
if (T) {
  g <- split.data.frame(vl.data, vl.data$celltype)
  lapply(g, function(x) {
    
    summary(as.data.frame(x))
    
    
    #summary(aov(data = as.data.frame(x), formula = value~Group))
  })
}





temp <- split.data.frame(x, x$Group)
xx <- as.vector(temp[[1]]$value)
yy <- as.vector(temp[[2]]$value)
xy <- c(xx, yy)
a <- factor(c(rep(1, length(xx)), rep(2, length(yy))))
leveneTest(xy~a)
t.test(xx, yy)
























