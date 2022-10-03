library(readxl)
library(stringr)

#Open and read the file
alldata <- read.delim("counts.txt", stringsAsFactors = FALSE, row.names = 1)
View(alldata)

#Create a dataframe using the columns from 6 to 11
df <- as.data.frame(alldata[,(6:11)])
View(df)

#To make it easier, we reduce the name of the columns focusing on the type of sample
new_name <- str_extract(names(df), "[:alpha:]{3,4}_[:digit:]")
colnames(df) <- new_name
View(df)

###############
#NORMALIZATION#
###############

library(limma) #package limma for normalization
library(edgeR)

#We observed genes with 0 in all the columns, we remove them
df <- df[apply(df, 1, function(x) any(x != 0)),]
View(df)

data_inf <- read.delim("data_inf.txt") #information about the samples

group <- factor(data_inf$type) #condition group

design <- model.matrix(~0+group) #matrix design for limma
colnames(design) <- c("CTRL", "STR") #rename the columns
design

norm <- calcNormFactors(df) #normalization factors
y <- voom(df, design, lib.size=colSums(df)*norm) #function voom for normalization
df_norm <- y$E #here we have our normalized data
View(df_norm)

############
#ClUSTERING#
############

library(ggplot2)
library(factoextra)

#elbow method
#Kmeans for a range of values to determine the value from which
#the reduction in the total sum of variance within the cluster isn't substantial

df_norm1 <- as.data.frame(df_norm)

#It takes a few minutes
fviz_nbclust(x = df_norm1, FUNcluster = kmeans, method = "wss", k.max = 15, 
             diss = get_dist(df_norm1, method = "euclidean"), nstart = 50)

#From 4 clusters it begins to stabilize
#For the cluster representation, we obtain a new df with the means of the gene expression for each condition

df_norm1$CTRL <- rowMeans(df_norm1[,c("CTRL_1", "CTRL_2", "CTRL_3")])
df_norm1$STR <- rowMeans(df_norm1[,c("STR_1", "STR_2", "STR_3")])
df_norm1 <- df_norm1[,(7:8)]
  
kdata <- kmeans(df_norm1,4)

ggplot(data = df_norm1, aes(x = CTRL, y = STR, col = as.factor(kdata$cluster), legend("Clusters"))) + geom_point() +
  labs(color = "Clusters") + labs(title = "KMeans Clustering")


##################################
#DIFfERENTIAL EXPRESSION ANALYSIS#
##################################

fit <- lmFit(y,design) #we fit the linear model for each gene

cont.matrix <- makeContrasts(STR-CTRL, levels = design) #contrast matrix
cont.matrix

fit <- contrasts.fit(fit, cont.matrix) #we calculate estimated coefficients and standard errors for a given set of tests

#moderated t-statistics of differential expression using empirical Bayesian moderation of standard errors
fit <- eBayes(fit)
options(digits=3)

#set pvalue adjusted threshold
mypval=0.05

colnames(fit$coefficients)
mycoef="STR - CTRL"

#output table of the 10 most important genes for this comparison
topTable(fit,coef=mycoef)

#we can export it for the report
write.table(topTable(fit,coef=mycoef),file="Genes_difference_10.txt",row.names=T,quote=F,sep="\t")

#full table ("n = number of genes in the fit")
limma.res <- topTable(fit,coef=mycoef,n=dim(fit)[1])

#only significant genes (adjusted p-value < mypval)
limma.res.pval <- topTable(fit,coef=mycoef,n=dim(fit)[1],p.val=mypval)
dim(limma.res.pval)

#Significant genes with low adjusted p-value and high FC
myfc=2

limma.res.pval.FC <- limma.res.pval[which(abs(limma.res.pval$logFC)>myfc),]
dim(limma.res.pval.FC)

#to order the genes based on their code
library(tibble)

difference_df <- limma.res.pval.FC
difference_df <- tibble::rownames_to_column(difference_df, "Gene")
difference_df <- difference_df[order(difference_df$Gene),]

require(multtest)
library(RankProd)

#Our normalized data
View(df_norm)

#To obtain the genes
geneNames <- row.names(df_norm)

#Class
#CTRL 1, STR 0
data.cl <- as.data.frame(design)
data.cl <- data.cl[,2] #now it's numeric with only the values of class

#Statistical test
#Function mt.teststat(X,classlabel,test="t",nonpara="n")
#test="t", tests are based on two-sample Welch's t-statistics (unequal variances)
#nonpara="n", original data is used

teststat<-mt.teststat(df_norm, data.cl, test="t",nonpara="n")

#Correction by Benjamini & Hochberg (BH) (FDR < 0.05)

DF=3 #degrees of freedom default
rawp0<-2*(1-pt(abs(teststat), DF))
res.multtest <- mt.rawp2adjp(rawp0, "BH")
res.multtest <-res.multtest$adjp[order(res.multtest$index),]
res.multtest <-cbind(geneNames,res.multtest)

#geneNames as row names
multtestd <- as.data.frame(res.multtest)
row.names(multtestd) <- multtestd$geneNames
multtestd <- multtestd[,-1]

length.BH.significant <-length(which(multtestd[,2]<0.05))
length.BH.significant
BH.significant <-subset(multtestd, multtestd[,2]<0.05)
dim(BH.significant)

#We join the table of previous results with this new data

BHdata <- tibble::rownames_to_column(BH.significant, "Gene")
BHdata <- BHdata[order(BHdata$Gene),]

significant_genes  <- merge(x = difference_df, y = BHdata , by = c("Gene"))

#we write the limma output table for significant genes to a tab-delimited file
write.table(significant_genes,file="Genes_difference.txt",row.names=F,quote=F,sep="\t")


#########
#HEATMAP#
#########

library("gplots")

#Our normalized data
View(df_norm)

#We create another dataframe to select the significant genes

dfgenes <- as.data.frame(df_norm)
dfgenes <- tibble::rownames_to_column(dfgenes, "Gene")
dfgenes <- dfgenes[order(dfgenes$Gene),]

genedata <- dfgenes[dfgenes$Gene %in% significant_genes$Gene,] #genes in common

#geneNames as row names
row.names(genedata) <- genedata$Gene
genedata <- genedata[,-1]

data <- as.matrix(genedata)

heatmap.2(data, scale = "row", col = bluered(100), 
          trace = "none", density.info = "none", Colv = TRUE,
          colsep=c(3:3, 6:6), cexCol = 1.2)


##############
#VOLCANO PLOT#
##############

library(tidyverse)
library(ggrepel)

#To perform the Volcano Plot, we are going to select the significant genes but
#without passing the Fold Change filter, in this case the representation is more complete
View(limma.res)

VPgenes <- tibble::rownames_to_column(limma.res, "Gene")

VPdata <- VPgenes[, c("Gene", "logFC", "P.Value", "adj.P.Val")]
colnames(VPdata) <- c("Genes", "logFC", "PValue", "FDR")

VPdata <- VPdata%>%
  mutate(
    Expression = case_when(logFC >= log(2) & FDR <= 0.05 ~ "Up-regulated",
                           logFC <= -log(2) & FDR <= 0.05 ~ "Down-regulated",
                           TRUE ~ "Unchanged")
  )

ggplot(VPdata, aes(logFC, -log(FDR,10))) +
  geom_point(aes(color = Expression), size = 2/5) +
  xlab(expression("log"[2]*"FC")) + 
  ylab(expression("-log"[10]*"FDR")) +
  scale_color_manual(values = c("green", "gray50", "red")) +
  guides(colour = guide_legend(override.aes = list(size=1.5))) 

###############
#GO ENRICHMENT#
###############

library(clusterProfiler) # GO enrichment, KEGG enrichment and plots
library(org.At.tair.db) # AGI codes and other symbols for A. thaliana genes
# keytypes(org.At.tair.db) # To see the different keys (symbol, entry, etc.) available in the database
# head(keys(org.At.tair.db, keytype = "TAIR"))

#For the enrichment, we are going to select the significant genes
View(significant_genes)


# Basic enrichment function:
GOBP <- 
  enrichGO(significant_genes$Gene, # list of genes (AGI vector - genes in TAIR format -)
           OrgDb = org.At.tair.db, # database / organism
           keyType = "TAIR", # key that you are using (symbol, entry, TAIR format, etc.)
           ont = "BP", # BP, CC or MF
           pvalueCutoff = 0.05,
           pAdjustMethod = "BY", #Bonferroni-Yuketielly
           minGSSize = 1, # min number of genes per ontology in the result
           readable = FALSE # if TRUE the program changes ATxGxxxx by names (i.e., ERF2), but if the name is not recorded you'll obtain NA
  )

barplot(GOBP)

#Data
dataBP <- data.frame(GOBP$ID, GOBP$Description, GOBP$GeneRatio, GOBP$BgRatio, GOBP$pvalue, GOBP$p.adjust, GOBP$qvalue, GOBP$geneID, GOBP$Count, GOBP@ontology)
dataBP$GOBP.geneID <- gsub("/",", ", dataBP$GOBP.geneID)
colnames(dataBP) <- c("ID","term","GeneRatio","BgRatio","pvalue","adj_pval","qvalue","genes","count", "category")

write.table(dataBP,file="GOBP.txt",sep="\t",row.names=F) #export to a file


######################
#GOPLOT VISUALIZATION#
######################

library(GOplot)

datosgenes <- significant_genes[, c("Gene", "logFC")]
colnames(datosgenes) <- c("ID", "logFC")
datosGO <- dataBP[, c("ID", "term", "adj_pval", "genes", "count", "category")]

#We can select the processes of interest to facilitate visualization
process <- c("cellular response to jasmonic acid stimulus","glycosyl compound metabolic process",
             "response to toxic substance", "regulation of signal transduction",
             "response to ethylene", "response to auxin", "amine metabolic process",
             "alpha-amino acid metabolic process", "systemic acquired resistance", "response to oomycetes")

#Data (datosGO must contain ID, term, adj_pval, genes, count and category, datosgenes must contain: ID y logFC)
circ <- circle_dat(datosGO,datosgenes)

#We select the data from the category of Biological Process
dataGO <- circ[circ$category == "BP",]
genes <- unique(dataGO[,c("genes","logFC")])

chord <- chord_dat(circ, genes, process)
head(chord) #1 is assigned to the process, 0 means not

chord <- chord_dat(data = circ, genes = genes, process = process)
GOChord(chord, space = 0.02, gene.order = 'logFC', gene.space = 0.25, gene.size = 3.5)
