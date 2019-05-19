# Load the library.
library(DESeq2)

# Set up the conditions based on the experimental setup.
cond_1 = rep("cond1", 3)  ##check settings
cond_2 = rep("cond2", 3)  ##check settings

# Read the data from the standard input.
countData = read.table("rcSperm", header=TRUE, sep="\t", row.names=1)
compareGroup = countData[,c(1:3,13:15)]  ##check settings
# Build the dataframe from the conditions
samples = names(compareGroup)
samples
condition = factor(c(cond_1, cond_2))
colData = data.frame(samples=samples, condition=condition)

#coldata = read.table(coldata_file, header=TRUE, sep="\t", row.names=1 )

# Create DESEq2 dataset.
dds = DESeqDataSetFromMatrix(compareGroup, colData=colData, design = ~condition)

#Set the reference to be compared
dds$condition = relevel(dds$condition,"cond1")

# Run deseq
dds = DESeq(dds)

# Format the results.
res = results(dds)

# Turn it into a dataframe to have proper column names.
res_out = data.frame('id'=rownames(res),res)

# Get normalized counts and write this to a file
nc = counts(dds,normalized=TRUE)

# Turn it into a dataframe to have proper column names.
dt = data.frame("id"=rownames(nc),nc)

#merge normalized count and DElist
total <- merge(res_out,dt,by='id')
write.table(total,file='DESeq2_out.txt',sep='\t', row.names=FALSE,quote=FALSE)

#MA plot
par(mar=c(4,4,2,2))
svg(filename = 'MAplot.svg')
x=subset(total,total$baseMean>5)
col=ifelse(x$pvalue>0.05|abs(x$log2FoldChange)<=0,'gray', ifelse(x$pvalue<=0.05&x$log2FoldChange>0,'red','blue'))
plot(x$baseMean, x$log2FoldChange, pch=16, col=col, log='x',cex=0.5, xlab="baseMean", ylab="log2FoldChange", 
     main="MAplot",ylim=c(-8,8))
abline(h=0,col='gray66')
legend("topright", "P<=0.05, log2FC>0", cex = 0.8, col = 'red', pch=16)
legend("bottomright", "P<=0.05, log2FC<0", cex = 0.8, col = 'blue', pch=16)
legend("bottomleft", "P>0.05", cex = 0.8, col = 'gray', pch=16)
dev.off()

##mirna target MA plot
miRNA_target_list<-read.table('out.txt_gene_table_intersect_2_list',header = FALSE)

svg(filename = 'miRNA targeted MAplot.svg')
x=subset(total,total$id%in%miRNA_target_list$V1&total$baseMean>5)
col=ifelse(x$pvalue>0.05|abs(x$log2FoldChange)<=0,'gray', ifelse(x$pvalue<=0.05&x$log2FoldChange>0,'red','blue'))
plot(x$baseMean, x$log2FoldChange, pch=16, col=col, log='x',cex=0.6, xlab="baseMean", ylab="log2FoldChange", 
     main="miRNA targeted MAplot",ylim=c(-3,3))
abline(h=0,col='gray66')
legend("topright", "P<=0.05, log2FC>0", cex = 0.8, col = 'red', pch=16)
legend("bottomright", "P<=0.05, log2FC<0", cex = 0.8, col = 'blue', pch=16)
legend("bottomleft", "P>0.05", cex = 0.8, col = 'gray', pch=16)
dev.off()

sig_gene_df<-subset(total,total$pvalue<0.05)  
write.table(sig_gene_df[,1],file='sig_gene_list',sep='\t', row.names=FALSE,quote = FALSE,col.names = FALSE)  #for GO term
write.table(sig_gene_df,file='sig_gene_exp',sep='\t', row.names=FALSE,quote = FALSE)  #DE sig

paste('DE genes: ', nrow(sig_gene_df))
paste('DE genes up-regulated: ',nrow(subset(sig_gene_df,sig_gene_df$log2FoldChange>0)))
sig_gene_up<-nrow(subset(sig_gene_df,sig_gene_df$log2FoldChange>0))
paste('DE genes down-regulated: ',nrow(subset(sig_gene_df,sig_gene_df$log2FoldChange<0)))
sig_gene_down<-nrow(subset(sig_gene_df,sig_gene_df$log2FoldChange<0))
paste('KO sRNA targeted genes: ',nrow(miRNA_target_list))

##heatmap
library(pheatmap)
sig_gene_df_hm<-sig_gene_df[,8:13] #sample columns  ##check settings
sig_gene_log<-log2(sig_gene_df_hm+1)
rownames(sig_gene_log)<-sig_gene_df$id
sig_gene_log <- sig_gene_log - rowMeans(sig_gene_log)
sig_gene_log[sig_gene_log > 2] <- 2
sig_gene_log[sig_gene_log < -2] <- -2
pheatmap(sig_gene_log[,], cluster_cols=FALSE, filename = 'heatmap_sig.pdf')
dev.off()

##chi-square test plot
sig_gene<-subset(total,total$pvalue<0.05)[,1]
no_sig_gene<-subset(total,total$pvalue>0.05)[,1]

target_sig<-length(which(sig_gene%in%miRNA_target_list$V1))
notarget_sig<-length(which(!(sig_gene%in%miRNA_target_list$V1)))
target_nosig<-length(which(no_sig_gene%in%miRNA_target_list$V1))
notarget_nosig<-length(which(!(no_sig_gene%in%miRNA_target_list$V1)))

TargetSigExp<-sig_gene_df[sig_gene%in%miRNA_target_list$V1,]
paste('DE genes targeted up-regulated: ',nrow(subset(TargetSigExp,TargetSigExp$log2FoldChange>0)))
target_sig_up<-nrow(subset(TargetSigExp,TargetSigExp$log2FoldChange>0))
paste('DE genes targeted down-regulated: ',nrow(subset(TargetSigExp,TargetSigExp$log2FoldChange<0)))
target_sig_down<-nrow(subset(TargetSigExp,TargetSigExp$log2FoldChange<0))

write.table(TargetSigExp[,1],file='targeted_sig_gene_list',sep='\t', row.names=FALSE,quote = FALSE,col.names = FALSE)  #for GO term
write.table(TargetSigExp,file='targeted_sig_gene_exp',sep='\t', row.names=FALSE,quote = FALSE)  #DE sig targeted

chisq_table<-data.frame(Targeted=c(target_sig,target_nosig),Not_targeted=c(notarget_sig,notarget_nosig),row.names = c('P<0.05','P>0.05'))
chisq.test(chisq_table)
chisq<-chisq.test(chisq_table)
k_pval<-chisq.test(chisq_table)[3]
k_df<-chisq.test(chisq_table)[2]
X_squared<-chisq.test(chisq_table)[1]
k_method<-chisq.test(chisq_table)[4]
chisq_table
library(corrplot)
svg(filename = 'miRNA targeted chi-square.svg')
corrplot(chisq$residuals, is.cor = FALSE,las=2, main='miRNA_target', mar=c(2,2,2,2))
dev.off()


library(VennDiagram)
svg(filename = 'venn_miRNA_targeted_gene_vs_DE_gene.svg')
draw.pairwise.venn(nrow(sig_gene_df), nrow(miRNA_target_list), chisq_table[1,1], category = c('DE_gene','miRNA_targeted_gene'), lty = rep(2, 2), fill = c("#d32536ff", "#a135cfff"), 
                   alpha = rep(0.5, 2), cat.pos = c(0, 0), cat.dist = rep(0.025, 2), cex = 2)
dev.off()

cat(' DE genes: ', nrow(sig_gene_df),'\n',
    'DE genes up-regulated: ',sig_gene_up,'\n',
    'DE genes down-regulated: ',sig_gene_down,'\n\n',
    'KO sRNA targeted genes: ',nrow(miRNA_target_list),'\n',
    'DE genes targeted up-regulated: ',target_sig_up,'\n',
    'DE genes targeted down-regulated: ',target_sig_down,'\n\n\n',
    'method:',k_method[[1]],'\n\n',
    'p-value =',k_pval[[1]],'\n',
    'df =',k_df[[1]],'\n',
    'X-squared =',X_squared[[1]],'\n\n',
    '\t', colnames(chisq_table)[1],'\t',colnames(chisq_table)[2],'\n',
    rownames(chisq_table)[1],'\t',chisq_table[1,1],'\t',chisq_table[1,2],'\n',
    rownames(chisq_table)[2],'\t',chisq_table[2,1],'\t',chisq_table[2,2],'\n',
    file='stats.txt')

