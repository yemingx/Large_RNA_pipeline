exp_465<-read.table('targeted_sig_gene_exp_465', header = T)
exp_dko<-read.table('targeted_sig_gene_exp_dko', header = T)
exp_tko<-read.table('targeted_sig_gene_exp_tko', header = T)
exp_xw<-read.table('targeted_sig_gene_exp_xw', header = T)
exp_xr<-read.table('targeted_sig_gene_exp_xr', header = T)

library(gplots)
venn(list(ko465=exp_465$id,dko=exp_dko$id,tko=exp_tko$id,xw=exp_tko$id,xr=exp_xr$id))

list1<-Reduce(intersect,list(ko465=exp_465$id,dko=exp_dko$id,tko=exp_tko$id,xw=exp_tko$id,xr=exp_xr$id))
write.table(list1,'415_intersect_gene_id',quote = FALSE,col.names = FALSE,row.names = FALSE)
