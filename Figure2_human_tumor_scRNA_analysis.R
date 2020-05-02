
library(Matrix)

rm(list=ls())

#################
# # to download & load human counts and metadata
wd <- "~/Downloads/analysis_dir"

if(!dir.exists(wd)){
  dir.create(wd)
}

data_url <- "https://www.dropbox.com/s/i6907pbef4hqj1j/human_tumor_scdata.rd?dl=1"
if(Sys.info()["sysname"]=="Windows"){
  download.file(url = data_url,destfile = file.path(wd,"human_tumor_scdata.rd"),mode="wb")
}else{
  download.file(url = data_url,destfile = file.path(wd,"human_tumor_scdata.rd"))
}
download_check <- try(load(file.path(wd,"human_tumor_scdata.rd")))
if(download_check!="human_tumor_scdata"){stop("Error using download.file(), due to operating system?-- please manually download human DC dataset from https://www.dropbox.com/s/i6907pbef4hqj1j/human_tumor_scdata.rd?dl=1 and place in working directory and remove this chunk of code")}
###############

#load("/users/andrew leader/Dropbox/PAPER_DATA_REPO/Oh_et_al_Nature_Cancer_2020_files/human_tumor_scdata.rd")

message("Summing celltype expression over patients")
patient_counts <- list()
for(patient in unique(human_tumor_scdata$cell_to_patient)){
  message(patient)
  patient_mask <- human_tumor_scdata$cell_to_patient==patient
  s <- split(names(human_tumor_scdata$cell_to_patient)[patient_mask],
             human_tumor_scdata$cell_to_sub_lineage[patient_mask])
  s <- s[names(s)!=""]
  patient_counts[[patient]] <- do.call(cbind,lapply(s,function(x){rowSums(human_tumor_scdata$umitab[,x,drop=F])}))
}


PDL1_count_mat <- matrix(0,
                         nrow=length(patient_counts),
                         ncol=length(setdiff(unique(human_tumor_scdata$cell_to_sub_lineage),"")),
                         dimnames=list(names(patient_counts),
                                       setdiff(unique(human_tumor_scdata$cell_to_sub_lineage),"")))

CD80_count_mat <- PDL1_count_mat

for(patient in names(patient_counts)){
  PDL1_count_mat[patient,colnames(patient_counts[[patient]])] <- patient_counts[[patient]]["CD274",,drop=F]
  CD80_count_mat[patient,colnames(patient_counts[[patient]])] <- patient_counts[[patient]]["CD80",,drop=F]
}

PDL1_percent <- PDL1_count_mat/rowSums(PDL1_count_mat)*100
CD80_percent <- CD80_count_mat/rowSums(CD80_count_mat)*100
PDL1_percent[PDL1_percent < 1] <- 1
CD80_percent[CD80_percent < 1] <- 1

PDL1_expr_mat <- matrix(NA,
                                          nrow=length(patient_counts),
                                          ncol=length(setdiff(unique(human_tumor_scdata$cell_to_sub_lineage),"")),
                                          dimnames=list(names(patient_counts),
                                                        setdiff(unique(human_tumor_scdata$cell_to_sub_lineage),"")))
CD80_expr_mat <- PDL1_expr_mat

for(patient in names(patient_counts)){
  numi <- colSums(patient_counts[[patient]])
  PDL1_expr_mat[patient,colnames(patient_counts[[patient]])] <- patient_counts[[patient]]["CD274",,drop=F]/numi
  CD80_expr_mat[patient,colnames(patient_counts[[patient]])] <- patient_counts[[patient]]["CD80",,drop=F]/numi
}
PDL1_expr_mat[PDL1_expr_mat<1e-6] <- 1e-6
CD80_expr_mat[CD80_expr_mat<1e-6] <- 1e-6


tab <- table(human_tumor_scdata$cell_to_patient,human_tumor_scdata$cell_to_sub_lineage)

for(sub_lin in colnames(PDL1_expr_mat)){
  PDL1_expr_mat[tab[rownames(PDL1_expr_mat),sub_lin]<20,sub_lin] <- NA
  CD80_expr_mat[tab[rownames(CD80_expr_mat),sub_lin]<20,sub_lin] <- NA
}



#current
annots_ord <- c("AM","IM","MoMac","CD14 Mono","CD16 Mono","moDC","DC2","Mature DC","DC1")

plot_open <-function(fn){
  pdf(fn,height=4,width=4)
  par(mgp=c(2,1,0))
}

plot_open(fn=file.path(wd,"PDL1_fraction.pdf"))
boxplot.matrix(log10(PDL1_percent[,annots_ord]),
               boxwex=.5,las=2,range=0,frame=FALSE,border="red",ylim=c(0,2),xaxt="n",yaxt="n",pars=list(lty=1))
mtext(expression(atop("Fraction of "*italic("CD274")*" RNA contributed"),"by cell population"))

mtext(expression("% contribution of "*italic("CD274")*" RNA)"),side=2,line=2.5)
axis(side=2,at=0:2,labels=c("< 1",10,100),las=2)
mtext(annots_ord,side=1,las=2,at=seq(length(annots_ord)))
dev.off()

plot_open(fn=file.path(wd,"CD80_fraction.pdf"))
boxplot.matrix(log10(CD80_percent[,annots_ord]),
               boxwex=.5,las=2,range=0,frame=FALSE,border="red",ylim=c(0,2),xaxt="n",yaxt="n",pars=list(lty=1))
mtext(expression(atop("Fraction of "*italic("CD80")*" RNA contributed"),"by cell population"))

mtext(expression("% contribution of "*italic("CD80")*" RNA)"),side=2,line=2.5)
axis(side=2,at=0:2,labels=c("< 1",10,100),las=2)
mtext(annots_ord,side=1,las=2,at=seq(length(annots_ord)))
dev.off()

plot_open(fn=file.path(wd,"PDL1_exprs.pdf"))
boxplot.matrix(log10(PDL1_expr_mat[,annots_ord]),
               boxwex=.5,las=2,range=0,frame=FALSE,border="red",ylim=c(-6,-3),xaxt="n",yaxt="n",pars=list(lty=1))
mtext(expression("Log"["10"]*"("*italic("CD274")*" counts/total RNA)"),side=2,line=2.5)
axis(side=2,at=c(-6:-3),labels=c("< -6",-5:-3),las=2)
mtext(expression(italic("CD274")*" expression per cell type"))
mtext(annots_ord,side=1,las=2,at=seq(length(annots_ord)))
dev.off()


plot_open(fn=file.path(wd,"CD80_exprs.pdf"))
boxplot.matrix(log10(CD80_expr_mat[,annots_ord]),
               boxwex=.5,las=2,range=0,frame=FALSE,border="red",ylim=c(-6,-3),xaxt="n",yaxt="n",pars=list(lty=1))
mtext(expression("Log"["10"]*"("*italic("CD80")*" counts/total RNA)"),side=2,line=2.5)
axis(side=2,at=c(-6:-3),labels=c("< -6",-5:-3),las=2)
mtext(expression(italic("CD80")*" expression per cell type"))
mtext(annots_ord,side=1,las=2,at=seq(length(annots_ord)))
dev.off()

