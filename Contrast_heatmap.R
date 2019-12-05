#---------------------------------------------------------------------------------------------
# <<< Running contrast heatmaps >>>
#---------------------------------------------------------------------------------------------
# Gene expression for your conditions
# centered and scaled, log2 transformed, covariate corrected TPM:
AA <- AA.expr.correct
DA <- DA.expr.correct

# module genes from the diffcoEx output
mod.list=pw5$module_list[c(1,2)]

outpath="~/Dropbox/Glioma_project2/TCGA_AA_DA/IDH1_ASTROCYTOMAS/DIFFCOEX_R/"

pdf("module_expr.pdf")
library(WGCNA)
CoExH(cluster_list = mod.list,cases_expr = DA, controls_expr = AA, runName = 'G2', outName = outpath)
dev.off()

#---------------------------------------------------------------------------------------------
# The CoExH function (credit: Prashant Srivastava)
#---------------------------------------------------------------------------------------------
####
# cluster_list is a list of character vectors of geneid. Names will be used to name the pdf
# cases_expr is a matrix of expression values with samples in columns and geneid in rows
# control_expr is used (when you want to compare 2 conditions on the same heatmap with one condition in each triangle
#              is used as a reference to order genes depending on their similarity of co-expression
#              is also a matrix of expression values with samples in columns and geneid in rows
# outName is the path where the pdf will be
# runName
CoExH = function(cluster_list,cases_expr,controls_expr=NULL,runName="",outName=getwd()){
  for(i in 1:length(cluster_list))
  {
    tmp=intersect(cluster_list[[i]],rownames(cases_expr))
    Mname=names(cluster_list)[i]
    casex=t(cases_expr[tmp,])
    cases.cor=bicor(casex)
    if(length(controls_expr)!=0)
    {
      ctlex=t(controls_expr[tmp,])
      controls.cor=bicor(ctlex)
      t.dist=as.dist(1-controls.cor);
      t.hclust=hclust(t.dist,method='ward.D2');
      t.cutree=cutreeHybrid(t.hclust,as.matrix(t.dist),minClusterSize=2,deepSplit=0)
      plotC1C2Heatmap(t.cutree$labels,cases.cor,controls.cor,casex,ctlex,
                      ordering=NULL,Dname=runName,outName=outName,
                      Mname=Mname)
      print(paste(Mname,"cases vs controls in",runName))
    }else{
      t.dist=as.dist(1-cases.cor);
      t.hclust=hclust(t.dist,method='ward.D2');
      t.cutree=cutreeHybrid(t.hclust,as.matrix(t.dist),minClusterSize=2,deepSplit=0)
      plotC1C2Heatmap(t.cutree$labels,cases.cor,cases.cor,casex,casex,
                      ordering=NULL,Dname=runName,outName=outName,
                      Mname=Mname)
      print(paste(Mname,"in",runName))
    }
  }
}

########### inside functions for CoExH

# datRef is the expression dataset used to order the group of genes clustered depending on the similarity of their co-expression
# colorh1 is the assignement of each gene in one of the cluster of genes
getEigenGeneValues<-function(datRef,colorh1,datAll)
{
  eigenGenesCoef<-list()
  for (c in 1:length(unique(colorh1)))
  {
    eigenGenesCoef[[unique(colorh1)[c]]]<-prcomp(scale(datRef[,which(colorh1 == unique(colorh1)[c])]))$rotation[,1]
  }
  values<-NULL
  for( c in unique(colorh1))
  {
    v<-rbind(datAll)[,which(colorh1 == c)] %*%  eigenGenesCoef[[c]]
    values<-cbind(values,sign(mean(v))*v)
  }
  colnames(values)<-unique(colorh1)
  values
}


plotC1C2Heatmap<-function(colorh1C1C2,AdjMat1C1,AdjMat1C2, datC1, datC2,ordering=NULL,Mname="",Dname="",outName=getwd())
{
  if (is.null(ordering))
  {
    h<-hclust(as.dist(1-abs(bicor(getEigenGeneValues(datC2,colorh1C1C2,rbind(datC1,datC2))))))
    for (c in h$label[h$order])
    {
      ordering<-c(ordering,which(colorh1C1C2 ==c))
    }
  }
  mat_tmp<-(AdjMat1C1[ordering,ordering])
  mat_tmp[which(row(mat_tmp)>col(mat_tmp))]<-(AdjMat1C2[ordering,ordering][which(row(mat_tmp)>col(mat_tmp))])
  diag(mat_tmp)<-0

# >>>>>> COLORS!
   colR <- colorRampPalette(c("midnightblue", "white", "red3"))
   #colR <- colorRampPalette(c("dodgerblue4","white","darkred"))

# NB!!!!! >>> If you want the plot to go staight to the folder, unsilence the pdf!!!

  #pdf(file= paste0(outName,"/",Mname,"_HBicor_",Dname,".pdf"), height=24, width = 24)
  #corrplot(mat_tmp, col=colR(200), tl.pos = "n", method = 'color',title=Mname,diag=F,shade.lwd = 0.001)
  image(mat_tmp,col=colR(200),axes=F,asp=1,main=Mname, zlim=c(-1,1))
 # dev.off()
  # cat('\t',min(mat_tmp),max(mat_tmp),'\n')
}
