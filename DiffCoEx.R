# DiffCoEx
options(stringsAsFactors=F)

# Expression matrix of condition 1
#load("~/Dropbox/Glioma_project2/TCGA_AA_DA/anaplastic/AA.cov.corrected.expr.Rdata")
# [1,] "AA.expr.correct"

# Expression matrix of condition 2
#load("~/Dropbox/Glioma_project2/TCGA_AA_DA/diffuse/DA.cov.corrected.expr.Rdata")
# [1,] "DA.expr.correct"

# ------------------------------------------------------------------------------------------------------------
#  Preparing data:
# ------------------------------------------------------------------------------------------------------------
# Make sure there are no values with 0 sd

sd.check<-function(dat_mat,check_rows=F,check_cols=T,verbose=T,help=F){
if(help==T){
  cat("\nUSE: check that values in rows/columns of matrix/dataframe vary: sd>0 or factors are informative: have >1 level / not ids n.levels=n.rows\n")
}

if(is.matrix(dat_mat)){dat_class=apply(dat_mat,2,class)}
if(is.data.frame(dat_mat)){dat_class=unlist(lapply(dat_mat,class))}
  cat("\n+++++++++++++++++ data qc check +++++++++++++++++\n")
  cat("dat_mat contains",ncol(dat_mat),"variables, of which :\n")
  print(table(dat_class))
    dat_num=dat_mat[,dat_class==("numeric"),drop=F]
    dat_fac=dat_mat[,dat_class=="factor",drop=F]
    dat_otr=dat_mat[,!(dat_class%in%c("factor","numeric")),drop=F]

rowind=1:nrow(dat_mat)
colind=1:ncol(dat_mat)

if(check_cols==T){
  fac_col=apply(dat_fac,2,function(x) sum(table(x)!=0))
  num_col=apply(dat_num,2,sd)

  numvarcol=names(num_col)[num_col==0]
    if(length(numvarcol)){cat("\tcol - values do not vary (sd=0) :\t",paste(numvarcol,collapse=", "),"\n")}
  facvarcol=names(fac_col)[fac_col==1]
    if(length(facvarcol)){cat("\tcol - contains single value :\t\t",paste(facvarcol,collapse=", "),"\n")}
  facvaridc=names(fac_col)[fac_col==nrow(dat_fac)]
    if(length(facvaridc)){cat("\tcol - as many levels as rows :\t\t",paste(facvaridc,collapse=", "),"\n")}
    colind=!(colnames(dat_mat)%in% c(numvarcol,facvarcol,facvaridc))
}

if(check_rows==T){
  fac_row=apply(dat_fac,1,function(x) sum(table(x)!=0))
  num_row=apply(dat_num,1,sd)

  numvarrow=names(num_row)[num_row==0]
    if(length(numvarrow)){cat("\trow - values do not vary (sd=0) :\t",paste(numvarrow,collapse=", "),"\n")}
  facvarrow=names(fac_row)[fac_row==1]
    if(length(facvarrow)){cat("\trow - contains single value :\t\t",paste(facvarrow,collapse=", "),"\n")}
  facvaridr=names(fac_row)[fac_row==nrow(dat_fac)]
    if(length(facvaridr)){cat("\trow - as many levels as rows :\t\t",paste(facvaridr,collapse=", "),"\n")}
  rowind=!(rownames(dat_mat)%in% c(numvarrow,facvarrow,facvaridr))
}

if(ncol(dat_otr)>0){
    cat("\tcol - not factor nor numeric :\t",paste(names(dat_otr),collapse=", "),"\n")
}
  cat("\n")

  return(invisible(dat_mat[rowind,colind]))
}



AA.expr.correct=sd.check(AA.expr.correct, check_rows = T)
DA.expr.correct=sd.check(DA.expr.correct, check_rows = T)

# make sure that the genes are in the same order:
rows=intersect(rownames(AA.expr.correct), rownames(DA.expr.correct))
AA.expr.correct=AA.expr.correct[rownames(AA.expr.correct)%in%rows,]
DA.expr.correct=DA.expr.correct[rownames(DA.expr.correct)%in%rows,]
table(rownames(AA.expr.correct)==rownames(DA.expr.correct)) # True

data=list()
data[['condition']]=AA.expr.correct
data[['control']]=DA.expr.correct


# ------------------------------------------------------------------------------------------------------------
# DiffCoEx Function:
# ------------------------------------------------------------------------------------------------------------

wgcna.diffcoex_L <- function (list_expr, pow = 6, minModuleSize = 40, mergeHeight = 0.15,
                            datDescr = "", signType = "unsigned")
{
  print("  NOTE : Input data (list_expr) is expected as a list with each entry : rows = genes, columns = samples")
  library("WGCNA")
  enableWGCNAThreads()
  print("Can deal with only one softpower (pow) at a time ")
  set.seed(0)

  bicorL <- list() # correlation of gene expression
  for (ireg in 1:length(list_expr)) {
    print("■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■")
    print(paste(names(list_expr)[ireg], ireg, "of", length(names(list_expr))))
    t0 = Sys.time()

    COND <- list_expr[[ireg]]
    # CONDav <- scale(COND, scale = F)
    CONDav=COND
    bicorL[[ireg]] <- bicor(t(as.matrix(CONDav))) # bicor correlation of conditions

  }

  names(bicorL) <- names(list_expr)
  collectGarbage()
  t0 = Sys.time()

  softPower = pow
  print(softPower)

  dissTOM <- diffcoex_paper(softPower, bicorL, signtype = signType) # *
  print(Sys.time() - t0)
  collectGarbage()

  geneTree = hclust(as.dist(dissTOM), method = "average")
  collectGarbage()

  dynamicMods = cutreeDynamic(dendro = geneTree, distM = dissTOM,
                              method = "hybrid", cutHeight = 0.996, deepSplit = T,
                              pamRespectsDendro = FALSE, minClusterSize = minModuleSize)

  dynamicColors = labels2colors(dynamicMods)
  collectGarbage()

  Datall <- t(as.data.frame(list_expr)) # merge the conditions
  collectGarbage()

  mergedColor <- mergeCloseModules(Datall, dynamicColors, cutHeight = mergeHeight)$color

  print(paste("mergedColor =", length(unique(mergedColor)),
              unique(mergedColor)))
  print(Sys.time() - t0)
  collectGarbage()
  mstat = as.data.frame(table(mergedColor)) # using dynamicColors not Merged
  mstat = mstat[order(mstat[, 2], decreasing = T), ]
  msta0 = mstat[(mstat[, 1] == "grey"), ]
  msta0$module = "M0"
  mstat = mstat[!(mstat[, 1] == "grey"), ]
  mstat$module = paste0("M", 1:nrow(mstat))
  mstat = rbind(mstat, msta0)
  colnames(mstat) = c("color", "ngenes", "module")
  mstat$color = as.character(mstat$color)

  if (datDescr != "") {
    mstat$module = paste(mstat$module, dat_descr, sep = "_")
  }

  module_list = list()
  for (imod in 1:nrow(mstat)) {
    module_list[[mstat$module[imod]]] = rownames(list_expr[[1]])[mergedColor ==
                                                                   mstat$color[imod]]
  }

  module_list[["bkgrnd"]] = rownames(list_expr[[1]])

  module_expr = list()
  for (ilis in 1:length(list_expr)) {
    for (imod in 1:length(module_list)) {
      module_expr[[names(list_expr)[ilis]]][[names(module_list)[imod]]] = list_expr[[names(list_expr)[ilis]]][module_list[[names(module_list)[imod]]],
                                                                                                              ]
    }
  }
  mbg = as.data.frame("bkgrnd")
  mbg$length = nrow(list_expr[[1]])
  mbg$name = "bkgrnd"
  colnames(mbg) = colnames(mstat)
  mstat = rbind(mstat, mbg)
  readme = "\n\tModules are named based on size M1 - biggest, M0 - unclustered, bkgrnd - all input genes, output contains :\n    \t1. module_list - list containing names of genes in each module\n    \t2. module_expr - expression matrix of all genes in module / input dataset\n    \t3. mstat       - key used to name modules, includes module size\n    \t4. GeneTree     - object to plot the WGCNA style dendrogram\n    \n"
  cat(readme)
  return(invisible(list(module_list = module_list, module_expr = module_expr,
                        mstat = mstat, plotobj = geneTree, readme = readme)))
}

applydiffcoex <- function(beta2,corr_mat_list,signtype=signType) # for multiple conditions
{
  correl=vector(mode="list", length=length(corr_mat_list)); # empty (for adjacencies)
  compDij=vector(mode="list", length=length(corr_mat_list)); # empty

  #compute the cij0 (for all the conditions)
  for (k in 1:length(corr_mat_list))
  {
    datmat=as.matrix(corr_mat_list[[k]])
    correl[[k]]= sign(datmat)*(datmat)^2
    diag(correl[[k]])=0 # change the diagonal to 0 instead of 1
  }
  # correl holds information about adjacency: sign(corr)*(corr)^2


  cij0 = Reduce("+",correl)/length(correl); # Reduce adds two matrices (correl is a list of matrices)

  #compute Dij
  for (element in 1:length(correl))
  {
    compDij[[element]]=abs(correl[[element]]- cij0)/2;
  }

     Dij=(1/(length(corr_mat_list)-1))*Reduce("+", compDij)^(beta2/2) # compDij is a list of matrices

  dissTOM=TOMdist(Dij, TOMType = signtype); # WGCNA function
  collectGarbage()
  return(dissTOM)
}

diffcoex_paper <- function(beta2,bicorL,signtype=signType){
  AdjMatC1=sign(bicorL[[1]])*(bicorL[[1]])^2
  AdjMatC2=sign(bicorL[[2]])*(bicorL[[2]])^2
  message("Condition - Control")

  diag(AdjMatC1)<-0
  diag(AdjMatC2)<-0
  beta1=beta2

  dissTOMC1C2=TOMdist((abs(AdjMatC1-AdjMatC2)/2)^(beta1/2))
} # for 2 conditions (condition-control)



# ------------------------------------------------------------------------------------------------------------
# Run DiffCoEx Function
# ------------------------------------------------------------------------------------------------------------
# Anaplastic vs diffuse IDH1 mut astrocytoma

pw6=wgcna.diffcoex_L(list_expr = data, pow = 6, signType = "unsigned") ; pw6$mstat # 2 modules

cat(pw6$readme)
pw6$mstat

head(pw6$module_list)
head(pw6$module_list[[2]])


#===================================================================================================================================
# Significance of DiffCoEx: the functions
#===================================================================================================================================

permutationProcedureModule2Module<-function(permutation,d,c1,c2,colorh1C1C2)
{
  d1<-d[permutation,]
  d2<-d[-permutation,]
  dispersionModule2Module(c1,c2,d1,d2,colorh1C1C2)
}

dispersionModule2Module<-function(c1,c2,datC1,datC2,colorh1C1C2)
{
  if (c1==c2)
  {
    difCor<-(cor(datC1[,which(colorh1C1C2 == c1)],method="spearman")-
               cor(datC2[,which(colorh1C1C2 == c1)],method="spearman"))
    difCor=difCor^2

    n<-length(which(colorh1C1C2  ==c1))

   #(1/((n^2 -n)/2)*(sum(difCor)/2))^(0.5)
    (sum(difCor)/(n^2 -n))^(0.5)

    # NB! not defining a variable will result in returing the value

  }
  else if (c1!=c2)
  {
    difCor<-(cor(datC1[,which(colorh1C1C2 == c1)],datC1[,which(colorh1C1C2==c2)],method="spearman")-
               cor(datC2[,which(colorh1C1C2 == c1)],datC2[,which(colorh1C1C2==c2)],method="spearman"))^2

    n1<-length(which(colorh1C1C2  ==c1))
    n2<-length(which(colorh1C1C2  ==c2))
    #   (1/((n1*n2))*(sum(difCor)))^(0.5)
    (sum(difCor)/(n1*n2))^(0.5)

  }
}


mod.name.col <- function(module_ensg_list){
  mcol=module_ensg_list$bkgrnd
  for (i in 1:(length(module_ensg_list)-1)) {
    mcol[mcol%in%module_ensg_list[[i]]]=names(module_ensg_list[i])

  }
  return(mcol)
}


# The actual function:
DifCoEx.stats <- function(data.list=data,module_list=module_list, plot=T, OutPath=OutPath, datdesc=""){
  message("***Input data with genes as rows and as samples columns***")
  message("NB! Comparing 2 conditions only")

  datC1=t(data.list[[1]])
  datC2=t(data.list[[2]])
  colorh1C1C2 <- mod.name.col(module_list) # module labels

  permutations<-NULL
  for (i in 1:1000) # or repeat()
  {
    permutations<-rbind(permutations,sample(1:(nrow(datC1)+nrow(datC2)),nrow(datC1)))
  }

  #d<-rbind(scale(datC1),scale(datC2))
  d<-rbind((datC1),(datC2))

  dispersionMatrix<-matrix(nrow=length(unique(colorh1C1C2))-1,ncol=length(unique(colorh1C1C2))-1)

  nullDistrib<-list()
  i<-j<-0
  for (c1 in setdiff(unique(colorh1C1C2),"M0")) # setdiff(x,y) removes y from x (the gray from the modules)
  {
    i<-i+1
    j<-0

    nullDistrib[[c1]]<-list()
    for (c2 in setdiff(unique(colorh1C1C2),"M0"))
    {
      j<-j+1
      dispersionMatrix[i,j]<-dispersionModule2Module(c1,c2,datC1,datC2,colorh1C1C2)

      print("■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■")
      print(paste("Comparing",c1,"and",c2))

      nullDistrib[[c1]][[c2]]<- apply(
        X=permutations
        ,MARGIN = 1
        ,FUN = permutationProcedureModule2Module,d,c2,c1,colorh1C1C2)
    }
  }


  pvalMatrix=matrix(nrow=length(setdiff(unique(colorh1C1C2),"M0")),ncol=length(setdiff(unique(colorh1C1C2),"M0")))
  colnames(pvalMatrix)<-setdiff(unique(colorh1C1C2),"M0")
  rownames(pvalMatrix)<-setdiff(unique(colorh1C1C2),"M0")

  if(plot){

    pdf(paste0(OutPath,'/',"histograms_dispersion_",datdesc,".pdf"), height = 3.5)
  }

  for (i in 1:length(colnames(pvalMatrix))) {
    for (j in 1:length(colnames(pvalMatrix))) {
      pvalMatrix[i,j]<-length(which(nullDistrib[[i]][[j]] >= dispersionMatrix[i,j]))
      # how many random dispresions are larger than observed
      pval=(sum(abs(nullDistrib[[i]][[j]]) > abs(dispersionMatrix[i,j])) + 1) / (length(nullDistrib[[i]][[j]]) + 1)
      pvalMatrix[i,j]=pval

      diff=nullDistrib[[i]][[j]]
      obs=dispersionMatrix[i,j]

      if(plot==T){

        hist(diff, col="gray", xlab="Dispersion",xlim=c(min(diff),obs+0.01), ylab="Frequency"
             , main=paste('Comparing',names(nullDistrib[i]),'and',names(nullDistrib[[i]][j]))
             , breaks=100); abline(v=obs, col="red")

      }
    }
  }

  if(plot){dev.off()}

  return(pvalMatrix)
}


# ------------------------------------------------------------------------------------------------------------
# Significance of DifCoEx Run:
# ------------------------------------------------------------------------------------------------------------
OutPath="~/Dropbox/Glioma_project2/AnaplasticAstro/DIFCOEX_UNSIGNED/"

stats_pw5=DifCoEx.stats(data.list=data,module_list=pw5$module_list, plot=T, OutPath=OutPath, datdesc = "pw5")
Head(stats_pw5)

stats_pw6=DifCoEx.stats(data.list=data,module_list=stuffs$module_list, plot=F, OutPath=OutPath, datdesc = "pw5")
Head(stats_pw6)
