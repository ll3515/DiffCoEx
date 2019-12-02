# DiffCoEx
options(stringsAsFactors=F)

#load("~/Dropbox/Glioma_project2/TCGA_AA_DA/anaplastic/AA.cov.corrected.expr.Rdata")
# [1,] "AA.expr.correct"

#Load("~/Dropbox/Glioma_project2/TCGA_AA_DA/diffuse/DA.cov.corrected.expr.Rdata")
# [1,] "DA.expr.correct"
#Head(DA.expr.correct)

#Load("~/Dropbox/Glioma_project2/TCGA_AA_DA/wt/wt.corrected.expr.Rdata") # only G3
#[,1]
#[1,] "wt.correct"
#Head(wt.correct)

#Load("/Users/ll3515/Dropbox/Glioma_project2/TCGA_AA_DA/ODG/corrected_AO_and_DO_expression.Rdata")
#[,1]
#[1,] "AO.expr.correct"
#[2,] "DO.expr.correct"

Load("~/Dropbox/Glioma_project2/TCGA_AA_DA/IDH1_ASTROCYTOMAS/CovCor_IDH1Astro_expr_clin.Rdata")
#[1,] "idh1astro_clin"
#[2,] "idh1astro_expr_covcor"

# RM 1p/19q_codel
codel = c("TCGA-HT-7607","TCGA-S9-A6U5","TCGA-S9-A6WL","TCGA-CS-5394","TCGA-S9-A6WN","TCGA-P5-A781")
idh1astro_clin = idh1astro_clin[!rownames(idh1astro_clin)%in%codel,]; Head(idh1astro_clin)

  g3=rownames(idh1astro_clin)[idh1astro_clin$tumor_grade=="G3"] # 70
  g2=rownames(idh1astro_clin)[idh1astro_clin$tumor_grade=="G2"] # 49

AA.expr.correct=idh1astro_expr_covcor[,colnames(idh1astro_expr_covcor)%in%g3]
DA.expr.correct=idh1astro_expr_covcor[,colnames(idh1astro_expr_covcor)%in%g2]

# ------------------------------------------------------------------------------------------------------------
#  Preparing data:
# ------------------------------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------------------------------

# Make sure there are no values with 0 sd
AA.expr.correct=sd.check(AA.expr.correct, check_rows = T)
DA.expr.correct=sd.check(DA.expr.correct, check_rows = T)

rows=intersect(rownames(AA.expr.correct), rownames(DA.expr.correct))
AA.expr.correct=AA.expr.correct[rownames(AA.expr.correct)%in%rows,]
DA.expr.correct=DA.expr.correct[rownames(DA.expr.correct)%in%rows,]
table(rownames(AA.expr.correct)==rownames(DA.expr.correct)) # True

data=list()
data[['condition']]=AA.expr.correct
data[['control']]=DA.expr.correct




# ------------------------------------------------------------------------------------------------------------
# Run DifCoEx:
# ------------------------------------------------------------------------------------------------------------

# Anaplastic vs diffuse IDH1 mut astrocytoma

#pw10=wgcna.diffcoex_L(list_expr = data, pow = 10, signType = "unsigned") ; pw10$mstat # 0 modules
#pw9=wgcna.diffcoex_L(list_expr = data, pow = 9, signType = "unsigned") ; pw9$mstat # 0 modules
#pw8=wgcna.diffcoex_L(list_expr = data, pow = 8, signType = "unsigned") ; pw8$mstat # 2 modules
#pw7=wgcna.diffcoex_L(list_expr = data, pow = 7, signType = "unsigned") ; pw7$mstat # 2 modules
#pw6=wgcna.diffcoex_L(list_expr = data, pow = 6, signType = "unsigned") ; pw6$mstat # 2 modules
#pw5=wgcna.diffcoex_L(list_expr = data, pow = 5, signType = "unsigned") ; pw5$mstat #
#pw5=wgcna.diffcoex_L(list_expr = data, pow = 5.5, signType = "unsigned") ; pw5$mstat #
#pw5=wgcna.diffcoex_L(list_expr = data, pow = 5.9, signType = "unsigned") ; pw5$mstat # 1 module
pw5=wgcna.diffcoex_L(list_expr = data, pow = 5.8, signType = "unsigned") ; pw5$mstat #


altpw5 = pw5
#save(altpw5, file = "~/Dropbox/Glioma_project2/TCGA_AA_DA/IDH1_ASTROCYTOMAS/DIFFCOEX/pw5.5_covcor_data_rm_1p19q_rerun.Rdata")
# too big modules
#save(altpw5, file = "~/Dropbox/Glioma_project2/TCGA_AA_DA/IDH1_ASTROCYTOMAS/DIFFCOEX/pw5.8_covcor_data_rm_1p19q_rerun.Rdata")
Load("~/Dropbox/Glioma_project2/TCGA_AA_DA/IDH1_ASTROCYTOMAS/DIFFCOEX/pw5.8_covcor_data_rm_1p19q_rerun.Rdata")

# Head(pw6)
# cat(pw6$readme)
# pw6$mstat

# Head(pw6$module_list)
# Head(pw6$module_list[[2]])

# *** Batch corrected data (for both grades at the same model)
#save(pw6, file="~/Dropbox/Glioma_project2/TCGA_AA_DA/IDH1_ASTROCYTOMAS/DIFFCOEX/pw6_covcor_data.Rdata")
Load("~/Dropbox/Glioma_project2/TCGA_AA_DA/IDH1_ASTROCYTOMAS/DIFFCOEX/pw6_covcor_data.Rdata")


table(altpw5$module_list[[2]]%in%pw6$module_list[[2]])
# modules are lager, but core elements do not change

pw5.mods = altpw5$module_list

M1 = pw5.mods[[1]]
M2 = pw5.mods[[2]]
bkgrnd = rownames(AA.expr.correct)
write.delim(M1, file = "M1.txt", col.names = F, row.names = F)
write.delim(M2, file = "M2.txt", col.names = F, row.names = F)
write.delim(bkgrnd, file = "bkgrnd.txt", col.names = F, row.names = F)
