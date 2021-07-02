library("TCGAbiolinks")
library(GEOquery)
library(tidyr)
library(preprocessCore)
library(ROCR)
library(DESeq2)
library(biomaRt)

pacman::p_load("magrittr", "stringr", "readr", "stringr", "openxlsx", "readxl", "FactoMineR", "factoextra", "ggplot2",
               "ggpubr", "grid", "gridExtra", "ggdendro", "dendroextras", "ggthemes", "plyr", "dplyr", "tidyr", "venn", "pca3d", "survival",
               "survminer", "XML", "methods", "rvest", "scales")

source("C:/Users/raevs/Desktop/Oncobox/LTS_vs_STS_manuscript/scripts/preprocessing.R", chdir = TRUE)
source("C:/Users/raevs/Desktop/Oncobox/LTS_vs_STS_manuscript/scripts/visualization.R", chdir = TRUE)

dataFolder <- "C:/Users/raevs/Desktop/Oncobox/LTS_vs_STS_manuscript/gliomas_pathways_PFS/"
setwd(dataFolder)

hgnc <- read_tsv("C:/Users/raevs/Desktop/Oncobox/LTS_vs_STS_manuscript/gliomas_pathways_PFS/hgnc_complete_set.txt", col_names = TRUE)
our_pathways <- read.xlsx("C:/Users/raevs/Desktop/Oncobox/LTS_vs_STS_manuscript/gliomas_pathways_PFS/number_of_genes_in_pathways.xlsx")
dat_oncobox <- read_delim("C:/Users/raevs/Desktop/Oncobox/LTS_vs_STS_manuscript/exp_data/slovenia_exp_matrix.txt",
        "\t", escape_double = FALSE, trim_ws = TRUE)
our_genes <- dat_oncobox$SYMBOL

#####################
### TCGA analysis ###
#####################

# query <- GDCquery(project = "TCGA-LGG",
#                   data.category = "Clinical",
#                   data.type = "Clinical Supplement",
#                   data.format = "BCR Biotab")
# GDCdownload(query)
# clinical.BCRtab.all <- GDCprepare(query)
# colnames(clinical.BCRtab.all$clinical_radiation_lgg)
# names(clinical.BCRtab.all)
# names(clinical.BCRtab.all$clinical_patient_brca)


library("XML")
library("methods")
library("rvest")

files <- list.files('clinical_xml_lgg_gbm_TCGA')
files <- files[str_detect(files, 'xml')]
idh_statuses <- data.frame(matrix(ncol = 3, nrow = 0))
x <- c("bcr_patients_barcode", "idh_statuse", "type")
colnames(idh_statuses) <- x

for (file in files) {
  result <- xmlParse(file = paste(c("./clinical_xml_lgg_gbm_TCGA/", file), collapse = ''))
  name = str_sub(file, start = 34, end = -1L - 4)
  print(name)
  idh_statuses[name, 'bcr_patients_barcode'] = name
  rootnode <- xmlRoot(result)
  if (length(intersect(names(rootnode[[2]]), c("ldh1_mutation_found"))) == 1) {
    id_raw = rootnode[[2]][["ldh1_mutation_found"]]
    id_raw = xmlToList(id_raw, addAttributes = TRUE, simplify = FALSE)
    print(id_raw)
    if (length(intersect(names(id_raw), c("text"))) == 1) {
      idh_statuses[name, 'idh_statuse'] = id_raw$text
      print(id_raw$text)
    }
  }
  if (length(intersect(names(rootnode[[2]]), c("histological_type"))) == 1) {
    id_raw = rootnode[[2]][["histological_type"]]
    id_raw = xmlToList(id_raw, addAttributes = TRUE, simplify = FALSE)
    print(id_raw)
    if (length(intersect(names(id_raw), c("text"))) == 1) {
      idh_statuses[name, 'type'] = id_raw$text
      print(id_raw$text)
    }
  }


}
idh_statuses = na.omit(idh_statuses)
write.xlsx(idh_statuses, 'idh_statuses_LGG_TCGA.xlsx')

data = read_tsv('gdc_download_20201117_132227.209144/3bba81d3-f253-42ea-b9eb-72bc8641fd69/nationwidechildrens.org_clinical_follow_up_v1.0_lgg.txt')
data = read_tsv('C:/Users/raevs/Desktop/Oncobox/LTS_vs_STS_manuscript/gliomas_pathways_PFS/gdc_download_20201117_132243.033093/c9cdbc76-105d-429b-9fce-f000819716f9/nationwidechildrens.org_clinical_follow_up_v1.0_gbm.txt')
data = read_tsv('C:/Users/raevs/Desktop/Oncobox/LTS_vs_STS_manuscript/gliomas_pathways_PFS/gdc_download_20201117_132251.493810/86a66df6-0e51-4555-81ed-fbd2b245de94/nationwidechildrens.org_clinical_follow_up_v1.0_nte_gbm.txt')


data_lgg = data[3:nrow(data),]
data_lgg = data_lgg[, c('bcr_patient_barcode', 'bcr_followup_barcode', 'death_days_to')]
data_lgg = data_lgg[, c('bcr_patient_barcode', 'bcr_followup_barcode', 'new_tumor_event_dx_days_to')]
#data_alt=read_excel('../primary_vs_recurrent_new/miRNA_exprs_IDH_survival.xlsx')[,c("Patient","Days")]
#colnames(data_alt)[1]='bcr_patient_barcode'
#data_lgg=left_join(data_lgg,data_alt,by='bcr_patient_barcode')
#data_disc=data_all[!(data_all$death_days_to==data_all$Days),]
#data_lgg$Days=as.character(data_lgg$Days)
#data_lgg$death_days_to=as.character(data_lgg$death_days_to)
#data_lgg[is.na(data_lgg$Days),'Days']=''
#data_lgg[(data_lgg$death_days_to=="[Not Applicable]")|(data_lgg$death_days_to=="[Discrepancy]"),'death_days_to']=data_lgg[(data_lgg$death_days_to=="[Not Applicable]")|(data_lgg$death_days_to=="[Discrepancy]"),'Days']
### OVERALL SURVIVAL
data_lgg = data_lgg[(data_lgg$death_days_to != "[Not Applicable]") & (data_lgg$death_days_to != "[Discrepancy]") & (data_lgg$death_days_to != "[Not Available]"), c('bcr_patient_barcode', 'death_days_to')]
data_lgg = data_lgg[(data_lgg$new_tumor_event_dx_days_to != "[Not Applicable]") & (data_lgg$new_tumor_event_dx_days_to != "[Not Available]"), c('bcr_patient_barcode', 'new_tumor_event_dx_days_to')]
#idh_data=read_excel('C:/Users/raevs/Desktop/Oncobox/LTS_vs_STS_manuscript/primary_vs_recurrent_new/miRNA_exprs_IDH_survival.xlsx')
#idh_data=read_excel('C:/Users/raevs/Desktop/Oncobox/LTS_vs_STS_manuscript/primary_vs_recurrent_new/mutations.xlsx')
#idh_data_new=read_excel('C:/Users/raevs/Desktop/Oncobox/LTS_vs_STS_manuscript/gliomas_pathways_PFS/mutations_GBM_mgmt.xlsx')
#idh_data=read_excel('mutations_LGG.xlsx')
#patients_needed_idh_mut=idh_data[!is.na(idh_data$overall_status),'barcode']$barcode
#patients_needed_idh_wt=idh_data[is.na(idh_data$overall_status),'barcode']$barcode
#mgmt_data=read_xls('TCGA_GBM_MGMT_status.xls',col_names = c('bcr_patient_barcode','dataset','Subtype','CIMP','MGMT-STP27_response','MGMT-STP27_class'))
#mgmt_data=mgmt_data[-c(1),]
#patients_needed_met=mgmt_data[mgmt_data$`MGMT-STP27_class`=='M','bcr_patient_barcode']$bcr_patient_barcode
#patients_needed_un=mgmt_data[mgmt_data$`MGMT-STP27_class`=='U','bcr_patient_barcode']$bcr_patient_barcode
#patients_needed=intersect(patients_needed_met,patients_needed_idh_wt)
#patients_needed=na.omit(idh_data)$Patient
#patients_needed=idh_data[is.na(idh_data$IDH_mutation_status),'Patient']$Patient
#patients_needed=idh_statuses[idh_statuses$idh_statuse=='NO','bcr_patients_barcode']
#length(intersect(data_lgg$bcr_patient_barcode,patients_needed))
#length(intersect(idh_statuses$bcr_patients_barcode,str_sub(colnames(dat_TCGA),end=12)))
#data_lgg=data_lgg[!is.na(match(data_lgg$bcr_patient_barcode,patients_needed)),]

data_lgg = data_lgg[!duplicated(data_lgg),]
colnames(data_lgg)[2] = 'death_days_to'
#data_lgg$new_tumor_event_dx_days_to=as.numeric(as.character(data_lgg$new_tumor_event_dx_days_to))
data_lgg$death_days_to = as.numeric(as.character(data_lgg$death_days_to))
data_lgg = na.omit(data_lgg)
data_lgg = aggregate(. ~ bcr_patient_barcode, data = data_lgg, FUN = min)

#length(unique(data_lgg$bcr_patient_barcode))

dat_TCGA = read.csv('C:/Users/raevs/Desktop/Oncobox/LTS_vs_STS_manuscript/MSI/TCGA_read_counts.csv', row.names = 1)
colnames(dat_TCGA) = str_replace_all(colnames(dat_TCGA), pattern = '[.]', replacement = '-')
exprs_lgg = dat_TCGA[, !is.na(match(str_sub(colnames(dat_TCGA), end = 1L + 11), data_lgg$bcr_patient_barcode))]
exprs_annotation = data.frame(bcr_patient_barcode = str_sub(colnames(exprs_lgg), end = 1L + 11),
                            sample_id = colnames(exprs_lgg))
data_lgg = merge(data_lgg, exprs_annotation, by = 'bcr_patient_barcode')
counts_lgg = as.data.frame(dat_TCGA[, as.character(data_lgg$sample_id)])
counts_lgg = t(t(as.matrix(counts_lgg)) / estimateSizeFactorsForMatrix(as.matrix(counts_lgg)))
write.csv(counts_lgg, 'counts_gbm_TCGA.csv', row.names = TRUE)
#counts_lgg=read_csv('counts_lgg.csv')
colnames(counts_lgg) = paste0('tumour_', colnames(counts_lgg))
counts_lgg = counts_lgg + 1
counts_lgg = as.data.frame(counts_lgg)
counts_lgg$norm = c(10 ** (rowMeans(log10(counts_lgg))))

a = counts_lgg[, 120:127]
counts_lgg$SYMBOL = row.names(counts_lgg)
counts_lgg = counts_lgg[, c("SYMBOL", colnames(counts_lgg)[1:(ncol(counts_lgg) - 1)])]
write.csv(counts_lgg, 'counts_gbm_PFS_pathways.csv', row.names = FALSE)

pathways = as.data.frame(cbind(row.names(counts_lgg), 10 ** rowMeans(log10(counts_lgg[, data_lgg[data_lgg$death_days_to < median(data_lgg$death_days_to), "sample_id"]])),
                             10 ** rowMeans(log10(counts_lgg[, data_lgg[data_lgg$death_days_to > median(data_lgg$death_days_to), "sample_id"]]))))
colnames(pathways) = c("SYMBOL", "Tumour_Bad", "Norm_Good")
write.csv(pathways, 'counts_gbm_PFS_pathways_good_bad.csv', row.names = FALSE)



#survival analysis
pathways = counts_lgg
p = pathways$SYMBOL
pathways$SYMBOL = NULL
pathways$norm = NULL

pathways = read_excel('Mainz_primary_GBM_PFS.xlsx', sheet = 2)
pathways$Database = NULL

our_pathways = read.xlsx("C:/Users/raevs/Desktop/Oncobox/LTS_vs_STS_manuscript/gliomas_pathways_PFS/number_of_genes_in_pathways.xlsx")
pathways = pathways[!is.na(match(pathways$Pathway, our_pathways$pathway)),]
p = pathways$Pathway
pathways$Pathway = NULL
row.names(pathways) = p
colnames(pathways) = str_sub(colnames(pathways), 8, -1L)
pathways = as.data.frame(t(pathways))

#pathways=pathways-1
pathways$bcr_patient_barcode = c(str_sub(row.names(pathways), 1, 3))
pathways$bcr_patient_barcode = c(str_sub(row.names(pathways), 1, 12))
pathways = aggregate(. ~ bcr_patient_barcode, data = pathways, FUN = mean)
p = pathways$bcr_patient_barcode
pathways$bcr_patient_barcode = NULL
row.names(pathways) = p

intersect(row.names(pathways), as.character(data_lgg$sample_id))

library("survival")
library("survminer")

HRs = c()
log_rank = c()
HR_p = c()
i = 1
for (pathway in colnames(pathways)) {
  #pathway=colnames(pathways)[1]
  print(i)
  i = i + 1
  #pathway='CRNDE'
  dat2 = data.frame(bcr_patient_barcode = data_lgg$`Patient ID`,
                  time = data_lgg$death_days_to,
                  BES = pathways[as.character(data_lgg$`Patient ID`), pathway], status = 1,
                  Censored = 1)

  #dat2=data.frame(bcr_patient_barcode=row.names(pathways),
  #                time=TTP[row.names(pathways),][,2],
  #                BES=pathways[,pathway],status=1,
  #                Censored=1)
  dat2 = dat2[!duplicated(dat2),]
  colnames(dat2)[2] = 'time'
  dat2$time = as.numeric(as.character(dat2$time))
  dat2 = aggregate(. ~ bcr_patient_barcode, data = dat2, FUN = min)
  median = median(dat2$BES)
  dat2$time = dat2$time / 30
  dat2 = dat2[(dat2$BES >= quantile(dat2$BES, 0.6667)) | (dat2$BES <= quantile(dat2$BES, 0.33333)),]
  res.cox <- coxph(Surv(time = dat2$time, event = dat2$Censored) ~ dat2$BES > median, data = dat2)
  HRs[pathway] = res.cox$coefficients
  a = summary(res.cox)
  log_rank[pathway] = a$logtest[3]
  HR_p[pathway] = a$coefficients[5]
}



#p_adj=p.adjust(p_values,method = "fdr")
sorted_p = sort(log_rank)
#final=data.frame(gene=names(sorted_p),padj=sorted_p,pvalue=p_values[names(sorted_p)],HR=HRs[names(sorted_p)])
#final=data.frame(pathways=names(sorted_p),padj=sorted_p,pvalue=p_values[names(sorted_p)],HR=HRs[names(sorted_p)])


final = data.frame(gene = as.character(names(sorted_p)), log_rank = sorted_p, HR_p = HR_p[names(sorted_p)], HR = HRs[names(sorted_p)])
final$gene = as.character(final$gene)
final = data.frame(pathways = as.character(names(sorted_p)), log_rank = sorted_p, HR_p = HR_p[names(sorted_p)], HR = HRs[names(sorted_p)])
final$pathways = as.character(final$pathways)

library(venn)
venn(list(HR = final[final$HR_p < 0.05, 'gene'], log_rank = final[final$log_rank < 0.05, 'gene']))

#final=final[order(final$padj),]
#final=na.omit(final[final$padj<0.1,])
#final_old=read_excel('CGGA_recurrent_LGG_OS_HRs_pathways.xlsx')
#final$pathway=row.names(final)
write.xlsx(final, 'Mainz_PFS_HRs_pathways_third.xlsx')

final_old = read_excel('GBM_OS_HRs_pathways.xlsx')

fit <- survfit(Surv(time = dat2$time, event = dat2$Censored) ~ dat2$BES > median, data = dat2)

pdf(file = "CRNDE_Slovenia_KM.pdf")
ggsurvplot(fit,
           pval = TRUE,
           pval.method = TRUE,
           size = 2,
           legend.labs = c("lower than 1/3", "higher than 1/3"),
           pval.size = 7,
           font.x = 25,
           font.y = 25,
           font.tickslab = 23,
           xlab = "Time (months)",
           ylab = "Overall survival probability",
           risk.table = TRUE,
           pval.coord = c(30, 0.9),
           pval.method.coord = c(30, 1.0))
dev.off()
dat2$Status = as.numeric(dat2$BES > median)
model <- coxph(Surv(time = dat2$time, event = dat2$Censored) ~ dat2$Status, data = dat2[, c('time', 'Censored', 'Status')])
pdf(file = "CRNDE_Slovenia_forest.pdf", width = 7, height = 1.5)
ggforest(model)
dev.off()


#######################
### IVYGAP analysis ###
#######################

annotation = read.csv('IVYGAP_annotation.csv')
annotation = annotation[!is.na(annotation$survival_days),]
annotation$tumor_name = as.character(annotation$tumor_name)
bams = read.csv('IVYGAP_bams.csv')
bams$tumor_name = as.character(bams$tumor_name)
annotation = left_join(annotation[, c("tumor_name", "survival_days")], bams[, c("tumor_name", "block_name", "specimen_name", "rna_integrity_number",
                                                                          "structure_name", "files_contents", "structure_acronym", "bam_download_link", "bai_download_link")])
length(unique(annotation[annotation$structure_name == 'Cellular Tumor sampled by reference histology', 'tumor_name']))
length(unique(annotation[str_detect(annotation$structure_name, 'Cellular Tumor sampled by'), 'tumor_name']))


length(intersect(annotation$tumor_name, bams$tumor_name))
write.xlsx(annotation, 'IVYGAP_annotation_and_files.xlsx')

#####################
### CGGA analysis ###
#####################
CGGA_annotation <- read_tsv('CGGA.mRNAseq_693_clinical.20200506.txt')
# CGGA_annotation_C <- read_tsv('CGGA.mRNAseq_325_clinical.20200506.txt')

# CGGA_annotation <- rbind(CGGA_annotation_B, CGGA_annotation_C)
CGGA_annotation <- CGGA_annotation[!is.na(CGGA_annotation$OS),]
CGGA_annotation <- CGGA_annotation[!is.na(CGGA_annotation$PRS_type),]
CGGA_annotation <- CGGA_annotation[CGGA_annotation$PRS_type != 'Secondary',]

CGGA_annotation_short <- na.omit(CGGA_annotation[CGGA_annotation$Histology == 'GBM', c("CGGA_ID", "IDH_mutation_status", "MGMTp_methylation_status")])
nrow(CGGA_annotation_short[(CGGA_annotation_short$IDH_mutation_status == 'Wildtype') & (CGGA_annotation_short$MGMTp_methylation_status == 'un-methylated'),])

counts <- read_tsv("CGGA.mRNAseq_693.RSEM-genes.20200506.txt")
# counts_C <- read_tsv("CGGA.mRNAseq_325.RSEM-genes.20200506.txt")
counts <- as.data.frame(counts)

# counts <- merge(counts_B, counts_C, by = 'Gene_Name')
p <- counts$Gene_Name
row.names(counts) <- p
counts$Gene_Name <- NULL
counts <- counts[CGGA_annotation$CGGA_ID]
counts <- round(counts)

counts <- t(t(as.matrix(counts)) / estimateSizeFactorsForMatrix(as.matrix(counts)))
# write.csv(counts_lgg, "expression_CGGA_693.csv", row.names = TRUE)
#counts_lgg=read_csv("counts_lgg.csv")

###### primary LGG OS ######

data_lgg <- CGGA_annotation[(CGGA_annotation$PRS_type == 'Primary') & str_detect(CGGA_annotation$Histology, 'GBM'), c("CGGA_ID", "OS", "Censor (alive=0; dead=1)")]

#data_lgg=CGGA_annotation[(CGGA_annotation$Histology!='GBM')&(CGGA_annotation$IDH_mutation_status=='Wildtype'),c("CGGA_ID","OS","Censor (alive=0; dead=1)")]
#data_lgg=CGGA_annotation[(CGGA_annotation$Histology=='GBM')&(CGGA_annotation$MGMTp_methylation_status=='methylated')&(CGGA_annotation$IDH_mutation_status=='Wildtype'),c("CGGA_ID","OS","Censor (alive=0; dead=1)")]

data_lgg <- na.omit(data_lgg)

colnames(data_lgg) <- c('bcr_patient_barcode', 'death_days_to', 'Censored')
data_lgg$death_days_to <- as.numeric(as.character(data_lgg$death_days_to))
data_lgg <- na.omit(data_lgg)

counts_lgg <- as.data.frame(counts[, as.character(data_lgg$bcr_patient_barcode)])
write.csv(counts_lgg, 'counts_CGGA_693_GBM_OS.csv', row.names = TRUE)
a <- read_csv('counts_CGGA_693_GBM_OS.csv')

# colnames(counts_lgg) <- paste0("tumour_", colnames(counts_lgg))
counts_lgg <- counts_lgg + 1
counts_lgg <- as.data.frame(counts_lgg)
counts_lgg$norm <- c(10 ** (rowMeans(log10(counts_lgg))))
counts_lgg$SYMBOL <- row.names(counts_lgg)
counts_lgg <- counts_lgg[, c("SYMBOL", colnames(counts_lgg)[1:(ncol(counts_lgg) - 1)])]
write.csv(counts_lgg, "CGGA_693_counts_primary_GBM_MGMT_met_IDH_WT_OS_pathways.csv", row.names = FALSE)

data_lgg <- data_lgg[data_lgg$Censored == 1,]
pathways <- as.data.frame(cbind(row.names(counts_lgg), 10 ** rowMeans(log10(counts_lgg[, data_lgg[data_lgg$death_days_to < median(data_lgg$death_days_to), "bcr_patient_barcode"]$bcr_patient_barcode])),
                             10 ** rowMeans(log10(counts_lgg[, data_lgg[data_lgg$death_days_to > median(data_lgg$death_days_to), "bcr_patient_barcode"]$bcr_patient_barcode]))))


# x <- data_lgg[data_lgg$death_days_to < median(data_lgg$death_days_to), "bcr_patient_barcode"]
# counts_lgg_x <- counts_lgg[, x$bcr_patient_barcode]

# pathways <- as.data.frame(cbind(
#   row.names(counts_lgg), 10 ** rowMeans(log10(counts_lgg[, data_lgg[data_lgg$death_days_to < median(data_lgg$death_days_to), "bcr_patient_barcode"]$bcr_patient_barcode])),
#                              10 ** rowMeans(log10(counts_lgg[, data_lgg[data_lgg$death_days_to > median(data_lgg$death_days_to), "bcr_patient_barcode"]$bcr_patient_barcode]))))

colnames(pathways) <- c("SYMBOL", "Tumour_Bad", "Norm_Good")
write.csv(pathways, "CGGA_693_counts_primary_GBM_OS_pathways_good_bad.csv", row.names = FALSE)




##### Slovenia PFS #####
annotation = read_excel("KopijaAB_Samples for NGS 14.03.2019-2_ENG.xlsx", sheet = 2)
annotation = na.omit(annotation[, c("Oncobox ID", "Survival (months)")])
#annotation=na.omit(annotation[,c("Oncobox ID",'TTP (months)')])
colnames(annotation)[2] = 'survival months'
annotation = annotation[(annotation$`survival months` != '/') & (annotation$`survival months` != 'FALSE') & (annotation$`survival months` != 'survival months'),]
annotation = separate_rows(annotation, 'Oncobox ID', sep = ', ')
annotation$`Oncobox ID` = str_remove_all(annotation$`Oncobox ID`, "_")
row.names(annotation) = annotation$`Oncobox ID`

names = c('GB-1_S11_R1_001', 'GB_10_S5_R1_001', 'GB_11_S6_R1_001', 'GB_12_S7_R1_001', 'GB_16_S8_R1_001',
'GB_17_S9_R1_001', 'GB_18_S10_R1_001', 'GB_2_S4_R1_001', 'GB_3_S5_R1_001', 'GB_4_S6_R1_001', 'GB_5_S5_R1_001',
'GB_6_S9_R1_001', 'GB_7_S10_R1_001', 'GB_8_S4_R1_001')

# dat_oncobox = read_csv('C:/Users/raevs/Desktop/Oncobox/LTS_vs_STS_manuscript/CNV/oncobox_read_counts_ALL.csv')
dat_oncobox <- read_delim("C:/Users/raevs/Desktop/Oncobox/LTS_vs_STS_manuscript/exp_data/slovenia_exp_matrix.txt",
        "\t", escape_double = FALSE, trim_ws = TRUE)
colnames(dat_oncobox) <- c("SYMBOL", names)
p = dat_oncobox$SYMBOL
# dat_oncobox = dat_oncobox[, names]
dat_oncobox = as.data.frame(dat_oncobox)
dat_oncobox$SYMBOL <- NULL
row.names(dat_oncobox) = p
colnames(dat_oncobox) = str_sub(colnames(dat_oncobox), 1, -1L - 10)
colnames(dat_oncobox) = str_remove_all(colnames(dat_oncobox), "_")
colnames(dat_oncobox) = str_remove_all(colnames(dat_oncobox), "-")
genes = row.names(dat_oncobox)
annotation = annotation[intersect(annotation$`Oncobox ID`, colnames(dat_oncobox)),]
dat_oncobox = dat_oncobox[, annotation$`Oncobox ID`]
dat_oncobox = t(t(as.matrix(dat_oncobox)) / estimateSizeFactorsForMatrix(as.matrix(dat_oncobox)))
data_lgg = annotation
data_lgg$Censored = 1
colnames(data_lgg) = c('bcr_patient_barcode', 'death_days_to', 'Censored')
data_lgg$death_days_to = as.numeric(data_lgg$death_days_to) * 30

counts_lgg = as.data.frame(dat_oncobox)
row.names(counts_lgg) = p
colnames(counts_lgg) = paste0('tumour_', colnames(counts_lgg))
counts_lgg = counts_lgg + 1
counts_lgg = as.data.frame(counts_lgg)
counts_lgg$norm = c(10 ** (rowMeans(log10(counts_lgg))))
counts_lgg$SYMBOL = genes
counts_lgg = counts_lgg[, c("SYMBOL", colnames(counts_lgg)[1:(ncol(counts_lgg) - 1)])]
write.csv(counts_lgg, 'Slovenia_counts_primary_GBM_PFS_pathways.csv', row.names = FALSE)
pathways = as.data.frame(cbind(row.names(counts_lgg), 10 ** rowMeans(log10(counts_lgg[, data_lgg[data_lgg$death_days_to < median(data_lgg$death_days_to), "bcr_patient_barcode"]$bcr_patient_barcode])),
                             10 ** rowMeans(log10(counts_lgg[, data_lgg[data_lgg$death_days_to > median(data_lgg$death_days_to), "bcr_patient_barcode"]$bcr_patient_barcode]))))
colnames(pathways) = c("SYMBOL", "Tumour_Bad", "Norm_Good")
write.csv(pathways, 'Slovenia_counts_primary_GBM_OS_pathways_good_bad.csv', row.names = FALSE)




a = read_csv('Slovenia_counts_primary_GBM_PFS_pathways.csv')

##### Mainz #####
files = read_excel('Annotation_Mainz.xls', 1)
TTP = read_excel('Annotation_Mainz.xls', 2)
row.names(TTP) = TTP$`Patient ID`
files = files[files$`source name` == 'Glioblastoma, primary tumor', c("characteristics: Patient",
                                                                 "title")]
colnames(files) = c("Patient ID", "bcr_patient_barcode")
data_lgg = left_join(files, TTP, by = 'Patient ID')
#data_lgg$bcr_patient_barcode=str_sub(data_lgg$bcr_patient_barcode,end = -1L-9)
data_lgg = na.omit(data_lgg)
length(unique(data_lgg$`Patient ID`))
data_lgg$Censored = 1
data_lgg$death_days_to = as.numeric(data_lgg$TTP) * 30

counts_lgg = read_tsv('GSE139533_mergede_protein_coding_counts.tab')
p = counts_lgg$SYMBOL
counts_lgg = counts_lgg[, data_lgg$bcr_patient_barcode]
row.names(counts_lgg) = p
write.csv(counts_lgg, 'Mainz_counts_.csv', row.names = TRUE)

row.names(counts_lgg) = p
counts_lgg = t(t(as.matrix(counts_lgg)) / estimateSizeFactorsForMatrix(as.matrix(counts_lgg)))

colnames(counts_lgg) = paste0('tumour_', colnames(counts_lgg))
counts_lgg = counts_lgg + 1
counts_lgg = as.data.frame(counts_lgg)
counts_lgg$norm = c(10 ** (rowMeans(log10(counts_lgg))))
counts_lgg$SYMBOL = row.names(counts_lgg)
counts_lgg = counts_lgg[, c("SYMBOL", colnames(counts_lgg)[1:(ncol(counts_lgg) - 1)])]
write.csv(counts_lgg, 'Mainz_counts_primary_GBM_PFS_pathways.csv', row.names = FALSE)

pathways = as.data.frame(cbind(row.names(counts_lgg), 10 ** rowMeans(log10(counts_lgg[, data_lgg[data_lgg$death_days_to < median(data_lgg$death_days_to), "bcr_patient_barcode"]$bcr_patient_barcode])),
                             10 ** rowMeans(log10(counts_lgg[, data_lgg[data_lgg$death_days_to > median(data_lgg$death_days_to), "bcr_patient_barcode"]$bcr_patient_barcode]))))
colnames(pathways) = c("SYMBOL", "Tumour_Bad", "Norm_Good")
write.csv(pathways, 'Mainz_counts_primary_GBM_PFS_pathways_good_bad.csv', row.names = FALSE)



##### comparisons Glioblsatoma PFS
## IDH WT


# final_TCGA_OS = read_excel('GBM_PFS_HRs_pathways_third.xlsx')
# final_TCGA_OS = na.omit(final_TCGA_OS)
final_TCGA_OS = read_excel('GBM_OS_HRs_pathways_third.xlsx')
final_TCGA_OS = na.omit(final_TCGA_OS)
final_TCGA_OS <- as.data.frame(final_TCGA_OS)
final_CGGA_OS = read_excel('CGGA_primary_GBM_OS_HRs_pathways_third.xlsx')
final_CGGA_OS = na.omit(final_CGGA_OS)
final_CGGA_OS <- as.data.frame(final_CGGA_OS)
# final_Slovenia_OS = read_excel('Slovenia_primary_GBM_PFS_HRs_pathways_third.xlsx')
# final_Slovenia_OS = na.omit(final_Slovenia_OS)
final_Slovenia_OS = read_excel('Slovenia_primary_GBM_OS_HRs_pathways_third.xlsx')
final_Slovenia_OS = na.omit(final_Slovenia_OS)
final_Slovenia_OS <- as.data.frame(final_Slovenia_OS)
final_Mainz_OS = read_excel('Mainz_PFS_HRs_pathways_third.xlsx')
final_Mainz_OS = na.omit(final_Mainz_OS)
final_Mainz_OS <- as.data.frame(final_Mainz_OS)

# a = read_excel('counts_gbm_PFS_pathways_good_bad.xlsx')
a = read_excel('counts_gbm_OS_pathways_good_bad.xlsx')
a <- as.data.frame(a)
colnames(a)[1] = 'pathways'
final_TCGA_OS = left_join(final_TCGA_OS, a, by = 'pathways')
final_TCGA_OS = final_TCGA_OS[(final_TCGA_OS$log_rank < 0.05) & (final_TCGA_OS$HR_p < 0.05),]
final_TCGA_OS = na.omit(final_TCGA_OS[order(final_TCGA_OS$Tumour_Bad),])
final_TCGA_OS = rbind(final_TCGA_OS[1:10,], final_TCGA_OS[(nrow(final_TCGA_OS) - 9):nrow(final_TCGA_OS),])

# a = read_excel('CGGA_counts_primary_GBM_PFS_pathways_good_bad.xlsx')
# a = read_excel('CGGA_counts_primary_GBM_OS_pathways_good_bad.xlsx')
a = read_excel('CGGA_693_counts_primary_GBM_OS_pathways_good_bad.xlsx')
a <- as.data.frame(a)
colnames(a)[1] = 'pathways'
final_CGGA_OS = left_join(final_CGGA_OS, a, by = 'pathways')
final_CGGA_OS = final_CGGA_OS[(final_CGGA_OS$log_rank < 0.05) & (final_CGGA_OS$HR_p < 0.05),]
final_CGGA_OS = final_CGGA_OS[order(final_CGGA_OS$Tumour_Bad),]
final_CGGA_OS = rbind(final_CGGA_OS[1:10,], final_CGGA_OS[(nrow(final_CGGA_OS) - 9):nrow(final_CGGA_OS),])

# a = read_excel('Slovenia_counts_primary_GBM_PFS_pathways_good_bad.xlsx')
a = read_excel('Slovenia_counts_primary_GBM_OS_pathways_good_bad.xlsx')
a <- as.data.frame(a)
colnames(a)[1] = 'pathways'
final_Slovenia_OS = left_join(final_Slovenia_OS, a, by = 'pathways')
final_Slovenia_OS = final_Slovenia_OS[(final_Slovenia_OS$log_rank < 0.05) & (final_Slovenia_OS$HR_p < 0.05),]
final_Slovenia_OS = na.omit(final_Slovenia_OS[order(final_Slovenia_OS$Tumour_Bad),])
final_Slovenia_OS = rbind(final_Slovenia_OS[1:10,], final_Slovenia_OS[(nrow(final_Slovenia_OS) - 9):nrow(final_Slovenia_OS),])



genes_TCGA = final_TCGA_OS$gene
genes_CGGA = final_CGGA_OS$gene
genes_Slovenia = final_Slovenia_OS$gene
# genes_Mainz = final_Mainz_OS$gene
pathways = final_Mainz_OS$pathways

library(venn)
# my_list_minus = list(TCGA = final_TCGA_OS[(final_TCGA_OS$pvalue < 0.05) & (final_TCGA_OS$HR < 0), 'gene']$pathways,
#                    CGGA = final_CGGA_OS[(final_CGGA_OS$pvalue < 0.05) & (final_CGGA_OS$HR < 0), 'pathways']$pathways)
my_list_minus = list(TCGA = final_TCGA_OS[(final_TCGA_OS$pvalue < 0.05) & (final_TCGA_OS$HR < 0), 'gene']$pathways,
                   CGGA = final_CGGA_OS[(final_CGGA_OS$pvalue < 0.05) & (final_CGGA_OS$HR < 0), 'gene']$pathways)

my_list_plus = list(TCGA = final_TCGA_OS[(final_TCGA_OS$pvalue < 0.05) & (final_TCGA_OS$HR > 0), 'pathways']$pathways,
                  CGGA = final_CGGA_OS[(final_CGGA_OS$pvalue < 0.05) & (final_CGGA_OS$HR > 0), 'pathways']$pathways)


my_list_plus = list(TCGA_GBM_OS = final_TCGA_OS[(final_TCGA_OS$Tumour_Bad > 0), 'pathways']$pathways,
                  CGGA_GBM_OS = final_CGGA_OS[(final_CGGA_OS$Tumour_Bad > 0), 'pathways']$pathways,
                  Slovenia_GBM_OS = final_Slovenia_OS[(final_Slovenia_OS$Tumour_Bad > 0), 'pathways']$pathways)
my_list_minus = list(TCGA_GBM_OS = final_TCGA_OS[(final_TCGA_OS$Tumour_Bad < 0), 'pathways']$pathways,
                   CGGA_GBM_OS = final_CGGA_OS[(final_CGGA_OS$Tumour_Bad < 0), 'pathways']$pathways,
                   Slovenia_GBM_OS = final_Slovenia_OS[(final_Slovenia_OS$Tumour_Bad < 0), 'pathways']$pathways)

# my_list_plus = list(TCGA_GBM_OS = final_TCGA_OS[(final_TCGA_OS$log_rank < 0.05) & (final_TCGA_OS$HR_p < 0.05) & (final_TCGA_OS$HR > 0), 'gene']$gene,
#                   Mainz_GBM_OS = final_Mainz_OS[(final_Mainz_OS$log_rank < 0.05) & (final_Mainz_OS$HR_p < 0.05) & (final_Mainz_OS$HR > 0), 'gene']$gene,
#                    Slovenia_GBM_OS = final_Slovenia_OS[(final_Slovenia_OS$log_rank < 0.05) & (final_Slovenia_OS$HR_p < 0.05) & (final_Slovenia_OS$HR > 0), 'gene']$gene)
# my_list_minus = list(TCGA_GBM_OS = final_TCGA_OS[(final_TCGA_OS$log_rank < 0.05) & (final_TCGA_OS$HR_p < 0.05) & (final_TCGA_OS$HR < 0), 'gene']$gene,
#                    Mainz_GBM_OS = final_Mainz_OS[(final_Mainz_OS$log_rank < 0.05) & (final_Mainz_OS$HR_p < 0.05) & (final_Mainz_OS$HR < 0), 'gene']$gene,
#                    Slovenia_GBM_OS = final_Slovenia_OS[(final_Slovenia_OS$log_rank < 0.05) & (final_Slovenia_OS$HR_p < 0.05) & (final_Slovenia_OS$HR < 0), 'gene']$gene)
my_list_plus = list(TCGA_GBM_OS = final_TCGA_OS[(final_TCGA_OS$log_rank < 0.05) & (final_TCGA_OS$HR_p < 0.05) & (final_TCGA_OS$HR > 0), 'gene'],
                  CGGA_GBM_OS = final_CGGA_OS[(final_CGGA_OS$log_rank < 0.05) & (final_CGGA_OS$HR_p < 0.05) & (final_CGGA_OS$HR > 0), 'gene'],
                   Slovenia_GBM_OS = final_Slovenia_OS[(final_Slovenia_OS$log_rank < 0.05) & (final_Slovenia_OS$HR_p < 0.05) & (final_Slovenia_OS$HR > 0), 'gene'])
my_list_minus = list(TCGA_GBM_OS = final_TCGA_OS[(final_TCGA_OS$log_rank < 0.05) & (final_TCGA_OS$HR_p < 0.05) & (final_TCGA_OS$HR < 0), 'gene'],
                   CGGA_GBM_OS = final_CGGA_OS[(final_CGGA_OS$log_rank < 0.05) & (final_CGGA_OS$HR_p < 0.05) & (final_CGGA_OS$HR < 0), 'gene'],
                   Slovenia_GBM_OS = final_Slovenia_OS[(final_Slovenia_OS$log_rank < 0.05) & (final_Slovenia_OS$HR_p < 0.05) & (final_Slovenia_OS$HR < 0), 'gene'])

my_list_plus <- list(TCGA_GBM_OS = final_TCGA_OS[(final_TCGA_OS$log_rank < 0.05) & (final_TCGA_OS$HR_p < 0.05) & (final_TCGA_OS$HR > 0), "pathways"],
                  CGGA_GBM_OS = final_CGGA_OS[(final_CGGA_OS$log_rank < 0.05) & (final_CGGA_OS$HR_p < 0.05) & (final_CGGA_OS$HR > 0), "pathways"],
                   Slovenia_GBM_OS = final_Slovenia_OS[(final_Slovenia_OS$log_rank < 0.05) & (final_Slovenia_OS$HR_p < 0.05) & (final_Slovenia_OS$HR > 0), "pathways"])
my_list_minus <- list(TCGA_GBM_OS = final_TCGA_OS[(final_TCGA_OS$log_rank < 0.05) & (final_TCGA_OS$HR_p < 0.05) & (final_TCGA_OS$HR < 0), "pathways"],
                   CGGA_GBM_OS = final_CGGA_OS[(final_CGGA_OS$log_rank < 0.05) & (final_CGGA_OS$HR_p < 0.05) & (final_CGGA_OS$HR < 0), "pathways"],
                   Slovenia_GBM_OS = final_Slovenia_OS[(final_Slovenia_OS$log_rank < 0.05) & (final_Slovenia_OS$HR_p < 0.05) & (final_Slovenia_OS$HR < 0), "pathways"])


my_list_all = list(TCGA = c(my_list_plus$TCGA_GBM_OS, my_list_minus$TCGA_GBM_OS),
                 CGGA = c(my_list_plus$CGGA_GBM_OS, my_list_minus$CGGA_GBM_OS),
                 Slovenia = c(my_list_plus$Slovenia_GBM_OS, my_list_minus$Slovenia_GBM_OS))

venn(x = my_list_plus, ilab = TRUE, zcolor = "style")
venn(x = my_list_minus, ilab = TRUE, zcolor = "style")
Reduce(intersect, my_list_plus)

cat(intersect(my_list_minus$TCGA_GBM_OS, my_list_minus$Mainz_GBM_OS), sep = "\n")
genes = c(intersect(my_list_minus$TCGA_GBM_OS, my_list_minus$Mainz_GBM_OS),
        intersect(my_list_plus$TCGA_GBM_OS, my_list_plus$Mainz_GBM_OS))


intersect(my_list_minus$TCGA_GBM_OS, my_list_minus$Mainz_GBM_OS)
write(intersect(my_list_plus$TCGA_GBM_OS, my_list_plus$CGGA_GBM_OS), "GBM_common_plus_pathways_TCGA_CGGA.txt")

final_TCGA_OS_short = final_TCGA_OS[!is.na(match(final_TCGA_OS$gene, our_genes)),]
final_CGGA_OS_short = final_CGGA_OS[!is.na(match(final_CGGA_OS$gene, our_genes)),]
final_Slovenia_OS_short = final_Slovenia_OS[!is.na(match(final_Slovenia_OS$gene, our_genes)),]

Reduce(intersect, my_list_plus)


test_randomness = function(genes_1, genes_2, number_1, number_2, number) {
  n = c()
  for (i in c(1:1000)) {
    n = c(n, length(intersect(sample(genes_1, number_1), sample(genes_2, number_2))))
    print(length(intersect(sample(genes_1, number_1), sample(genes_2, number_2))))
  }
  print(sum(n >= number) / length(n))
  return(n)
}

test_randomness_3 = function(genes_1, genes_2, genes_3, number_1, number_2, number_3, number) {
  n = c()
  for (i in c(1:1000)) {
    n = c(n, length(intersect(sample(genes_3, number_3), intersect(sample(genes_1, number_1), sample(genes_2, number_2)))))
    print(length(intersect(sample(genes_3, number_3), intersect(sample(genes_1, number_1), sample(genes_2, number_2)))))
  }
  return(sum(n >= number) / length(n))
}

test_randomness_3(genes_TCGA, genes_CGGA, genes_Slovenia, 1270, 1355, 958, 1)

a = test_randomness(genes_TCGA, genes_Mainz, 547, 332, 13)
a = test_randomness(pathways, pathways, 17, 14, 3)

length(final_TCGA_OS[(final_TCGA_OS$pvalue < 0.05) & (final_TCGA_OS$HR < 0), 'gene']$gene)

ggplot() + aes(as.numeric(c(dat2$BES))) + geom_histogram(binwidth = 2, colour = "black", fill = "white")



###########
### PCA ###
###########

TCGA = read_csv('counts_gbm_TCGA.csv')
CGGA = read_csv('counts_CGGA_693_GBM_OS.csv')
Mainz = read_csv('Mainz_counts_.csv')
# dat_oncobox = read_csv('C:/Users/raevs/Desktop/Oncobox/LTS_vs_STS_manuscript/CNV/oncobox_read_counts_ALL.csv')
dat_oncobox <- read_delim("C:/Users/raevs/Desktop/Oncobox/LTS_vs_STS_manuscript/exp_data/slovenia_exp_matrix.txt",
        "\t", escape_double = FALSE, trim_ws = TRUE)
p = dat_oncobox$SYMBOL
dat_oncobox = as.data.frame(dat_oncobox)
dat_oncobox$SYMBOL <- NULL
row.names(dat_oncobox) = p

metadata = read_excel('C:/Users/raevs/Desktop/Oncobox/LTS_vs_STS_manuscript/gliomas_pathways_PFS/our_glioblastoma.xlsx')
names = metadata[str_detect(metadata$`Sample source`, 'glioblastoma') | str_detect(metadata$`Sample source`, 'Glioblastoma'), 'RNAseq file name']$`RNAseq file name`
names = unique(na.omit(names))
names = str_sub(names, end = -1L - 9)
names = intersect(colnames(dat_oncobox), names)
names_slovenia = c('GB-1_S11_R1_001', 'GB_10_S5_R1_001', 'GB_11_S6_R1_001', 'GB_12_S7_R1_001', 'GB_16_S8_R1_001',
                 'GB_17_S9_R1_001', 'GB_18_S10_R1_001', 'GB_2_S4_R1_001', 'GB_3_S5_R1_001', 'GB_4_S6_R1_001', 'GB_5_S5_R1_001',
                 'GB_6_S9_R1_001', 'GB_7_S10_R1_001', 'GB_8_S4_R1_001')
# names = setdiff(names, names_slovenia)

dat_our = dat_oncobox[, colnames(dat_oncobox) %in% names]
dat_our$X1 = p

Slovenia = dat_oncobox[, colnames(dat_oncobox) %in% names_slovenia]
Slovenia$X1 = p

# all = merge(merge(merge(merge(TCGA, CGGA, by = 'X1'), Mainz, by = 'X1'), Slovenia, by = 'X1'), dat_our, by = 'X1')
all = merge(merge(TCGA, CGGA, by = 'X1'), Slovenia, by = 'X1')
#a=merge(merge(TCGA,CGGA,by='X1'),Mainz,by='X1')

p = all$X1
row.names(all) = p
all$X1 = NULL
# metadata = data.frame(name = colnames(all),
#                         dataset = c(rep('TCGA', ncol(TCGA) - 1),
#                                   rep('CGGA', ncol(CGGA) - 1),
#                                   rep('Mainz', ncol(Mainz) - 1),
#                                   rep('Slovenia', ncol(Slovenia) - 1),
#                                   rep('Oncobox', ncol(dat_our) - 1)))
# metadata <- metadata[!metadata$dataset %in% c("Mainz", "Oncobox"),]
metadata = data.frame(name = colnames(all),
                        dataset = c(rep('TCGA', ncol(TCGA) - 1),
                                  rep('CGGA', ncol(CGGA) - 1),
                                  rep('Slovenia', ncol(Slovenia) - 1)))

all = t(t(as.matrix(all)) / estimateSizeFactorsForMatrix(as.matrix(all)))
colnames(all) = paste0('tumour_', colnames(all))
all = all + 1
all = as.data.frame(all)
all$norm = c(10 ** (rowMeans(log10(all))))
all$SYMBOL = row.names(all)
all = all[, c("SYMBOL", colnames(all)[1:(ncol(all) - 1)])]
write.csv(all, 'all_counts_primary_GBM_OS_pathways.csv', row.names = FALSE)

all_normalized_quantiles = as.data.frame(normalize.quantiles(as.matrix(all[, !colnames(all) %in% c("SYMBOL")] + 1), copy = FALSE))

all_pathways = read_excel('all.xlsx', sheet = 2)
all_pathways$Database = NULL

our_pathways = read.xlsx("C:/Users/raevs/Desktop/Oncobox/LTS_vs_STS_manuscript/gliomas_pathways_PFS/number_of_genes_in_pathways.xlsx")
all_pathways = all_pathways[!is.na(match(all_pathways$Pathway, our_pathways$pathway)),]
p = all_pathways$Pathway
all_pathways$Pathway = NULL
all_pathways <- as.data.frame(all_pathways)
row.names(all_pathways) = p
colnames(all_pathways) = str_sub(colnames(all_pathways), 8, -1L)

# library(pca3d)
# a = all_normalized_quantiles[intersect(our_genes, row.names(all_normalized_quantiles)),]
# pca <- prcomp(t(log10(all_normalized_quantiles)))
# pca <- prcomp(t(all_pathways))


# Run PCA - Gene Expression
all_normalized_quantiles <- all_normalized_quantiles[, intersect(colnames(all_normalized_quantiles), paste0("tumour_", metadata$name))]
pca_plot_by_gene_exp <- plot_pca(data_for_pca = t(log10(all_normalized_quantiles)))

export_analysis_plot(filename = "pca_plot_by_gene_exp",
                     plot = pca_plot_by_gene_exp,
                     path = "C:/Users/raevs/Desktop/Oncobox/LTS_vs_STS_manuscript/plots/pca",
                     scale = 1,
                     width = 210 / 2,
                     height = 297 / 2.75,
                     units = "mm",
                     load_fonts = FALSE, cairo_pdf_device = TRUE)

# Run PCA - Pathways Scores
all_pathways <- all_pathways[, intersect(colnames(all_pathways), metadata$name)]
pca_plot_by_pathways_scores <- plot_pca(data_for_pca = t(all_pathways))

export_analysis_plot(filename = "pca_plot_by_pathways_scores",
                     plot = pca_plot_by_pathways_scores + xlim(-100, 100) + ylim(-100, 100),
                     path = "C:/Users/raevs/Desktop/Oncobox/LTS_vs_STS_manuscript/plots/pca",
                     scale = 1,
                     width = 210 / 2,
                     height = 297 / 2.75,
                     units = "mm",
                     load_fonts = FALSE, cairo_pdf_device = TRUE)

#png("pca.png")
# pca3d(pca, group = factor(metadata$dataset), legend = "topleft", palette = c("#F0A0FF", "#0075DC", "#993F00", "#4C005C", "#191919"), shape = NULL)

# pdf("pca_plot_pathways.pdf", width = 10, height = 50)
# pca2d(pca, group = factor(metadata$dataset), legend = "topright", palette = c("#F0A0FF", "#0075DC", "#993F00", "#2BCE48", "#191919"), shape = NULL)
# dev.off()

# library(rgl)
# plot3d(pca$scores[, 1:3], col = c(1)))


# "#005C31", "#2BCE48", "#FFCC99", "#808080", "#94FFB5",
# "#8F7C00", "#9DCC00", "#C20088", "#003380", "#FFA405",
# "#FFA8BB", "#426600", "#FF0010", "#5EF1F2", "#00998F",
# "#E0FF66", "#740AFF", "#990000", "#FFFF80", "#FFE100",
# "#FF5005"

#### Volcano plots and GSEA viz ###

library(enrichplot)
library(EnhancedVolcano)

final_TCGA_OS = read_excel('GBM_PFS_HRs_genes_third.xlsx')
final_TCGA_OS$log10_log_rank = log10(final_TCGA_OS$log_rank)
final_TCGA_OS = na.omit(final_TCGA_OS)
# final_Mainz_OS = read_excel('Mainz_PFS_HRs_genes_third.xlsx')
# final_Mainz_OS$log10_log_rank = log10(final_Mainz_OS$log_rank)
# final_Mainz_OS = na.omit(final_Mainz_OS)
final_Slovenia_OS = read_excel('Slovenia_primary_GBM_PFS_HRs_genes_third.xlsx')
final_Slovenia_OS$log10_log_rank = log10(final_Slovenia_OS$log_rank)
final_Slovenia_OS = na.omit(final_Slovenia_OS)
final_Slovenia_OS <- as.data.frame(final_Slovenia_OS)

EnhancedVolcano(final_Slovenia_OS,
                x = 'HR',
                y = 'log_rank',
                lab = final_Slovenia_OS$gene,
                pCutoff = 0.05, FCcutoff = 30)

# pdf("Volcano_Slovenia_PFS_short.pdf", height = 8, width = 8)
# EnhancedVolcano(final_Slovenia_OS,
#                 x = 'HR',
#                 y = 'log_rank',
#                 lab = final_Slovenia_OS$gene,
#                 pCutoff = 0.05, FCcutoff = 30, ylim = c(0, 3.5), xlim = c(-3, 3))
# dev.off()

library("openxlsx")
library(clusterProfiler)
library(enrichplot)

genes_CGGA_TCGA = intersect(my_list_plus$TCGA_GBM_OS, my_list_plus$CGGA_GBM_OS)
genes_CGGA_Slovenia = intersect(my_list_plus$Slovenia_GBM_OS, my_list_plus$CGGA_GBM_OS)
genes_TCGA_Mainz = intersect(my_list_minus$TCGA_GBM_OS, my_list_minus$Mainz_GBM_OS)

genes_TCGA_Mainz = c(intersect(my_list_minus$TCGA_GBM_OS, my_list_minus$Mainz_GBM_OS),
 intersect(my_list_plus$TCGA_GBM_OS, my_list_plus$Mainz_GBM_OS))

entrez = as.character(hgnc[!is.na(match(hgnc$symbol, genes_TCGA_Mainz)), "entrez_id"]$entrez_id)
ego2 <- enrichGO(entrez, OrgDb = "org.Hs.eg.db", ont = "ALL", readable = TRUE)


pdf("enrich_plot_TCGA_Mainz_PFS.pdf", height = 6, width = 9)
dotplot(ego1, showCategory = 100)
dev.off()

pdf("enrich_net_CGGA_TCGA.pdf", height = 5, width = 5)
emapplot(ego1, showCategory = 100, layout = 'circle', min_edge = 1)
dev.off()

cat(genes_CGGA_TCGA, sep = "\n")



genes_CGGA_TCGA_OS <- c("STC1", "CD248", "ESM1", "COL5A1", "LRRFIP1", "TMEM110", "COL10A1", "MAP2K3", "MMP11",
"EN2", "PLAUR", "CD276", "FAM20C", "SERPINA5", "SLC43A3", "KDELR2", "ANO6", "HSP90B1",
"MARVELD1", "IL4I1", "CLEC5A", "CUL7", "THBD", "BACE2", "ESYT1", "ACTA2", "MYO1B",
"LBH", "CCL2", "AEBP1", "MXRA5", "LAX1", "DNAJC22", "ARMC10", "SLC35C1",
"RGS16", "LOXL1", "STEAP3", "PDAP1", "MCAM", "RGS3", "ITGA3", "SPHK1",
"TMEM2", "PDIA4", "FN1", "CKAP4", "COL22A1", "OSBPL5", "C1QTNF6", "MYO1C",
"RELT", "PXDN", "CNPY4", "SOCS3", "LUCAT1", "MPZL3", "IKBIP", "CRNDE")

entrez = as.character(hgnc[!is.na(match(hgnc$symbol, genes_CGGA_TCGA_OS)), "entrez_id"]$entrez_id)
go_enrichment <- enrichGO(entrez, OrgDb = "org.Hs.eg.db", ont = "ALL", readable = TRUE)