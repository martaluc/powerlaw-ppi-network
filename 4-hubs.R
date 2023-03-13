source('functions.R')
library(data.table)
library(formattable)
library(clusterProfiler) # version 4.4.4


#---------------------------------------------
# select preys in mass spec studies
#---------------------------------------------
intact <- read.csv('databases/IntAct_afterFiltering.csv')
prey_table <- get_prey_massSpec(intact)


#-------------------------------------------------
# select Huri study = 32296183 (yeast two hybrid)
#-------------------------------------------------
s <- '32296183'
degree_32296183 <- calculate_degree_y2h_study(intact,s)

#-------------------------------------------------------------
# select hubs: ordering the degree
#--------------------------------------------------------------

dir_out <- 'output'

# prey-based hubs
prey_table <- prey_table[order(-as.numeric(prey_table$degree_prey)),]
write.csv(prey_table[,c(1,2,4)], file = paste0(dir_out,'/prey_hubs.csv'),row.names = F)

# hubs Y2h study
degree_32296183 <- degree_32296183[order(-degree_32296183$degree),]
write.csv(degree_32296183, file = paste0(dir_out,'/Y2H_hubs.csv'), row.names = F)

# hubs normalized bait degree
bait_table <- read.csv('output/bait_usage_intact2022.csv')
hippie_intact <- read.csv('databases/HIPPIE_union_Intact2022_afterReviewed_mapping.csv')
g <- unique(hippie_intact[,c('IDs_interactor_A','IDs_interactor_B')])
degree_hippie_intact <- degree_wo_bidirEdges(g)
degree_hippie_intact$proteins <- rownames(degree_hippie_intact)
bait_table$degree <- degree_hippie_intact$degree[match(bait_table$bait_uniprot,degree_hippie_intact$proteins)]
bait_table$normaliz_degree <- bait_table$degree/bait_table$bait_usage
bait_table <- bait_table[order(-bait_table$normaliz_degree),]

write.csv(bait_table, file = paste0(dir_out,'/hubs_normal_bait.csv'), row.names = F)

# hubs of HIPPIEunionIntact
# hippie_intact <- read.csv('databases/HIPPIE_union_Intact2022_afterReviewed_mapping.csv')
# g <- unique(hippie_intact[,c('IDs_interactor_A','IDs_interactor_B')])
# degree_hippie_intact <- degree_wo_bidirEdges(g)
# degree_hippie_intact$proteins <- rownames(degree_hippie_intact)
degree_hippie_intact <- degree_hippie_intact[order(-degree_hippie_intact$degree),]
symbol <- as.data.frame(mapIds(org.Hs.eg.db, keys = degree_hippie_intact$proteins, keytype = "UNIPROT", column= "SYMBOL"))
colnames(symbol)[1] <- 'symbol'
degree_hippie_intact$symbol <- symbol$symbol[match(degree_hippie_intact$proteins,rownames(symbol))]
write.csv(degree_hippie_intact, file = paste0(dir_out,'/AggregatedNetwork_hubs.csv'), row.names = F)

#------------------------------------------------------------------------------------
# clusterProfile - enrichment analyses
#------------------------------------------------------------------------------------

if(!dir.exists('output/enrichDO_results')){
  dir.create('output/enrichDO_results')
}
# if(!dir.exists('output/enrichKEGG_results')){
#   dir.create('output/enrichKEGG_results')
# }
if(!dir.exists('output/enrichPathway_results')){
  dir.create('output/enrichPathway_results')
}
if(!dir.exists('output/enrichGO_results')){
  dir.create('output/enrichGO_results')
}

prey_table <- read.csv('output/prey_hubs.csv')
bait_table <- read.csv('output/hubs_normal_bait.csv')
y2h <- read.csv('output/Y2H_hubs.csv')
agg_net <- read.csv('output/AggregatedNetwork_hubs.csv')

entrez_prey <- as.data.frame(mapIds(org.Hs.eg.db, keys = prey_table$prey_uniprot, keytype = 'UNIPROT', column= "ENTREZID"))
colnames(entrez_prey)[1] <- 'entrezID'
entrez_bait <- as.data.frame(mapIds(org.Hs.eg.db, keys = bait_table$bait_uniprot, keytype = 'UNIPROT', column= "ENTREZID"))
colnames(entrez_bait)[1] <- 'entrezID'
entrez_y2h <- as.data.frame(mapIds(org.Hs.eg.db, keys = y2h$uniprot, keytype = 'UNIPROT', column= "ENTREZID"))
colnames(entrez_y2h)[1] <- 'entrezID'
entrez_agg <- as.data.frame(mapIds(org.Hs.eg.db, keys = agg_net$proteins, keytype = 'UNIPROT', column= "ENTREZID"))
colnames(entrez_agg)[1] <- 'entrezID'

n <- 50 # number of top hubs
clusterProfiler_analyses(entrez_prey$entrezID[seq(1,n)],unique(entrez_prey$entrezID),p = 0.05, n,'prey','BP','output')
clusterProfiler_analyses(entrez_bait$entrezID[seq(1,n)],unique(entrez_bait$entrezID),p = 0.05, n,'normalized_hubs','BP','output')
clusterProfiler_analyses(entrez_y2h$entrezID[seq(1,n)],unique(entrez_y2h$entrezID),p = 0.05, n,'Y2H_hubs','BP','output')
clusterProfiler_analyses(entrez_agg$entrezID[seq(1,n)],unique(entrez_agg$entrezID),p = 0.05, n,'Aggregated_network','BP','output')

# run to know how many genes are annotated in databases for those sets where no significant enrichment has been found 
# clusterProfiler_analyses(entrez_y2h$entrezID[seq(1,n)],unique(entrez_y2h$entrezID),p = 1, n,'Y2H_hubs','BP','output')
# clusterProfiler_analyses(entrez_bait$entrezID[seq(1,n)],unique(entrez_bait$entrezID),p = 1, n,'normalized_hubs','BP','output')
#---------------------------------------------  
# plot the results of the enrichment analyses
#---------------------------------------------

dirs <- c('output/enrichGO_results','output/enrichDO_results','output/enrichPathway_results')
for(d in dirs){
  print(d)
  n <- 50
  p <- 0.05
  if(file.exists(paste0(d,'/table_',unlist(strsplit(strsplit(d, split = '/')[[1]][2],split = '_'))[1],'_prey_',n,'_qvalue',p,'.csv'))){
    go_prey <- read.csv(paste0(d,'/table_',unlist(strsplit(strsplit(d, split = '/')[[1]][2],split = '_'))[1],'_prey_',n,'_qvalue',p,'.csv'))
    go_prey$Cluster <- rep('Prey hubs',nrow(go_prey))
  }
  if(file.exists(paste0(d,'/table_',unlist(strsplit(strsplit(d, split = '/')[[1]][2],split = '_'))[1],'_normalized_hubs_',n,'_qvalue',p,'.csv'))){
    go_normal <- read.csv(paste0(d,'/table_',unlist(strsplit(strsplit(d, split = '/')[[1]][2],split = '_'))[1],'_normalized_hubs_',n,'_qvalue',p,'.csv'))
    go_normal$Cluster <- rep('Normalized hubs',nrow(go_normal))
  }
  if(file.exists(paste0(d,'/table_',unlist(strsplit(strsplit(d, split = '/')[[1]][2],split = '_'))[1],'_Y2H_hubs_',n,'_qvalue',p,'.csv'))){
    go_y2h <- read.csv(paste0(d,'/table_',unlist(strsplit(strsplit(d, split = '/')[[1]][2],split = '_'))[1],'_Y2H_hubs_',n,'_qvalue',p,'.csv'))
    go_y2h$Cluster <- rep('Y2H hubs',nrow(go_y2h)) 
  }
  if(file.exists(paste0(d,'/table_',unlist(strsplit(strsplit(d, split = '/')[[1]][2],split = '_'))[1],'_Aggregated_network_',n,'_qvalue',p,'.csv'))){
    go_agg <- read.csv(paste0(d,'/table_',unlist(strsplit(strsplit(d, split = '/')[[1]][2],split = '_'))[1],'_Aggregated_network_',n,'_qvalue',p,'.csv'))
    go_agg$Cluster <- rep('Uncorrected aggregated \n network hubs',nrow(go_agg))
  }
  
  if(exists('go_y2h') & exists('go_normal')){
    t <- rbind(go_prey,go_normal,go_y2h,go_agg)
  }else{
    if(!exists('go_y2h') & exists('go_normal')){
      t <- rbind(go_prey,go_normal,go_agg)
    }else{
      if(!exists('go_y2h') & !exists('go_normal')){
        t <- rbind(go_prey,go_agg)
      }else{
        if(exists('go_y2h') & !exists('go_normal')){
          t <- rbind(go_prey,go_y2h,go_agg)
        }
      }
    }
  }
  
  t$Cluster <- factor(t$Cluster, levels = c('Uncorrected aggregated \n network hubs','Prey hubs','Normalized hubs','Y2H hubs'))
  t$Description <-sub("(.)", "\\U\\1",t$Description,perl=TRUE)
  
  res <- new("compareClusterResult", compareClusterResult = t, .call = match.call(expand.dots = TRUE))
  dotplot(res, font.size = 7, label_format = 30, color = 'qvalue') + theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=0.5))
  ggsave(filename = paste0(d,'/dotplot_',unlist(strsplit(strsplit(d, split = '/')[[1]][2],split = '_'))[1],'_n',n,'_qvalue',p,'.pdf'),height = 15, width = 10,units = 'cm')
  rm(list = ls(all.names = TRUE))
}


#----------------------------------------------------------------------
# enrichment between schizophrenia/psychotic disorder and chaperones
#----------------------------------------------------------------------
library(DOSE)
.initial <- function() {
  pos <- 1
  envir <- as.environment(pos)
  assign(".DOSEEnv", new.env(), envir = envir)
  .DOSEEnv <- get(".DOSEEnv", envir = envir)
  
  tryCatch(utils::data(list="dotbl",
                       package="DOSE"))
  dotbl <- get("dotbl")
  assign("dotbl", dotbl, envir = .DOSEEnv)
  rm(dotbl, envir = .GlobalEnv)
  
  tryCatch(utils::data(list="DOIC",
                       package="DOSE"))
  DOIC <- get("DOIC")
  assign("DOIC", DOIC, envir = .DOSEEnv)
  rm(DOIC, envir = .GlobalEnv)
  
}

uniprot_reviewed <- read.csv('databases/uniprot-compressed_true_download_true_fields_accession_2Creviewed_2C-2022.12.13-09.55.05.21.tsv', sep = '\t')
all_genes <- uniprot_reviewed$Entry
all_symbol <- mapIds(org.Hs.eg.db, keys = all_genes, keytype = 'UNIPROT', column= "SYMBOL")
all_symbol <- all_symbol[-which(is.na(all_symbol) == T)]
all_symbol <- unique(all_symbol)
length(all_symbol)

#retrieve DO term genes
if (!exists(".DOSEEnv")) {
  .initial()
}
DOSEEnv <- get(".DOSEEnv", envir = .GlobalEnv)
if (!exists("DO2ALLEG", envir=DOSEEnv)) {
  tryCatch(utils::data(list="DO2ALLEG", package="DOSE"))
  assign("DO2ALLEG", DO2ALLEG, envir = DOSEEnv)
  DO2ALLEG <- get("DO2ALLEG")
  rm(DO2ALLEG, envir = .GlobalEnv)
}
PATHID2EXTID <- get("DO2ALLEG", envir = DOSEEnv)

term <- c('Schizophrenia','Psychotic disorder')

for(t in term){
  print(t)
  if(t == 'Schizophrenia'){
    genes <- PATHID2EXTID$'DOID:5419'
  }
  if(t == 'Psychotic disorder'){
    genes <- PATHID2EXTID$'DOID:2468'
  }
  
  gene_symbol <- mapIds(org.Hs.eg.db, keys = genes, keytype = 'ENTREZID', column= "SYMBOL")
  gene_symbol <- gene_symbol[which(gene_symbol %in% all_symbol)]
  length(gene_symbol)
  
  chaperone <- read.csv('databases/chaperone_uniprot_function.csv')
  chaperone <- unique(chaperone$Entry)
  chaperone <- mapIds(org.Hs.eg.db, keys = chaperone, keytype = 'UNIPROT', column= "SYMBOL")
  chaperone <- unique(chaperone[which(chaperone %in% all_symbol)])
  length(chaperone)
  
  library(ggvenn)
  l <- list(gene_symbol,chaperone)
  names(l) <- c(t,'Chaperones')
  ggvenn(l, show_percentage = F, text_size = 4, set_name_size = 3)
  ggsave(paste0('plots/VennDiagram_',t,'-chaperone.pdf'),units = 'cm', height = 8,width = 10,device = 'pdf')
  print(sort(intersect(chaperone,gene_symbol)))
  
  #fisher test
  disease <- ifelse(all_symbol %in% gene_symbol,'yes','no')
  chaperone <-  ifelse(all_symbol %in% chaperone,'yes','no')
  d <- data.frame(all_symbol,disease,chaperone)
  d$disease <- factor(d$disease,levels = c("yes", "no"))
  d$chaperone <- factor(d$chaperone,levels = c("yes", "no"))
  m <- table(d$chaperone,d$disease)
  print(m)
  f <- fisher.test(m)
  f_greater <- fisher.test(m, alternative = 'greater')
  #print(f$p.value)
  print('Fisher test p-value')
  print(f_greater$p.value)
}

#-------------------------------------------------------------------------
# GTEx abundance - protein level
# performing the GO analysis for the most abundant proteins
#--------------------------------------------------------------------------

library(readxl)
# 'Table_S1_gene_info_at_protein_level.xlsx' has been downloaded by https://gtexportal.org/home/datasets (eGTEx section)
gtex_mapping <- read_excel('databases/Table_S1_gene_info_at_protein_level.xlsx', sheet = 11)

gtex_protein_relative_abun <- read_excel('databases/Table_S1_gene_info_at_protein_level.xlsx', sheet = 8)
colnames(gtex_protein_relative_abun) <- gtex_protein_relative_abun[1,]
gtex_protein_relative_abun <- gtex_protein_relative_abun[-1,] 

# calculate the median abundance for each protein
Median <- c()
Mean <- c()
n_NA <- c()
for(i in seq(1,nrow(gtex_protein_relative_abun))){
  row <- as.numeric(gtex_protein_relative_abun[i,seq(3,ncol(gtex_protein_relative_abun))])
  n_NA <- c(n_NA,length(which(is.na(row) == T))/length(row))
  Median <- c(Median,median(row, na.rm = T))
  Mean <- c(Mean,mean(row,na.rm = T))
}

length(Median)
dim(gtex_protein_relative_abun)
table_median_protein <- data.frame(gtex_protein_relative_abun$gene.id,Median,Mean,n_NA)
colnames(table_median_protein)[1] <- 'Gene_ID'
table_median_protein$entrez_id <- gtex_mapping$entrez_id[match(table_median_protein$Gene_ID,gtex_mapping$ensembl_id)]
# remove genes with more than 50% of NAs
table_median_protein <- table_median_protein[-which(table_median_protein$n_NA > 0.5),]
table_median_protein <- table_median_protein[order(-table_median_protein$Median),]
table_median_protein <- table_median_protein[-which(table_median_protein$entrez_id == 'NA'),]

# perform GO analysis for different set size
for(n in c(200,500,1000,2000,3000)){
  ego <- enrichGO(table_median_protein$entrez_id[seq(1,n)],universe=table_median_protein$entrez_id,OrgDb= org.Hs.eg.db,ont= 'BP',pAdjustMethod = "fdr",qvalueCutoff  = 0.05,pvalueCutoff =0.05,readable= TRUE,minGSSize = 2)
  write.csv(as.data.frame(ego),file = paste0('output/GO_Gtex_protein_relativeAbund',n,'.csv'), row.names = F)
}

#plot GO terms
go_200 <- read.csv('output/GO_Gtex_protein_relativeAbund200.csv')
go_200$Cluster <- rep('top-200',nrow(go_200)) 
go_500 <- read.csv('output/GO_Gtex_protein_relativeAbund500.csv')
go_500$Cluster <- rep('top-500',nrow(go_500))
go_1000 <- read.csv('output/GO_Gtex_protein_relativeAbund1000.csv')
go_1000$Cluster <- rep('top-1000',nrow(go_1000))
go_2000 <- read.csv('output/GO_Gtex_protein_relativeAbund2000.csv')
go_2000$Cluster <- rep('top-2000',nrow(go_2000))
go_3000 <- read.csv('output/GO_Gtex_protein_relativeAbund3000.csv')
go_3000$Cluster <- rep('top-3000',nrow(go_3000))
go_final <- rbind(go_200,go_500,go_1000,go_2000,go_3000)
go_final$qvalue <- as.numeric(go_final$qvalue)
go_final$Cluster <- factor(go_final$Cluster, levels = c('top-200','top-500','top-1000','top-2000','top-3000'))
res <- new("compareClusterResult", compareClusterResult = go_final, .call = match.call(expand.dots = TRUE))
dotplot(res, color = 'qvalue',font.size = 10, label_format = 35) + theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=0.5))
ggsave(filename = 'plots/dotplot_GO_Gtex_protein_relativeAbund.pdf',height = 18, width = 16,units = 'cm')

