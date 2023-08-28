source('functions.R')

intact <- read.csv('databases/IntAct_afterFiltering.csv')
#intact with bait/prey annotations
i <- which(intact$Experimental_roles_interactor_A == 'bait' & intact$Experimental_roles_interactor_B == 'prey')
j <- which(intact$Experimental_roles_interactor_B == 'bait' & intact$Experimental_roles_interactor_A == 'prey')
intact <- intact[union(i,j),]

#---------------------
# AP-MS
#----------------------
#select AP_MS studies
index <- grep('coimmunoprecipitation|pull down|tandem affinity purification|mass spectrometry|copurification|affinity chromatography technology',intact$Interaction_detection_methods)
intact_MS <- intact[index,]
#calculate degree
g <- unique(intact_MS[,c('IDs_interactor_A','IDs_interactor_B')])
degree_intact_MS <- degree_wo_bidirEdges(g)
degree_intact_MS$proteins <- rownames(degree_intact_MS)
final_MS <- get_bait_prey_total_degree(intact_MS,degree_intact_MS)
write.csv(final_MS,'output/table_bait_prey_total_degree_AP-MS.csv',row.names = F)

# test PL property of degree_bait, degree_prey and total_degree
p_bait <- check_powerLaw(as.numeric(final_MS$degree_bait[-which(final_MS$degree_bait == 0)]), plot = T, plot_name = 'plots/plot_bait_degree_AP_MS', t = 20, 'Degree')
p_prey <- check_powerLaw(as.numeric(final_MS$degree_prey[-which(final_MS$degree_prey == 0)]), plot = T, plot_name = 'plots/plot_prey_degree_AP_MS', t = 20, 'Degree')
# Supplementary Figure 5A
p_total <- check_powerLaw(as.numeric(final_MS$total_degree), plot = T, plot_name = 'plots/plot_total_degree_AP_MS', t = 20, 'Degree')

pvalue <- c(p_bait$p,p_prey$p,p_total$p)
xmin <- c(xmin_estimated(as.numeric(final_MS$degree_bait[-which(final_MS$degree_bait == 0)])),xmin_estimated(as.numeric(final_MS$degree_prey[-which(final_MS$degree_prey == 0)])),
          xmin_estimated(as.numeric(final_MS$total_degree)))
alpha <- c(alpha_estimated(as.numeric(final_MS$degree_bait[-which(final_MS$degree_bait == 0)])),alpha_estimated(as.numeric(final_MS$degree_prey[-which(final_MS$degree_prey == 0)])),
           alpha_estimated(as.numeric(final_MS$total_degree)))

PL_table_MS <- data.frame(rep('AP_MS',3),c('bait_degree','prey_degree','total_degree'),pvalue,xmin,alpha)
colnames(PL_table_MS)[c(1,2)] <- c('method','type')

#----------------
# Y2H
#---------------
# select Y2H studies
intact_y2h <- intact[grep('two hybrid',intact$Interaction_detection_methods),]
g <- unique(intact_y2h[,c('IDs_interactor_A','IDs_interactor_B')])
# calculate degree
degree_intact_y2h <- degree_wo_bidirEdges(g)
degree_intact_y2h$proteins <- rownames(degree_intact_y2h)
final_y2h <- get_bait_prey_total_degree(intact_y2h,degree_intact_y2h)
write.csv(final_y2h,'output/table_bait_prey_total_degree_Y2H.csv',row.names = F)

# test PL property of degree_bait, degree_prey and total_degree
p_bait <- check_powerLaw(as.numeric(final_y2h$degree_bait[-which(final_y2h$degree_bait == 0)]), plot = T, plot_name = 'plots/plot_bait_degree_Y2H', t = 20, 'Degree')
p_prey <- check_powerLaw(as.numeric(final_y2h$degree_prey[-which(final_y2h$degree_prey == 0)]), plot = T, plot_name = 'plots/plot_prey_degree_Y2H', t = 20, 'Degree')
# Supplementary Figure 5B
p_total <- check_powerLaw(as.numeric(final_y2h$total_degree), plot = T, plot_name = 'plots/plot_total_degree_Y2H', t = 20, 'Degree')

pvalue <- c(p_bait$p,p_prey$p,p_total$p)
xmin <- c(xmin_estimated(as.numeric(final_y2h$degree_bait[-which(final_y2h$degree_bait == 0)])),xmin_estimated(as.numeric(final_y2h$degree_prey[-which(final_y2h$degree_prey == 0)])),
          xmin_estimated(as.numeric(final_y2h$total_degree)))
alpha <- c(alpha_estimated(as.numeric(final_y2h$degree_bait[-which(final_y2h$degree_bait == 0)])),alpha_estimated(as.numeric(final_y2h$degree_prey[-which(final_y2h$degree_prey == 0)])),
           alpha_estimated(as.numeric(final_y2h$total_degree)))

PL_table_Y2H <- data.frame(rep('Y2H',3),c('bait_degree','prey_degree','total_degree'),pvalue,xmin,alpha)
colnames(PL_table_Y2H)[c(1,2)] <- c('method','type')

# union of the two final table about PL
final_table <- rbind(PL_table_MS,PL_table_Y2H)
write.csv(final_table,'output/table_PL_pvalue_AP-MS_Y2H.csv',row.names = F)

#-------------------------------------------------------------------------------
# same analysis but filtering studies with more than 100 interactions
#------------------------------------------------------------------------------

# intact <- read.csv('databases/IntAct_afterFiltering.csv')
# #intact with bait/prey annotations
# i <- which(intact$Experimental_roles_interactor_A == 'bait' & intact$Experimental_roles_interactor_B == 'prey')
# j <- which(intact$Experimental_roles_interactor_B == 'bait' & intact$Experimental_roles_interactor_A == 'prey')
# intact <- intact[union(i,j),]
# 
# # select studies with more than 100 interactions
# pubmedIDs <- unique(intact$Publication_Identifiers)
# length(pubmedIDs) # 7599
# table_num_inter <- get_numInter_studies(intact,pubmedIDs)
# table_num_inter <- as.data.frame(table_num_inter)
# table_num_inter$num_inter <- as.numeric(table_num_inter$num_inter)
# n <- 100
# subset_pubmed <- table_num_inter$pubmedID[which(table_num_inter$num_inter >= n)]
# length(subset_pubmed) # 204
# 
# intact <- intact[which(intact$Publication_Identifiers %in% subset_pubmed),]
# 
# #---------------------
# # AP-MS
# #----------------------
# #select AP_MS studies
# index <- grep('coimmunoprecipitation|pull down|tandem affinity purification|mass spectrometry|copurification|affinity chromatography technology',intact$Interaction_detection_methods)
# intact_MS <- intact[index,]
# length(unique(intact_MS$Publication_Identifiers))#176
# #calculate degree
# g <- unique(intact_MS[,c('IDs_interactor_A','IDs_interactor_B')])
# degree_intact_MS <- degree_wo_bidirEdges(g)
# degree_intact_MS$proteins <- rownames(degree_intact_MS)
# final_MS <- get_bait_prey_total_degree(intact_MS,degree_intact_MS)
# write.csv(final_MS,paste0('output/table_bait_prey_total_degree_AP-MS_',n,'.csv'),row.names = F)
# 
# # test PL property of degree_bait, degree_prey and total_degree
# p_bait <- check_powerLaw(as.numeric(final_MS$degree_bait[-which(final_MS$degree_bait == 0)]), plot = T, plot_name = paste0('plots/plot_bait_degree_AP_MS_',n), t = 20, 'Degree')
# p_prey <- check_powerLaw(as.numeric(final_MS$degree_prey[-which(final_MS$degree_prey == 0)]), plot = T, plot_name = paste0('plots/plot_prey_degree_AP_MS_',n), t = 20, 'Degree')
# p_total <- check_powerLaw(as.numeric(final_MS$total_degree), plot = T, plot_name = paste0('plots/plot_total_degree_AP_MS_',n), t = 20, 'Degree')
# 
# pvalue <- c(p_bait$p,p_prey$p,p_total$p)
# xmin <- c(xmin_estimated(as.numeric(final_MS$degree_bait[-which(final_MS$degree_bait == 0)])),xmin_estimated(as.numeric(final_MS$degree_prey[-which(final_MS$degree_prey == 0)])),
#           xmin_estimated(as.numeric(final_MS$total_degree)))
# alpha <- c(alpha_estimated(as.numeric(final_MS$degree_bait[-which(final_MS$degree_bait == 0)])),alpha_estimated(as.numeric(final_MS$degree_prey[-which(final_MS$degree_prey == 0)])),
#            alpha_estimated(as.numeric(final_MS$total_degree)))
# 
# PL_table_MS <- data.frame(rep('AP_MS',3),c('bait_degree','prey_degree','total_degree'),pvalue,xmin,alpha)
# colnames(PL_table_MS)[c(1,2)] <- c('method','type')
# 
# #----------------
# # Y2H
# #---------------
# # select Y2H studies
# intact_y2h <- intact[grep('two hybrid',intact$Interaction_detection_methods),]
# length(unique(intact_y2h$Publication_Identifiers)) #53
# g <- unique(intact_y2h[,c('IDs_interactor_A','IDs_interactor_B')])
# # calculate degree
# degree_intact_y2h <- degree_wo_bidirEdges(g)
# degree_intact_y2h$proteins <- rownames(degree_intact_y2h)
# final_y2h <- get_bait_prey_total_degree(intact_y2h,degree_intact_y2h)
# write.csv(final_y2h,paste0('output/table_bait_prey_total_degree_Y2H_',n,'.csv'),row.names = F)
# 
# # test PL property of degree_bait, degree_prey and total_degree
# p_bait <- check_powerLaw(as.numeric(final_y2h$degree_bait[-which(final_y2h$degree_bait == 0)]), plot = T, plot_name = paste0('plots/plot_bait_degree_Y2H_',n), t = 20, 'Degree')
# p_prey <- check_powerLaw(as.numeric(final_y2h$degree_prey[-which(final_y2h$degree_prey == 0)]), plot = T, plot_name = paste0('plots/plot_prey_degree_Y2H_',n), t = 20, 'Degree')
# p_total <- check_powerLaw(as.numeric(final_y2h$total_degree), plot = T, plot_name = paste0('plots/plot_total_degree_Y2H_',n), t = 20, 'Degree')
# 
# pvalue <- c(p_bait$p,p_prey$p,p_total$p)
# xmin <- c(xmin_estimated(as.numeric(final_y2h$degree_bait[-which(final_y2h$degree_bait == 0)])),xmin_estimated(as.numeric(final_y2h$degree_prey[-which(final_y2h$degree_prey == 0)])),
#           xmin_estimated(as.numeric(final_y2h$total_degree)))
# alpha <- c(alpha_estimated(as.numeric(final_y2h$degree_bait[-which(final_y2h$degree_bait == 0)])),alpha_estimated(as.numeric(final_y2h$degree_prey[-which(final_y2h$degree_prey == 0)])),
#            alpha_estimated(as.numeric(final_y2h$total_degree)))
# 
# PL_table_Y2H <- data.frame(rep('Y2H',3),c('bait_degree','prey_degree','total_degree'),pvalue,xmin,alpha)
# colnames(PL_table_Y2H)[c(1,2)] <- c('method','type')
# 
# # union of the two final table about PL
# final_table <- rbind(PL_table_MS,PL_table_Y2H)
# write.csv(final_table,paste0('output/table_PL_pvalue_AP-MS_Y2H_',n,'.csv'),row.names = F)

#------------------------------------------------------------------------------
# AP-MS and Y2H combined
#------------------------------------------------------------------------------

intact <- read.csv('databases/IntAct_afterFiltering.csv')
#intact with bait/prey annotations
i <- which(intact$Experimental_roles_interactor_A == 'bait' & intact$Experimental_roles_interactor_B == 'prey')
j <- which(intact$Experimental_roles_interactor_B == 'bait' & intact$Experimental_roles_interactor_A == 'prey')
intact <- intact[union(i,j),]
#select AP-MS studies
index1 <- grep('coimmunoprecipitation|pull down|tandem affinity purification|mass spectrometry|copurification|affinity chromatography technology',intact$Interaction_detection_methods)
#select Y2H studies
index2 <- grep('two hybrid',intact$Interaction_detection_methods)
intact <- intact[union(index1,index2),]

g <- unique(intact[,c('IDs_interactor_A','IDs_interactor_B')])
# calculate degree
degree_intact <- degree_wo_bidirEdges(g)
degree_intact$proteins <- rownames(degree_intact)
final <- get_bait_prey_total_degree(intact,degree_intact)
write.csv(final,'output/table_bait_prey_total_degree_AP-MS-Y2H_combined.csv',row.names = F)

# test PL property of degree_bait, degree_prey and total_degree
p_bait <- check_powerLaw(as.numeric(final$degree_bait[-which(final$degree_bait == 0)]), plot = T, plot_name = 'plots/plot_bait_degree_AP-MS_Y2H_combined', t = 20, 'Degree')
p_prey <- check_powerLaw(as.numeric(final$degree_prey[-which(final$degree_prey == 0)]), plot = T, plot_name = 'plots/plot_prey_degree_AP-MS_Y2H_combined', t = 20, 'Degree')
p_total <- check_powerLaw(as.numeric(final$total_degree), plot = T, plot_name = 'plots/plot_total_degree_AP-MS_Y2H_combined', t = 20, 'Degree')

pvalue <- c(p_bait$p,p_prey$p,p_total$p)
xmin <- c(xmin_estimated(as.numeric(final$degree_bait[-which(final$degree_bait == 0)])),xmin_estimated(as.numeric(final$degree_prey[-which(final$degree_prey == 0)])),
          xmin_estimated(as.numeric(final$total_degree)))
alpha <- c(alpha_estimated(as.numeric(final$degree_bait[-which(final$degree_bait == 0)])),alpha_estimated(as.numeric(final$degree_prey[-which(final$degree_prey == 0)])),
           alpha_estimated(as.numeric(final$total_degree)))

PL_table <- data.frame(rep('AP-MS_Y2H_combined',3),c('bait_degree','prey_degree','total_degree'),pvalue,xmin,alpha)
colnames(PL_table)[c(1,2)] <- c('method','type')
write.csv(PL_table,'output/table_PL_pvalue_AP-MS_Y2H_combined.csv',row.names = F)

#-------------------------------------------------------------------------------------
# test degree distribution from simulation (Figure 4)
#--------------------------------------------------------------------------------

f <- read.table('databases/Deg_math-simulation.txt')
bs <- check_powerLaw(f$V1,plot=T,plot_name = paste0('plots/plot_Deg'), t = 20, 'Degree')
bs$p
xmin <- xmin_estimated(f$V1)
alpha <- alpha_estimated(f$V1)

# d <- readLines('ergebnis.txt')
# length(d)
# 
# names <- c()
# pvalue <- c()
# xmin <- c()
# alpha <- c()
# n_NA <- c()
# 
# for(i in c(9,13,17,21,25)){
#   print(i)
#   if(length(grep('-',d[i])) != 0){
#     print(i)
#     f <- gsub('\\]','',d[i+1])
#     f <- gsub('\\[','',f)
#     f <- as.numeric(unlist(strsplit(f,', ')))
#     if(length(which(f==0)) != 0){
#       bs <- check_powerLaw(f[-which(f==0)],plot=T,plot_name = paste0('plots/plot_',d[i-1]), t = 20, 'Degree')
#       xmin <- c(xmin,xmin_estimated(f[-which(f==0)]))
#       alpha <- c(alpha,alpha_estimated(f[-which(f==0)]))
#     }else{
#       bs <- check_powerLaw(f,plot=T,plot_name = paste0('plots/plot_',d[i-1]), t = 20, 'Degree')
#       xmin <- c(xmin,xmin_estimated(f))
#       alpha <- c(alpha,alpha_estimated(f))
#     }
#     
#     names <- c(names,d[i-1])
#     pvalue <- c(pvalue,bs$p)
#     n_NA <- c(n_NA,length(which(is.na(bs$bootstraps$xmin) == T)))
#   }
# }
# 
# table <- data.frame(names,pvalue,xmin,alpha,n_NA)
# table$n_NA <- NULL
# write.csv(table,'output/table_PL_pvalue_ergebnis.csv',row.names = F)

#------------------------------------------------------------------------
# correlation between degree and protein mass (Supplementary Figure 7A)
#------------------------------------------------------------------------
hippie_intact <- read.csv('databases/HIPPIE_union_Intact2022_afterReviewed_mapping.csv')
method <- 'Mass_spec'
#method <- 'Y2H'

if(method == 'Mass_spec'){
  index <- grep('coimmunoprecipitation|pull down|tandem affinity purification|mass spectrometry|copurification|affinity chromatography technology',hippie_intact$Interaction_detection_methods)
}else{
  index <- grep('two hybrid',hippie_intact$Interaction_detection_methods)
}

hippie_intact <- hippie_intact[index,]
#degree
g <- unique(hippie_intact[,c('IDs_interactor_A','IDs_interactor_B')])
degree <- degree_wo_bidirEdges(g)
degree$proteins <- rownames(degree)

uniprot_table <- read.csv('databases/uniprot_idmapping_mass.tsv', sep = '\t')
mass_table <- uniprot_table[,c('From','Mass')]

degree$mass <- mass_table$Mass[match(degree$proteins,mass_table$From)]
# remove proteins with NA mass
degree <- degree[-which(is.na(degree$mass) == T),]

library("ggpubr")
ggscatter(degree, x = "mass", y = "degree", 
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "spearman",
          xlab = "Mass(Da)", ylab = "Degree",cor.coef.size = 5) + scale_x_continuous(trans='log2')+
  theme(text = element_text(size = 10), axis.text = element_text(size=10), axis.title = element_text(size=12)) 
ggsave(paste0('plots/scatterplot_spearman_correlation_mass-degree_',method,'.pdf'),height = 8, width = 8,units = 'cm')

#----------------------------------------------------------------------------------
# correlation between degree and protein abundance (Supplementary Figure 7B)
#----------------------------------------------------------------------------------
gtex_protein_relative_abun <- read_excel('databases/Table_S1_gene_info_at_protein_level.xlsx', sheet = 8)
colnames(gtex_protein_relative_abun) <- gtex_protein_relative_abun[1,]
gtex_protein_relative_abun <- gtex_protein_relative_abun[-1,] 

gtex_protein_relative_abun[gtex_protein_relative_abun == "NA"] <- NA

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
# remove genes with more than 50% of NAs
table_median_protein <- table_median_protein[-which(table_median_protein$n_NA > 0.5),]
dim(table_median_protein)

uniprot <- as.data.frame(mapIds(org.Hs.eg.db, keys = table_median_protein$Gene_ID, keytype = "ENSEMBL", column= "UNIPROT"))
colnames(uniprot)[1] <- 'uniprot'
uniprot$ensembl <- rownames(uniprot)

table_median_protein$uniprot <- uniprot$uniprot[match(table_median_protein$Gene_ID,uniprot$ensembl)]

# select specific studies
hippie_intact <- read.csv('databases/HIPPIE_union_Intact2022_afterReviewed_mapping.csv')
method <- 'Mass_spec'
#method <- 'Y2H'

if(method == 'Mass_spec'){
  index <- grep('coimmunoprecipitation|pull down|tandem affinity purification|mass spectrometry|copurification|affinity chromatography technology',hippie_intact$Interaction_detection_methods)
}else{
  index <- grep('two hybrid',hippie_intact$Interaction_detection_methods)
}

hippie_intact <- hippie_intact[index,]
#degree
g <- unique(hippie_intact[,c('IDs_interactor_A','IDs_interactor_B')])
degree <- degree_wo_bidirEdges(g)
degree$proteins <- rownames(degree)

degree$Median <- table_median_protein$Median[match(degree$proteins,table_median_protein$uniprot)]
# remove proteins with median = NA
degree <- degree[-which(is.na(degree$Median) == T),]

ggscatter(degree, x = "Median", y = "degree", 
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "spearman",
          xlab = "Protein abundance median", ylab = "Degree", cor.coef.size = 5) +
  theme(text = element_text(size = 10), axis.text = element_text(size=10), axis.title = element_text(size=12))
ggsave(paste0('plots/scatterplot_spearman_correlation_proteomics-degree_',method,'.pdf'),height = 8, width = 8,units = 'cm')

