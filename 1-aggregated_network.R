source('functions.R')
library(data.table)

#------------------
# read HIPPIE
#------------------

hippie <- read.table('databases/HIPPIE_14Feb2019.txt')
hippie <- add_pubmed_method_hippie(hippie)
uniprot_reviewed <- read.csv('databases/uniprot-compressed_true_download_true_fields_accession_2Creviewed_2C-2022.12.13-09.55.05.21.tsv', sep = '\t')
dim(uniprot_reviewed)
hippie <- filter_hippie(hippie,uniprot_reviewed)

#----------------
# read IntAct
#----------------

intact <- fread('databases/intact.txt', quote = '')
intact <- intact_parsing(intact)
intact <- filter_intact(intact,uniprot_reviewed)
length(unique(intact$Publication_Identifiers))
write.csv(intact,'databases/IntAct_afterFiltering.csv',row.names = F)

intact_bait_prey <- intact[which((intact$Experimental_roles_interactor_A == 'bait' & intact$Experimental_roles_interactor_B == 'prey') | (intact$Experimental_roles_interactor_A == 'prey' & intact$Experimental_roles_interactor_B == 'bait')),]
# frequency of interactions with bait and prey annotation
nrow(unique(intact_bait_prey))/nrow(unique(intact))

#------------------------
# merge HIPPIE and Intact
#------------------------

hippie <- hippie[,c('IDs_interactor_A','IDs_interactor_B','Publication_Identifiers', 'Interaction_detection_methods')]
hippie <- unique(hippie)
hippie$Interaction_detection_methods <- tolower(hippie$Interaction_detection_methods)
hippie$Interaction_detection_methods <- gsub('-',' ',hippie$Interaction_detection_methods)

intact <- intact[,c('IDs_interactor_A','IDs_interactor_B','Publication_Identifiers','Interaction_detection_methods')]
hippie_intact <- as.data.frame(rbind(hippie,intact))
print(paste0('Number of nodes (aggregated network): ',length(union(hippie_intact$IDs_interactor_A,hippie_intact$IDs_interactor_B)))) #17865
#retrieve number of edges
hippie_intact_unique <- unique(hippie_intact[,c('IDs_interactor_A','IDs_interactor_B')])
d_index <- index_edges_duplicates(hippie_intact_unique)
if(length(d_index) > 0){
  hippie_intact_unique <- hippie_intact_unique[-d_index,]
}
print(paste0('Number of edges (aggregated network): ',nrow(hippie_intact_unique)))# 471693
write.csv(hippie_intact, file = 'databases/HIPPIE_union_Intact2022_afterReviewed_mapping.csv', row.names = F)

#-------------------------------------------
# degree distribution of aggregated network
#-------------------------------------------
if(!dir.exists('plots')){
  dir.create('plots')
}

print('Testing power-law hypotesis for degree distribution ...')
g <- unique(hippie_intact[,c('IDs_interactor_A','IDs_interactor_B')])
degree_hippie_intact <- degree_wo_bidirEdges(g)
degree_hippie_intact$proteins <- rownames(degree_hippie_intact)
write.csv(degree_hippie_intact, file = 'output/degree_HIPPIEunionIntact2022.csv',row.names = F)

p <- check_powerLaw(degree_hippie_intact$degree, plot = T, plot_name = 'plots/plot_degree_hippieIntact2022', t = 20, 'Degree')

print(paste0('the p-value of degree distribution is ',p$p)) #0.35
print(paste0('xmin: ',xmin_estimated(degree_hippie_intact$degree)))
print(paste0('alpha is: ',alpha_estimated(degree_hippie_intact$degree)))
print(paste0('ntail is: ',length(which(degree_hippie_intact$degree >= xmin_estimated(degree_hippie_intact$degree)))))

#--------------------------------------
# bait usage distribution (Figure 3B)
#--------------------------------------
print('Testing power-law hypotesis for bait usage distribution ...')

intact <- read.csv('databases/IntAct_afterFiltering.csv')
baits <- get_bait_prey(intact,'bait')
bait_usage <- calculate_bait_usage(intact,baits,'bait')
#check the distribution of the bait usage
p_bait <- check_powerLaw(as.numeric(bait_usage$bait_usage),plot = T,'plots/plot_bait_usage_Intact2022', t = 20, 'Bait usage')

print(paste0('the p-value of bait usage distribution is ',p_bait$p))
print(paste0('xmin: ',xmin_estimated(as.numeric(bait_usage$bait_usage))))
print(paste0('alpha is: ',alpha_estimated(as.numeric(bait_usage$bait_usage))))
print(paste0('ntail is: ',length(which(as.numeric(bait_usage$bait_usage) >= xmin_estimated(as.numeric(bait_usage$bait_usage))))))

# add proteins with bait_usage = 0 and prey_usage >0
preys <- get_bait_prey(intact,'prey')
preys_noBaits <- setdiff(preys,baits)
#retrieve symbol of the preys_noBaits
preys_noBaits_symbol <- as.data.frame(mapIds(org.Hs.eg.db, keys = preys_noBaits, keytype = "UNIPROT", column= "SYMBOL"))
colnames(preys_noBaits_symbol)[1] <- 'symbol'
preys_noBaits_symbol <- preys_noBaits_symbol$symbol[match(preys_noBaits,rownames(preys_noBaits_symbol))]
preys_noBaits_table <- data.frame(preys_noBaits,preys_noBaits_symbol,rep(0,length(preys_noBaits)),rep(0,length(preys_noBaits)))
colnames(preys_noBaits_table) <- colnames(bait_usage)
# add the preys_noBaits to the bait_usage table
bait_usage <- rbind(bait_usage,preys_noBaits_table)
# add degree in the aggregated network (hippie+intact)
bait_usage$full_degree <- degree_hippie_intact$degree[match(bait_usage$bait_uniprot,degree_hippie_intact$proteins)]
write.csv(bait_usage, file = 'output/bait_usage_intact2022.csv',row.names = F)

#-----------------------------------------------------------------------
# check number of interactions for each pubmedID
#-----------------------------------------------------------------------

pubmed <- unique(hippie_intact$Publication_Identifiers)
pubmed <- pubmed[-which(is.na(pubmed) == T)]
table <- get_numInter_studies(hippie_intact,pubmed)
write.csv(table, file = 'output/table_singleStudy_numInter_HIPPIEunionIntact2022.csv',row.names = F)

#----------------------------------------------------
# calculate the degree distribution of single study
#----------------------------------------------------
table <- read.csv('output/table_singleStudy_numInter_HIPPIEunionIntact2022.csv', header = T)
n <- 2
studies <- table$pubmed[which(table$num_inter >= n)]
hippie_intact <- read.csv('databases/HIPPIE_union_Intact2022_afterReviewed_mapping.csv')
calculate_degree_singleStudy(hippie_intact,studies,table,label = 'HIPPIEunionIntact2022',dir = 'output',n,10,20)

#---------------------------------------------------------------------
# calculate the number of non-power and power law studies (Figure 2B)
#---------------------------------------------------------------------

table <- read.csv('output/degree_distr_singleStudy_HIPPIEunionIntact2022_ninter_2_noNA_10.csv')
table_numInter <- read.csv('output/table_singleStudy_numInter_HIPPIEunionIntact2022.csv')
#number of interactions
num <- c(8,20,seq(50,500,50))
table_ratio <- get_nPower_nonPower_table(table,table_numInter,num)
write.csv(table_ratio, file = 'output/table_ratio_HIPPIEunionIntact2022.csv',row.names = F)

library(ggplot2)
table_ratio$num_inter <- as.numeric(gsub('>=','',table_ratio$num_inter))
ggplot(data=table_ratio, aes(x=num_inter, y=ratio)) +
  geom_point() + xlab("Number of PPIs") + ylab('Non-PL/PL') +
  scale_x_continuous(labels = table_ratio$num_inter,breaks = table_ratio$num_inter) +
  theme(axis.text = element_text(size = 10), axis.title = element_text(size = 10)) + theme_bw()
ggsave(paste0('plots/plot_numInter_ratio_HIPPIEunionIntact2022.pdf'),height = 8, width = 13,units = 'cm')

#------------------------------------------------------------------------
# histogram of PPIs number of studies
#------------------------------------------------------------------------

table_study <- read.csv('output/table_singleStudy_numInter_HIPPIEunionIntact2022.csv')

NotFancy <- function(l) {
  l <- format(l, scientific = F)
  parse(text=l)
}

ggplot(table_study, aes(x=num_inter)) + geom_histogram(binwidth = 1) + scale_x_continuous(trans='log10', labels = NotFancy) + xlab('Number of PPIs') + ylab('Number of studies') +
  theme(axis.text = element_text(size = 1), axis.title = element_text(size = 1)) + theme_bw()
ggsave('plots/plot_ppis_numStudies.pdf', height = 10, width = 10, units = 'cm')

#---------------------------------------------------------
### correlation between degree and bait usage
#---------------------------------------------------------

bait_usage <- read.csv('output/bait_usage_intact2022.csv')
#remove rows with bait_usage = 0
bait_usage <- bait_usage[-which(bait_usage$bait_usage == 0),]

print('Pearson correlation between degree and bait usage')
cor.test(bait_usage$bait_usage,bait_usage$full_degree, method = 'pearson') # 0.57

#----------------------------------------------------------------------
# distribution of number of publications
# using data from Pubmed ( https://ftp.ncbi.nih.gov/gene/DATA/)
#----------------------------------------------------------------------
gene2pubmed <- read.csv('databases/gene2pubmed', sep = '\t')
gene2pubmed_human <- gene2pubmed[which(gene2pubmed$X.tax_id == '9606'),]
# count
t <- as.data.frame(table(gene2pubmed_human$GeneID))
p_t <- check_powerLaw(t$Freq, plot = T,'plots/plot_num_studies_gene2pubmed', t = 20,'Number of publications')
print(paste0('the p-value of publication number distribution is ',p_t$p))
print(paste0('xmin: ',xmin_estimated(t$Freq)))
print(paste0('alpha: ',alpha_estimated(t$Freq)))
