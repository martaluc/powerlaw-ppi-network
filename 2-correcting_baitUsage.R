source('functions.R')

#-------------------------------------------------------------------------------------------------------------------------------
# ATTENTION! Before running this script, please do the following:
# 1) Check the number of CPUs (we used 10), set in the following functions: 'calculate_degree_singleStudy' and 'calculate_new_degree_table'
# 2) Eventually change the number of CPUs according to the available resources on your machine.
#-------------------------------------------------------------------------------------------------------------------------------

intact <- read.csv('databases/IntAct_afterFiltering.csv')
pubmed <- unique(intact$Publication_Identifiers)
length(pubmed)
table_intact <- as.data.frame(get_numInter_studies(intact,pubmed))
table_intact <- add_methods(table_intact,intact)
write.csv(table_intact, file = paste0('output/','table_singleStudy_numInter_Intact.csv'),row.names = F)
ninter <- 2
table_intact <- read.csv('output/table_singleStudy_numInter_Intact.csv')
studies_intact <- table_intact$pubmed[which(table_intact$num_inter >= ninter)]
calculate_degree_singleStudy(intact,studies_intact,table_intact,label = 'Intact',dir = 'output',ninter,nRemove=10,tr=10)

#--------------------------------------------------------
# add the number of baits and preys in the previous table
#---------------------------------------------------------
final_table <- read.csv(paste0('output/degree_distr_singleStudy_Intact_ninter_',ninter,'_noNA_10.csv'))
final_table <- add_number_baits_preys(intact,final_table)
write.csv(final_table, file = paste0('output/degree_distr_singleStudy_Intact_ninter_',ninter,'_noNA_10.csv'),row.names = F)

#------------------------------------------------------------------
# degree recalculation after bias correction (asymmetric design) 
#------------------------------------------------------------------

ninter <- 2
table_study <- read.csv(paste0('output/degree_distr_singleStudy_Intact_ninter_',ninter,'_noNA_10.csv'))
#select only power law studies
table_power <- table_study[which(table_study$pvalue >= 0.1),]
# eliminate studies with 0 baits and 0 preys
table_power <- table_power[-which(table_power$n_baits == 0 & table_power$n_preys == 0),]
# eliminate studies with 1 bait and 1 prey
table_power <- table_power[-which(table_power$n_baits == 1 & table_power$n_preys == 1),]
calculate_new_degree_table(intact,table_power,ninter,nRemove=10,tr=10,'output')


#---------------------------------------------------------------
# wilcoxon test for testing size balance difference (Figure 4C)
#---------------------------------------------------------------
ninter <- 2
final_table_noNA <- read.csv(paste0('output/degree_bis_table_ninter',ninter,'_noNA_10.csv'))
i <- 200
final_table_noNA <- final_table_noNA[which(final_table_noNA$n_inter >= i),]
wilcox_test_ratio(final_table_noNA,'boxplot_ratio_bait_prey_or_prey_bait_noNA_10',i,'plots')

#---------------------------------------------------------------------------
# distribution of the ratio between bait and prey in PL studies (Figure 4B)
#---------------------------------------------------------------------------
ninter <- 2
final_table_noNA <- read.csv(paste0('output/degree_bis_table_ninter',ninter,'_noNA_10.csv'))
final_table_200 <- final_table_noNA[which(final_table_noNA$n_inter >= 200 ),]


ggplot(final_table_200, aes(ratio_bait_prey)) +
  geom_bar() + scale_x_binned() + ylab('Number of studies') + xlab('Size balance') +
  theme(axis.text = element_text(size = 10), axis.title = element_text(size = 10)) +theme_bw()
ggsave('plots/plot_distribution_ratio_PLstudies.pdf',height = 10.3, width = 10, units = 'cm')

