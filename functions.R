library(poweRlaw)
library(linkcomm)
library(ggplot2)
library(ggpubr)
library(plyr)
library(org.Hs.eg.db)
library(tidyr)
library(DOSE)
library(clusterProfiler)
library(ReactomePA)
library(igraph)

##########################################################################################################
# function that adds publication IDs and method for each row in HIPPIE database
# input: hippie database (downloaded by  http://cbdm-01.zdv.uni-mainz.de/~mschaefer/hippie/download.php)
# output: hippie with two more columns ('Publication_Identifiers' and 'Interaction_detection_methods')
###########################################################################################################

add_pubmed_method_hippie <- function(hippie){
  
  # add pubmed column 
  length(grep('pmids',hippie$info))
  pubmed <- strsplit(hippie$info, split = ';')
  pubmed <- unlist(lapply(pubmed,FUN = function(x) {ifelse(length(grep('pmids',x)) != 0,x[grep('pmids',x)],NA)}))
  length(pubmed)
  hippie$Publication_Identifiers <- pubmed
  hippie$Publication_Identifiers <- gsub('pmids:','',hippie$Publication_Identifiers)
  
  # add method column
  method <- strsplit(hippie$info, split = ';')
  method <- unlist(lapply(method,FUN = function(x) {ifelse(length(grep('experiments',x)) != 0,x[grep('experiments',x)],NA)}))
  hippie$Interaction_detection_methods <- method
  hippie$Interaction_detection_methods <- gsub('experiments:','',hippie$Interaction_detection_methods)
  
  hippie_parsing <- mutate(hippie,Publication_Identifiers = strsplit(as.character(Publication_Identifiers), ",")) 
  hippie_parsing <- as.data.frame(unnest(hippie_parsing,Publication_Identifiers))
  length(unique(hippie_parsing$Publication_Identifiers))
  
  hippie_parsing <- mutate(hippie_parsing,Interaction_detection_methods = strsplit(as.character(Interaction_detection_methods), ",")) 
  hippie_parsing <- as.data.frame(unnest(hippie_parsing,Interaction_detection_methods))
  return(hippie_parsing)
}

########################################################################
# function that retrieves only reviewed uniprotID from HIPPIE database
# hippie: dataframe
# uniprot_reviewed: dataframe that includes reviewed uniprotIDs
#########################################################################

filter_hippie <- function(hippie,uniprot_reviewed){
  
  # retrieve the reviewed ID when there are multiple entry_names
  multiple_ID_a <- unique(hippie$UniprotID_A[grep(',',hippie$UniprotID_A)])
  multiple_ID_b <- unique(hippie$UniprotID_B[grep(',',hippie$UniprotID_B)])
  multiple_ID <- union(multiple_ID_a,multiple_ID_b)
  s <- strsplit(multiple_ID,split = ',')
  Entry.Name <- c()
  Entry.Multiple <- c()
  for(i in seq(1,length(s))){
    index <- which(s[[i]] %in% uniprot_reviewed$Entry.Name)
    if(length(index) >= 1){
      Entry.Name <- c(Entry.Name,s[[i]][index[1]])
    }else{
      Entry.Name <- c(Entry.Name,NA)
    }
    Entry.Multiple <- c(Entry.Multiple,paste0(s[[i]], collapse = ','))
  }
  uniprot_reviewed_multiple <- as.data.frame(cbind(Entry.Name,Entry.Multiple))
  uniprot_reviewed_multiple$Entry <- uniprot_reviewed$Entry[match(uniprot_reviewed_multiple$Entry.Name,uniprot_reviewed$Entry.Name)]
  uniprot_reviewed_multiple$Entry.Name <- NULL
  colnames(uniprot_reviewed_multiple)[1] <- 'Entry.Name'
  
  tot_uniprot_reviewed <- rbind(uniprot_reviewed[,c('Entry','Entry.Name')],uniprot_reviewed_multiple[,c('Entry','Entry.Name')])
  
  # mapping
  hippie$IDs_interactor_A <- tot_uniprot_reviewed$Entry[match(hippie$UniprotID_A,tot_uniprot_reviewed$Entry.Name)]
  hippie$IDs_interactor_B <- tot_uniprot_reviewed$Entry[match(hippie$UniprotID_B,tot_uniprot_reviewed$Entry.Name)]
  
  # eliminate NA's
  indexNA1 <- which(is.na(hippie$IDs_interactor_A) == T)
  indexNA2 <- which(is.na(hippie$IDs_interactor_B) == T)
  nonRev_ID <- union(hippie$UniprotID_A[indexNA1],hippie$UniprotID_B[indexNA2])
  hippie <- hippie[-union(indexNA1,indexNA2),]
  return(hippie)
}

#####################################################################################################
# function that parses the IntAct database retrieving only
# human interactions
# d: dataframe (IntAct downloaded by https://ftp.ebi.ac.uk/pub/databases/intact/2022-02-03/psimitab/)
######################################################################################################

intact_parsing <- function(d){
  
  # 1) only human interactions
  d_human <- d[grep('taxid:9606',d$`Taxid interactor A`),]
  d_human <- d_human[grep('taxid:9606',d_human$`Taxid interactor B`),]
  dim(d_human)
  unique(d_human$`Taxid interactor A`)
  unique(d_human$`Taxid interactor B`)
  
  # 2) retrieve interaction only between proteins
  d_human <- d_human[grep('(protein)',d_human$`Type(s) interactor A`),]
  d_human <- d_human[grep('(protein)',d_human$`Type(s) interactor B`),]
  dim(d_human)
  unique(d_human$`Type(s) interactor A`)
  unique(d_human$`Type(s) interactor B`)
  
  # 3) retrieve only the role in the parenthesis
  d_human$`Experimental role(s) interactor A` <- gsub("[\\(\\)]", "", regmatches(d_human$`Experimental role(s) interactor A`, gregexpr("\\(.*?\\)", d_human$`Experimental role(s) interactor A`)))
  d_human$`Experimental role(s) interactor B` <- gsub("[\\(\\)]", "", regmatches(d_human$`Experimental role(s) interactor B`, gregexpr("\\(.*?\\)", d_human$`Experimental role(s) interactor B`)))
  
  # 4) retain only the rows that include 'uniprotkb'
  indexa <- grep('uniprotkb',d_human$`#ID(s) interactor A`)
  indexb <- grep('uniprotkb',d_human$`ID(s) interactor B`)
  d_human <- d_human[intersect(indexa,indexb),]
  dim(d_human)
  grep('intact',d_human$`#ID(s) interactor A`)
  grep('intact',d_human$`ID(s) interactor B`)
  grep('refseq',d_human$`#ID(s) interactor A`)
  length(grep('uniprot',d_human$`#ID(s) interactor A`))
  length(grep('uniprot',d_human$`ID(s) interactor B`))
  
  # 5) eliminate the string 'uniprotkb'
  d_human$`#ID(s) interactor A` <- gsub('uniprotkb:','',d_human$`#ID(s) interactor A`)
  d_human$`ID(s) interactor B` <- gsub('uniprotkb:','',d_human$`ID(s) interactor B`)
  
  # 6) eliminate -1 from the name of the uniprot ID
  d_human$`#ID(s) interactor A` <- gsub('-1','',d_human$`#ID(s) interactor A`)
  d_human$`ID(s) interactor B` <- gsub('-1','',d_human$`ID(s) interactor B`)
  
  # 7) remove '-' and considering them as original protein
  d_human$`#ID(s) interactor A` <- gsub('-.*','',d_human$`#ID(s) interactor A`)
  d_human$`ID(s) interactor B` <- gsub('-.*','',d_human$`ID(s) interactor B`)
  dim(d_human)
  grep('-',d_human$IDs_interactor_A)
  grep('-',d_human$IDs_interactor_B)
  
  # 8) retrieve only pubmedID
  d_human$`Publication Identifier(s)` <- gsub("^.*pubmed:|\\|.*$","",d_human$`Publication Identifier(s)`)
  
  colnames(d_human)[c(1,2)] <- c('IDs_interactor_A','IDs_interactor_B') 
  
  colnames(d_human) <- gsub('\\(','',colnames(d_human))
  colnames(d_human) <- gsub('\\)','',colnames(d_human))
  colnames(d_human) <- gsub(' ','_',colnames(d_human))
  
  return(d_human)
}

###############################################################
# function that retrieves only the reviwed UniprotIDs from IntAct
# intact: dataframe
# uniprot_reviewed: dataframe that includes reviewed uniprotIDs
###############################################################
filter_intact <- function(intact,uniprot_reviewed){
  
  intact <- intact[,c('IDs_interactor_A','IDs_interactor_B','Publication_Identifiers','Interaction_detection_methods', 'Experimental_roles_interactor_A', 'Experimental_roles_interactor_B','Interaction_types')]
  dim(intact)
  
  #retain only proteins with reviewed uniprotID
  index <- which(intact$IDs_interactor_A %in% uniprot_reviewed$Entry & intact$IDs_interactor_B %in% uniprot_reviewed$Entry)
  intact <- intact[index,]
  dim(intact)
  length(which(intact$IDs_interactor_B %in% uniprot_reviewed$Entry))
  length(union(intact$IDs_interactor_A,intact$IDs_interactor_B)) #17220
  method <- gsub('psi-mi:\"MI:','',intact$Interaction_detection_methods)
  method <- gsub('[0-9]','',method)
  method <- gsub('\\"\\(|\\)','',method)
  method <- gsub('"','',method)
  intact$Interaction_detection_methods <- method
  intact$Interaction_detection_methods <- gsub('-',' ',intact$Interaction_detection_methods)
  interaction <- gsub('psi-mi:\"\"MI:','',intact$Interaction_types)
  interaction <- gsub('[0-9]','',interaction)
  interaction <- gsub('\\"\\(|\\)','',interaction)
  interaction <- gsub('"','',interaction)
  intact$Interaction_types <- interaction
  
  return(intact)
}

#########################################################################################
# function that retrieves the indices of duplicated bidirectional edges (ex: A--B; B--A)
# input: network database with two columns
# output: indices of duplicated edges
########################################################################################
index_edges_duplicates <- function(network,verbose = F){
  xx <- cbind(as.character(network[, 1]), as.character(network[,2]))
  edges <- integer.edgelist(network)$edges
  ne <- nrow(edges)
  loops <- rep(0, ne)
  dups <- rep(0, ne)
  out <- .C("edgeDuplicates", as.integer(edges[, 1]), as.integer(edges[,2]), as.integer(ne), loops = as.integer(loops), dups = as.integer(dups), 
            as.logical(verbose))
  dups <- which(out$dups == 1)
  return(dups)
}

#####################################################################
# function to calculate the degree removing the bidirectional edges
# input: network with two columns
# output: dataframe with the degree as column and uniprotIDs as rows
#####################################################################

degree_wo_bidirEdges <- function(g){
  
  g_index <- index_edges_duplicates(g)
  #eliminate bidirectional edges
  if(length(g_index) > 0){
    g <- g[-g_index,]
    #print(paste0('there are ',length(g_index),' bidirectional edges'))
  }
  graph <- graph_from_data_frame(g)
  degree <- as.data.frame(degree(graph,normalized = F))
  colnames(degree) <- 'degree'
  return(degree)
}

#########################################################################################################################
# function to check the power-law distribution
# input:
# - data = degree vector
# - plot_name = name of the plot file
# - t = number of threads
# - xlabel = 'Degree' if the degree distribution is calculated, 'Baiu usage' is the bait usage distribution is calculated
# - myseed = seed used to have same random simulated distribution in bootstrapping
# - n_sims = number of simulations in bootstrapping (bootstrap_p function); 100 is the default suggested by poweRlaw developers
# output: bootstrap_p output which includes p-value, goodness-of-fit, and bootstraps
##########################################################################################################################

check_powerLaw <- function(data, plot = F, plot_name ='', t = 1, xlabel, myseed=1, n_sims=100){
  
  m_pl = displ$new(data)
  est = estimate_xmin(m_pl)
  m_pl$setXmin(est)
  #print(est)
  bs_p = bootstrap_p(m_pl,threads = t, seed = myseed, no_of_sims = n_sims)
  #print(summary(bs_p$bootstraps$xmin))
  p <- bs_p$p
  #plot
  if(plot == T){
    xy <- plot(m_pl)
    xy$y <- signif(xy$y, digits = 3)
    l <- lines(m_pl, col = 'green')
    scaleFUN <- function(x) signif(x,digits = 3)
    ggplot(xy, aes(x = x, y = y)) + geom_point() + geom_line(data = l, color = "red") +
      scale_x_continuous(trans='log2') + scale_y_continuous(trans='log2', labels = scaleFUN) + xlab(xlabel) + ylab('Frequency') +
      theme(axis.text = element_text(size = 10), axis.title = element_text(size = 10)) + theme_bw()
    ggsave(paste0(plot_name,'.pdf'),height = 6, width = 6,units = 'cm')
  }
  
  return(bs_p)
}

###############################################################
# functions that estimates the xmin in the power law fitting
# data: numeric vector (ex: degree vector)
# output: estimated xmin
###############################################################

xmin_estimated <- function(data){
  
  m_pl = displ$new(data)
  est = estimate_xmin(m_pl)
  
  return(est$xmin)
}

###########################################################################
# functions that estimates the alpha (exponent) in the power law fitting
# data: numeric vector (ex: degree vector)
# output: estimated alpha
##########################################################################

alpha_estimated <- function(data){
  
  m_pl = displ$new(data)
  est = estimate_xmin(m_pl)
  
  return(est$pars)
}

###########################################################################
# function that gets baits or preys
# input:
# - intact = intact database (since it has info about experimental role)
# - label = 'bait' or 'prey' if we want to get baits or preys, respectively
# output: vector of baits or preys
############################################################################

get_bait_prey <- function(intact,label){
  a <- intact$IDs_interactor_A[which(intact$Experimental_roles_interactor_A == label)]
  b <- intact$IDs_interactor_B[which(intact$Experimental_roles_interactor_B == label)]
  tot <- unique(union(a,b))
  return(tot)
}

#######################################################################################################
# it calculates the number of times each protein in IntAct is used as bait or prey (bait/prey usage)
# input: 
# - intact = intact dataframe
# - baits = vector of baits or preys
# - label = 'bait' or 'prey'
# output: dataframe including uniprotIDs, symbols, bait/prey usage, degree of proteins 
#         considering only the interactions where the protein has been tested as a bait or prey
########################################################################################################

calculate_bait_usage <- function(intact,baits,label){
  bait_usage <- c()
  degree_bait <- c()
  for(b in baits){
    table_a <- intact[which(intact$IDs_interactor_A == b & intact$Experimental_roles_interactor_A == label),]
    table_b <- intact[which(intact$IDs_interactor_B == b & intact$Experimental_roles_interactor_B == label),]
    table <- rbind(table_a,table_b)
    bait_usage <- c(bait_usage,length(unique(table$Publication_Identifiers)))
    g <- unique(table[,c('IDs_interactor_A','IDs_interactor_B')])
    degree <- degree_wo_bidirEdges(g)
    degree$proteins <- rownames(degree)
    degree_bait <- c(degree_bait,degree$degree[which(degree$proteins == b)])
  }
  bait_symbol <- as.data.frame(mapIds(org.Hs.eg.db, keys = baits, keytype = "UNIPROT", column= "SYMBOL"))
  colnames(bait_symbol)[1] <- 'symbol'
  bait_symbols <- bait_symbol$symbol[match(baits,rownames(bait_symbol))]
  f <- as.data.frame(cbind(baits,bait_symbols,bait_usage,degree_bait))
  if(label == 'bait'){
    colnames(f) <- c('bait_uniprot','bait_symbol','bait_usage','degree_bait')
  }
  if(label == 'prey'){
    colnames(f) <- c('prey_uniprot','prey_symbol','prey_usage','degree_prey')
  }
  return(f)
}

##############################################################################################################
# calculate the number of interactions (PPIs) for each study
# input:
# - network: dataframe
# - pubmed: vector of pubmedIDs
# output: matrix including for each study (pubmedID) the number of interactions and the number of proteins
##############################################################################################################

get_numInter_studies <- function(network,pubmed){
  
  pubmedID <- c()
  num_inter <- c()
  n_proteins <- c()
  
  for(i in pubmed){
    p <- unique(network[which(network$Publication_Identifiers == i),c('IDs_interactor_A','IDs_interactor_B')])
    d_index <- index_edges_duplicates(p)
    #eliminate bidirectional edges with the same pubmedID
    if(length(d_index) > 0){
      p <- p[-d_index,]
      print(i)
    }
    pubmedID <- c(pubmedID,i)
    num_inter <- c(num_inter,nrow(p))
    n_proteins <- c(n_proteins,length(unique(union(p$IDs_interactor_A,p$IDs_interactor_B))))
  }
  
  table <- cbind(pubmedID,num_inter,n_proteins)
  return(table)
}

######################################################################################################
# function that tests if the degree distribution of each study (identified by pubmedID) is power-law
# input:
# - network= dataframe
# - studies = vector of pubmedIDs
# - table = dataframe with pubmedIDs and number of PPIs
# - label = name of output file
# - dir = output directory
# - n = number of minimum PPIs
# - nRemove = threshold used to remove studies with NAs in random distributions: we used nRemove=10 (=10% of the 100 bootstrapping simulations)
# - tr = number of threads
# output: dataframe including for each study (pubmedID), p-value, estimated xmin, 
# number of values > xmin (ntail), number of proteins, number of unique degree values,
# goodness-of-fit (gof), and estimated alpha. The output is saved in a csv file.
######################################################################################################

calculate_degree_singleStudy <- function(network,studies,table,label,dir,n,nRemove,tr){
  
  pvalue <- c()
  pubmedID <- c()
  xmin_est <- c()
  xmax <- c()
  unique_degree <- c()
  ntail <- c()
  len_degree <- c()
  gof <- c()
  alpha <- c()
  studies_na_10 <- c()
  
  # for each study it tests the power-law property
  for(i in studies){
    print(i)
    g <- unique(network[which(network$Publication_Identifiers == i),c('IDs_interactor_A','IDs_interactor_B')])
    degree <- degree_wo_bidirEdges(g)
    if(is.na(xmin_estimated(degree$degree))){
      print(paste0('Unable to estimate xmin for ',i,' study \n'))
      next
    }
    p <- check_powerLaw(degree$degree, plot = F, plot_name = paste0(dir,'/plot_eachStudy/plot_pubmedID ',i), t = tr, 'Degree',myseed = 1)
    
    # studies_na_10 includes studies with more than 10 NAs
    s <- summary(degree$degree)
    if(length(which(is.na(p$bootstraps$xmin) == T)) >= nRemove){
      studies_na_10 <- c(studies_na_10,i)
    }
    xmin_est <- c(xmin_est,xmin_estimated(degree$degree))
    alpha <- c(alpha,alpha_estimated(degree$degree))
    xmax <- c(xmax,s[6])
    unique_degree <- c(unique_degree,length(unique(degree$degree)))
    ntail <- c(ntail,length(which(degree >= xmin_estimated(degree$degree))))
    len_degree <- c(len_degree,length(degree$degree))
    
    pvalue <- c(pvalue,p$p)
    gof <- c(gof,p$gof)
    pubmedID <- c(pubmedID,i)
    
  }
  
  final_table <- as.data.frame(cbind(pubmedID,pvalue,xmin_est,xmax,ntail,len_degree,unique_degree,gof,alpha))
  final_table$n_inter <- table$num_inter[match(final_table$pubmedID,table$pubmedID)]
  final_table$methods <-  table$methods[match(final_table$pubmedID,table$pubmedID)]
  
  if(length(studies_na_10) != 0){
    final_table_noNA_10 <- final_table[-which(final_table$pubmedID %in% studies_na_10),]
    write.csv(final_table_noNA_10, file = paste0(dir,'/degree_distr_singleStudy_',label,'_ninter_',n,'_noNA_',nRemove,'.csv'), row.names = F)
    write.table(studies_na_10, file = paste0(dir,'/studies_',nRemove,'NAbootstrapping_',label,'.txt'), quote = F, col.names = F, row.names = F)
  }
  
}

##############################################################################################################################
# function that calculates the number of PL and NPL studies and the ratio between NPL and PL studies
# input:
# - table= dataframe that includes the p-value (plausibility of the PL hypothesis) for single punmedID
# - table_numInter = dataframe of pubmedIDs and PPIs
# - num = number of PPIs threshold
# output: dataframe which includes the number of NPL, PL studies with more than a certain number of PPIs specified in 'num'
##############################################################################################################################

get_nPower_nonPower_table <- function(table,table_numInter,num){
  
  n_noPowerLaw_study <- c()
  n_powerLaw_study <- c()
  n_study <- c()
  n_powerLaw_ntail50_study <- c()
  n_powerLaw_ntail50_alpha_study <- c()
  
  for(i in num){
    subtable <- table[which(table$n_inter >= i),]
    n_noPowerLaw_study <- c(n_noPowerLaw_study, length(which(subtable$pvalue < 0.1)))
    n_powerLaw_study <- c(n_powerLaw_study, length(which(subtable$pvalue >= 0.1)))
    n_study <- c(n_study,length(table_numInter$pubmedID[which(table_numInter$num_inter >= i)]))
    n_powerLaw_ntail50_study <- c(n_powerLaw_ntail50_study, length(which(subtable$pvalue >= 0.1 & subtable$ntail >= 50)))
    n_powerLaw_ntail50_alpha_study <- c(n_powerLaw_ntail50_alpha_study,length(which(subtable$pvalue >= 0.1 & subtable$ntail >= 50 & subtable$alpha >2 & subtable$alpha<3)))
    
  }
  
  num <- paste('>=',num)
  final_table <- as.data.frame(cbind(num,n_study,n_noPowerLaw_study,n_powerLaw_study,n_powerLaw_ntail50_study,n_powerLaw_ntail50_alpha_study))
  colnames(final_table)[1] <- 'num_inter'
  final_table$ratio <- as.numeric(final_table$n_noPowerLaw_study)/as.numeric(final_table$n_powerLaw_study)
  return(final_table)
}

#########################################################################
# function that adds methods to each pubmedID
# it counts how many interactions are detected for each method
# and for each study it assigns the method with the highest frequency
# input:
# - final_table: table that includes the number of PPIs for each study
# - network: data frame of the network
# output: final_table with an extra column (methods)
#########################################################################

add_methods <- function(final_table,network, type = 'method'){
  
  methods <- c()
  for(j in final_table$pubmedID){
    if(j == 9e+07){
      j <- 	'90000000'
    }
    if(type == 'method'){
      method <- network$Interaction_detection_methods[which(network$Publication_Identifiers == j)]
    }
    if(type == 'interaction_type'){
      method <- network$Interaction_types[which(network$Publication_Identifiers == j)]
    }
    if(length(unique(method)) == length(unique(method[which(is.na(method) == T)]))){
      methods <- c(methods,NA)
    }else{
      counter <- as.data.frame(table(method))
      m <- as.character(counter$method[which(counter$Freq == max(counter$Freq))])
      #print(max(counter$Freq))
      if(max(counter$Freq) == -Inf){
        print(j)
      }
      if(length(m) > 1){
        #m <- m[sample(1:length(m),1)]
        m <- paste(m,collapse = ',')
      }
      methods <- c(methods,m)
    }
  }
  if(type == 'method'){
    final_table$methods <- methods
  }
  if(type == 'interaction_type'){
    final_table$interaction_type <- methods
  }
  
  return(final_table)
}

######################################################################################
# add number of baits and prey
# input: 
# - intact: data frame of intact database
# - final_table: table including the p-value of each study
# output: final_table with three extra columns (number of baits, number of preys 
# and number of common baits and preys)
######################################################################################

add_number_baits_preys <- function(intact,final_table){
  n_baits <- c()
  n_preys <- c()
  common_bait_prey <- c()
  for(k in final_table$pubmedID){
    oneStudy <- intact[which(intact$Publication_Identifiers == k),]
    bait <- get_bait_prey(oneStudy,'bait')
    n_baits <- c(n_baits,length(bait))
    prey <- get_bait_prey(oneStudy,'prey')
    n_preys <- c(n_preys,length(prey))
    common_bait_prey <- c(common_bait_prey,length(intersect(prey,bait)))
  }
  final_table$n_baits <- n_baits 
  final_table$n_preys <- n_preys
  final_table$n_common_bait_prey <- common_bait_prey
  return(final_table)
}

#################################################
# check if there are more baits than preys
# s = pubmedID
# output: 'bait' if the number of baits <= prey or 
# 'prey' if the number of preys < baits
#################################################

check_numberBaits_preys <- function(table_input,s){
  nbaits <- table_input$n_baits[which(table_input$pubmedID == s)]
  npreys <- table_input$n_preys[which(table_input$pubmedID == s)]
  
  if(nbaits <= npreys){
    label <- 'bait'
  }else{
    label <- 'prey'
  }
  return(label)
}

#####################################################################################################
# function that calculates the degree without considering interactions where a protein has been
# tested as a bait or prey (see Methods of the manuscript)
# input:
# - human_data_study: network for a specific study
# - degree: vector of degree
# - label: 'bait' or 'prey'
# output: vector of the degree recalculation after correcting for bait usage
####################################################################################################

calculate_degree_bis <- function(human_data_study,degree,label){
  
  degree_bis <- c()
  for (p in degree$proteins) {
    
    t <- human_data_study[which(human_data_study$IDs_interactor_A == p | human_data_study$IDs_interactor_B == p),]
    index_a <- which(t$IDs_interactor_A == p & t$Experimental_roles_interactor_A == label)
    index_b <- which(t$IDs_interactor_B == p & t$Experimental_roles_interactor_B == label)
    index <- union(index_a,index_b)
    if(length(index)!= 0){
      t <- t[-index,]
    }
    if(nrow(t) != 0){
      g <- unique(t[,c('IDs_interactor_A','IDs_interactor_B')])
      g_index <- index_edges_duplicates(g)
      if(length(g_index) > 0){
        g <- g[-g_index,]
        print(paste0('there are ',length(g_index),' bidirectional edges'))
      }
      graph <- graph_from_data_frame(g)
      degree_bis <- c(degree_bis,degree(graph,v = p,normalized = F))
    }
    else{
      degree_bis <- c(degree_bis,0)
    }
  }
  return(degree_bis)
}


################################################################################################################
# add some information to the final table that includes the degree recalculation
# input:
# - human_data: dataframe of the network
# - final_table: dataframe that includes the p-value of the degree for each powerlaw study after the correction
# - table_input: dataframe that includes only power-law studies with the associated p-value
# output: final_table that includes (besides its original columns), the p-value before degree correction,
# number of baits, number of preys, number of interactions
#################################################################################################################

add_info_degreeBis_table <- function(human_data,final_table,table_input){
  
  final_table$pvalue_original <- table_input$pvalue[match(final_table$pubmedID,table_input$pubmedID)]
  # add the methods
  final_table$methods <- table_input$methods[match(final_table$pubmedID,table_input$pubmedID)]
  final_table$n_baits <- table_input$n_baits[match(final_table$pubmedID,table_input$pubmedID)]
  final_table$n_preys <- table_input$n_preys[match(final_table$pubmedID,table_input$pubmedID)]
  final_table$n_inter <- table_input$n_inter[match(final_table$pubmedID,table_input$pubmedID)]
  final_table$gof_original <- table_input$gof[match(final_table$pubmedID,table_input$pubmedID)]
  final_table$xmin_est_original <- table_input$xmin_est[match(final_table$pubmedID,table_input$pubmedID)]
  final_table$n_common_bait_prey <- table_input$n_common_bait_prey[match(final_table$pubmedID,table_input$pubmedID)]
  
  final_table <- final_table[,c('pubmedID','pvalue_original','pvalue','gof_original','gof','xmin_est_original','xmin','methods','n_inter','n_baits','n_preys','n_common_bait_prey','n_changeDegree_proteins')]
  final_table$ratio_bait_prey <- ifelse(final_table$n_baits < final_table$n_preys,final_table$n_baits/final_table$n_preys,final_table$n_preys/final_table$n_baits)
  
  return(final_table)
}
###############################################################################################
# function that calculates the degree after correction (see Methods in the paper)
# input:
# - human_data: dataframe of the network (IntAct)
# - table_input: dataframe that includes only power-law studies with the associated p-values
# - ninter: 2
# - nRemove: threshold used to remove studies with NAs in random distributions: we used nRemove=10 (=10% of the 100 bootstrapping simulations)
# - tr: number of threads
# - dir_out: directory where the output is saved
###############################################################################################

calculate_new_degree_table <- function(human_data,table_input,ninter,nRemove,tr,dir_out){
  pvalue <- c()
  pubmedID <- c()
  gof <- c()
  xmin <- c()
  studies_na_10 <- c()
  n_changeDegree_proteins <- c()

  
  studies <- table_input$pubmedID
  for(s in studies){
    print(s)
    human_data_study <- human_data[which(human_data$Publication_Identifiers == s),c('IDs_interactor_A','IDs_interactor_B','Publication_Identifiers','Experimental_roles_interactor_A','Experimental_roles_interactor_B','Interaction_detection_methods')]
    # calculate degree
    g <- unique(human_data_study[,c('IDs_interactor_A','IDs_interactor_B')])
    degree <- degree_wo_bidirEdges(g)
    degree$proteins <- rownames(degree)
    
    label <- check_numberBaits_preys(table_input,s)
    degree_bis <- calculate_degree_bis(human_data_study,degree,label)
    
    degree$degree_bis <- degree_bis
    #write.csv(degree, file = paste0('degree_bis/degree_bis_',s,'.csv'),row.names = F)
    
    if(length(which(degree$degree_bis == 0)) != 0){
      if(is.na(xmin_estimated(degree$degree_bis[-which(degree$degree_bis == 0)]))){
        print(paste0('Unable to estimate xmin for ',s,' study \n'))
        next
      }
      pvalue_s <- check_powerLaw(degree$degree_bis[-which(degree$degree_bis == 0)], plot = F, plot_name = paste0('plot_degree_bis/plot_',s),t=tr,xlabel = 'Degree',myseed = as.numeric(s))
      pvalue <- c(pvalue,pvalue_s$p)
      gof <- c(gof,pvalue_s$gof)
      xmin <- c(xmin,xmin_estimated(degree$degree_bis[-which(degree$degree_bis == 0)]))

      if(length(which(is.na(pvalue_s$bootstraps$xmin) == T)) >= nRemove)
        studies_na_10 <- c(studies_na_10,s)
    }
    else{
      if(is.na(xmin_estimated(degree$degree_bis))){
        print(paste0('Unable to estimate xmin for ',s,' study \n'))
        next
      }
      pvalue_s <- check_powerLaw(degree$degree_bis, plot = F, plot_name = paste0('plot_degree_bis/plot_',s),t=tr,xlabel = 'Degree', myseed = as.numeric(s))
      pvalue <- c(pvalue,pvalue_s$p)
      gof <- c(gof,pvalue_s$gof)
      xmin <- c(xmin,xmin_estimated(degree$degree_bis))
      if(length(which(is.na(pvalue_s$bootstraps$xmin) == T)) >= nRemove)
        studies_na_10 <- c(studies_na_10,s)
    }
    pubmedID <- c(pubmedID,s)
    n_changeDegree_proteins <- c(n_changeDegree_proteins,length(which(degree$degree_bis < degree$degree)))
  }
  
  final_table <- as.data.frame(cbind(pubmedID,pvalue,gof,xmin,n_changeDegree_proteins))
  final_table <- add_info_degreeBis_table(human_data,final_table,table_input)
  #write.csv(final_table, file = paste0('degree_bis_table_ninter',ninter,'.csv'),row.names = F)
  if(length(studies_na_10) > 0){
    final_table_noNA_10 <- final_table[-which(final_table$pubmedID %in% studies_na_10),]
    write.csv(final_table_noNA_10, file = paste0(dir_out,'/degree_bis_table_ninter',ninter,'_noNA_',nRemove,'.csv'),row.names = F)
  }
  #return(final_table)
}

###################################################################################################
# function that uses wilcoxon to test the difference of size balance between studies that switch
# from power-law to non-power law and those that still remain power-law
# input:
# - final_table = dataframe that includes p-value for each study
# - plot_name = name of the output plot
# - ninter= threshold of the number of interactions
# - dir_out= directory where the plot is saved
####################################################################################################

wilcox_test_ratio <- function(final_table,plot_name,ninter,dir_out){
  
  data_switch <- cbind(final_table$ratio_bait_prey[which(final_table$pvalue < 0.1)],rep('PL-nonPL',length(final_table$ratio_bait_prey[which(final_table$pvalue < 0.1)])))
  data_noswitch <- cbind(final_table$ratio_bait_prey[which(final_table$pvalue >= 0.1)],rep('PL-PL',length(final_table$ratio_bait_prey[which(final_table$pvalue >= 0.1)])))
  print(paste0('number of studies powerlaw-nonPowerlaw ',nrow(data_switch)))
  print(paste0('number of studies powerlaw-powerlaw ',nrow(data_noswitch)))
  
  data <- as.data.frame(rbind(data_switch,data_noswitch))  
  colnames(data) <- c('ratio_bait_prey','type')
  data$ratio_bait_prey <- as.numeric(data$ratio_bait_prey)
  data$type <- as.factor(data$type)
  
  library(ggpubr)
  my_comparison <- list(c('PL-nonPL','PL-PL'))
  ggplot(data, aes(x=type, y=ratio_bait_prey)) + 
    geom_boxplot() + theme(axis.title.x = element_blank())+ stat_compare_means(comparisons = my_comparison, method = 'wilcox.test', method.args = list(alternative = 'less')) +
    theme(axis.text=element_text(size=10), axis.title = element_text(size = 10)) + ylab ('Size balance') + theme_bw()
  #ggsave(paste0(plot_name,'_ninter',ninter,'.png'))
  ggsave(paste0(dir_out,'/',plot_name,'_ninter',ninter,'.pdf'), height = 10, width = 10, units = 'cm')
}


###########################################################################################################
# function that calculates the prey degree in Mass Spectrometry studies (used for detecting the prey hubs)
# input: dataframe of network (IntAct)
# output: dataframe including uniprotIDs, symbols, prey usage, degree of proteins 
#         considering only the interactions where the protein has been tested as a prey
###########################################################################################################

get_prey_massSpec <- function(human_data){
  
  human_data_pubmed <- unique(human_data[,c('IDs_interactor_A','IDs_interactor_B','Publication_Identifiers','Experimental_roles_interactor_A','Experimental_roles_interactor_B','Interaction_detection_methods')])
  # filter for specific methods
  index <- grep('coimmunoprecipitation|pull down|tandem affinity purification|mass spectrometry|copurification|affinity chromatography technology',human_data_pubmed$Interaction_detection_methods)
  human_data_pubmed_method <- human_data_pubmed[index,]
  # filter studies with more than 100 interactions
  human_data_bis <- unique(human_data_pubmed_method[,c('IDs_interactor_A','IDs_interactor_B','Publication_Identifiers')])
  counter_inter <- plyr::count(human_data_bis, vars = 'Publication_Identifiers')
  colnames(counter_inter)[2] <- 'freq'
  n <- 100
  studies <- counter_inter$Publication_Identifiers[which(counter_inter$freq >= n)] 
  
  human_data_pubmed_method <- human_data_pubmed_method[which(human_data_pubmed_method$Publication_Identifiers %in% studies),]
  preys <- get_bait_prey(human_data_pubmed_method,'prey')
  prey_usage <- calculate_bait_usage(human_data_pubmed_method,preys,'prey')
  #write.csv(prey_usage, file = paste0(dir_out,'/prey_usage_massSpec.csv'),row.names = F)
  return(prey_usage)
}

#################################################################################################################
# calculate degree of proteins in Y2H study
# input:
# - intact: dataframe of the network
# - s = pubmedID of interest
# output: dataframe that includes the degree for each protein (both gene symbol and uniprotID) in the study 
###################################################################################################################

calculate_degree_y2h_study <- function(intact,s){
  
  intact_y2h <- intact[grep('two hybrid',intact$Interaction_detection_methods),]
  intact_32296183 <- intact_y2h[which(intact_y2h$Publication_Identifiers == s),]
  tot_bait <- get_bait_prey(intact_32296183,'bait')
  length(tot_bait)
  tot_prey <- get_bait_prey(intact_32296183,'prey')
  length(tot_prey)
  # calculate degree
  g <- unique(intact_32296183[,c('IDs_interactor_A','IDs_interactor_B')])
  degree <- degree_wo_bidirEdges(g)
  degree$proteins <- rownames(degree)
  degree <- degree[which(degree$proteins %in% union(tot_bait,tot_prey)),]
  # convert uniprot to gene symbol
  symbol <- as.data.frame(mapIds(org.Hs.eg.db, keys = degree$proteins, keytype = "UNIPROT", column= "SYMBOL"))
  colnames(symbol)[1] <- 'symbol'
  degree$symbol <- symbol$symbol[match(degree$proteins,rownames(symbol))]
  degree <- degree[,c('proteins','symbol','degree')]
  colnames(degree)[1] <- 'uniprot'
  #write.csv(degree,file = paste0(dir_out,'/degree_Y2H_',s,'.csv'), row.names = F)
  return(degree)
}

##############################################################################
# clusterProfiler analyses: GO, DO, KEGG and Reactome 
# input:
# - genes and background= list of entrezIDs
# - n= number of hubs
# - label= string of the kind of approach we used to detect hubs
# - o= ontology (used only in Gene Ontology analysis)
# - dir_out= directory where the output of the enrichment analysis is saved
##############################################################################

clusterProfiler_analyses <- function(genes,background,p = 0.05,n,label,o,dir_out){
  
  edo <- enrichDO(genes,ont = "DO",qvalueCutoff = p, pvalueCutoff =p, pAdjustMethod = "fdr",universe = background,minGSSize = 2,readable= TRUE)
  #ekegg <- enrichKEGG(genes,organism = "hsa",keyType = "kegg",qvalueCutoff = p,pvalueCutoff =p,pAdjustMethod = "fdr",universe=background,minGSSize = 2)
  epath <- enrichPathway(genes,organism = "human",qvalueCutoff = p,pvalueCutoff =p,pAdjustMethod = "fdr",universe=background,minGSSize = 2,readable= TRUE)
  ego <- enrichGO(genes,universe=background,OrgDb= org.Hs.eg.db,ont= o,pAdjustMethod = "fdr",qvalueCutoff  = p,pvalueCutoff =p,readable= TRUE,minGSSize = 2)
  
  if(nrow(edo) != 0){
    write.csv(as.data.frame(edo),file = paste0(dir_out,'/enrichDO_results/table_enrichDO_',label,'_',n,'_qvalue',p,'.csv'),row.names = F)
  }
  # if(!is.null(ekegg)){
  #   if(nrow(ekegg) != 0){
  #     write.csv(as.data.frame(ekegg),file = paste0(dir_out,'/enrichKEGG_results/table_enrichKEGG_',label,'_',n,'_qvalue',p,'.csv'),row.names = F)
  #   }
  # }
  
  if(nrow(epath) != 0){
    write.csv(as.data.frame(epath),file = paste0(dir_out,'/enrichPathway_results/table_enrichPathway_',label,'_',n,'_qvalue',p,'.csv'),row.names = F)
  }
  if(nrow(ego) != 0){
    write.csv(as.data.frame(ego),file = paste0(dir_out,'/enrichGO_results/table_enrichGO_',label,'_',n,'_qvalue',p,'.csv'),row.names = F)
  }
}

##########################################################################################################################
# input: 
# - network: IntAct (it has to include the Experimental_roles_interactor_A and Experimental_roles_interactor_B columns
#   with the bait/prey annotation)
# - degree: dataframe with the degree and proteins columns
# output:
# dataframe with uniprotIDs, symbols, degree_bait (= number of interactions where the protein has been studied as a bait),
# degree_prey (= number of interactions where the protein has been studied as a prey), 
# total_degree (= number of total interactions in the network)
##########################################################################################################################

get_bait_prey_total_degree <- function(network,degree){
  
  baits <- get_bait_prey(network,'bait')
  bait_table <- calculate_bait_usage(network,baits,'bait')
  bait_table$bait_usage <- NULL
  colnames(bait_table)[c(1,2)] <- c('uniprot','symbol')
  
  preys <- get_bait_prey(network,'prey')
  prey_table <- calculate_bait_usage(network,preys,'prey')
  prey_table$prey_usage <- NULL
  colnames(prey_table)[c(1,2)] <- c('uniprot','symbol')
  
  final <- merge(bait_table,prey_table, all=T)
  final$degree_bait[is.na(final$degree_bait)] <- 0
  final$degree_prey[is.na(final$degree_prey)] <- 0
  #add total degree
  final$total_degree <- degree$degree[match(final$uniprot,degree$proteins)]
  return(final)
}
