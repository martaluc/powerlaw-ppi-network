# powerlaw-ppi-network

Scripts to reproduce our analyses of the power-law property in single study and aggregated protein-protein interaction (PPI) networks.

#### Clone the repository

```
git clone https://github.com/martaluc/powerlaw-ppi-network.git
```

#### Download data from Zenodo
Download [here](https://doi.org/10.5281/zenodo.7695121) the required data (databases.zip), put it in the GitHub repository and unzip it.

#### Description
The code was run on Ubuntu (20.04.5 LTS) using R (4.2.2).

1-aggregated_network.R calculates the degree and the bait distribution of our aggregated PPI network (Figure 2A and Figure 3B). It also tests the power-law property of each single-study (Figure 2B).

2-correcting_baitUsage.R analyzes if the asymmetric design (i.e. number of baits and preys) affects the power-law property (Figure 4B and C).

3-merging_studies.R tests if the power-law property emerges from the aggregation of non-power law studies (Figure 3A).

4-hubs.R calculates the hubs and performs the functional and disease enrichment analyses (Figure 5).


###### Notes
Use the same number of threads as specified in the code to obtain the exact same p-values (when testing the power law property).
