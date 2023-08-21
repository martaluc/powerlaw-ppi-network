# Emergence of power-law distributions in protein-protein interaction networks through study bias

Scripts to reproduce our analyses of the power-law property in single study and aggregated protein-protein interaction (PPI) networks.

### Clone the repository

```
git clone https://github.com/martaluc/powerlaw-ppi-network.git
```

### Download data from Zenodo
Download [here](https://doi.org/10.5281/zenodo.7695121) the required data (databases.zip), move it in the GitHub repository and unzip it.

### Docker
1) Download Docker, following the instructions described [here](https://docs.docker.com/engine/). Keep in mind that you will need root permissions to both install and run Docker. Then you can pull the Docker image from [Docker Hub](https://hub.docker.com/r/ma10r02t90a/powerlaw-ppi):
```
docker pull ma10r02t90a/powerlaw-ppi
```
2) After pulling the `ma10r02t90a/powerlaw-ppi` image, you can run the Docker container interactively using the following command:
```
sudo docker run -it -e DISPLAY --name ppi -v /your/path/to/powerlaw-ppi-network/:/powerlaw-ppi-network ma10r02t90a/powerlaw-ppi /bin/bash
```
`-it` specifies an interactive terminal, allowing you to interact with the container.

`--name ppi` assigns the name "ppi" to the container. This flag can be omitted if you prefer.

`-v /your/path/to/powerlaw-ppi-network/:/powerlaw-ppi-network` mounts the directory from your local system (replace `/your/path/to/powerlaw-ppi-network/` with the actual path) into the container's `/powerlaw-ppi-network` directory. This allows you to share files between your local system and the container.

`ma10r02t90a/powerlaw-ppi` specifies the image name that you pulled from Docker Hub.

`/bin/bash` launches a Bash shell inside the container, providing you with an interactive terminal.


3) Once you are inside the container, move to the folder of interest with `cd /powerlaw-ppi-network` and run the R scripts with the following command:

```
Rscript <script_name.R>
```
Replace `<script_name.R>` with the actual name of the R script you want to run.

### Description
The code was run on Ubuntu (20.04.5 LTS) using R (4.2.2).

1-aggregated_network.R calculates the degree and the bait distribution of our aggregated PPI network (Figure 2A and Figure 3B) and the distribution of publication number (Figure 3C). It also tests the power-law property of each single study (Figure 2B).

2-correcting_baitUsage.R analyzes if the asymmetric design (i.e. number of baits and preys) affects the power-law property (Figure 4B and C).

3-merging_studies.R tests if the power-law property emerges from the aggregation of non-power law studies (Figure 3A).

4-hubs.R detects the hubs and performs the functional and disease enrichment analyses (Figure 5).

The results of the scripts are saved in the `output/` folder and the plots in the `plots/` folder.


##### Notes
• Before running the scripts, make sure to review the number of CPUs specified in the scripts. If needed, modify the CPUs settings based on the available computational resources on your machine.

• Whenever possible, it is recommended to use the exact number of threads specified in the code. This ensures the consistency of p-values, particularly when testing for the power law property.
