# Identification of influential subjects in a network using a causal framework and evaluation of congruence with centrality measures

Youjin Lee, Ashley Buchanan, Elizabeth Ogburn, Samuel R. Friedman, M. Elizabeth Halloran, Natallia V. Katenka, Jing Wu, and Georgios Nikolopoulos


## Overview

Researchers across a wide array of disciplines are interested in identifying the most influential subjects in a network. In a network setting, intervention effects and health outcomes can spill over from one node to another through network ties, and influential subjects are expected to have a greater impact than others. For this reason, network research in public health has attempted to maximize health and behavioral changes by intervening a subset of influential subjects in a network. 

Although influence is often defined only implicitly in most of the literature, we found that the operative notion of influence is inherently causal in many cases: influential subjects are those we should intervene on in order to achieve the greatest overall effect across the entire network. In this work, we formally define a causal notion of influence using a potential outcome framework. We review existing influence measures, such as centrality, that largely rely on the particular features of the network structure and/or on certain diffusion models that predict the pattern of information or diseases spreads through network ties.

We provide simulation studies to demonstrate when popular centrality measures can agree with our causal measure of influence.
As an illustrative example, we apply several popular centrality measures to the Transmission Reduction Intervention Project and demonstrate the assumptions under which each centrality can capture the causal influence of each participant in the study.


## Data 

We use the Transmission Reduction Intervention Project (TRIP) data for the application study. The study design is detailed at [here](https://www.nature.com/articles/srep38100).


We provide the sample data `Data/sample_net.RData`, which has the same data structure (e.g., the same number of nodes and edges) as the TRIP data we analyzed in the manuscript. This data was generated at random and is only meant to be used as an example, not to be used to replicate the results.


## Code for Reproducibility

* `sim.R`  
This `R` file generates the main simulation study in Section 5. We consider three data generating models: (i) homogeneous direct interference, (ii) traffic-dependent process, and (iii) homogeneous diffusion process.For each model, we calculate out-degree, betweenness, and diffusion centrality and compare them to our influence measure $\tau$. 

* `read_sim.R`
This file reads the simulation results from `sim.R` and creates the tables and figures that demonstrate the congruence of each centrality to the causal measure of influence. 


* `sim_missing.R`
This `R` file runs the same simulation as in `sim.R` with random missing edges.The simulation results are presented in Supplementary Material. 

* `sim_multi.R` 
This file explores the difference between identifying a single most influential subjects and multiple influential subjects with $N=10$ nodes. The result is presented in Supplementary Material.

