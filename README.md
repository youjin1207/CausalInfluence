# Finding influential subjects in a network using a causal framework

Youjin Lee, Ashley Buchanan, Elizabeth Ogburn, Samuel R. Friedman, M. Elizabeth Halloran, Natallia V. Katenka, Jing Wu, and Georgios Nikolopoulos


## Overview

Researchers across a wide array of disciplines are interested in finding the most influential subjects in a network. In a network setting, intervention effects and health outcomes can spill over from one node to another through network ties, and influential subjects are expected to have a greater impact than others. For this reason, network research in public health has attempted to maximize health and behavioral changes by intervening on a subset of influential subjects. 

Although influence is often defined only implicitly in most of the literature, the operative notion of influence is inherently causal in many cases: influential subjects are those we should intervene on to achieve the greatest overall effect across the entire network. In this work, we define a causal notion of influence using potential outcomes. We review existing influence measures, such as node centrality, that largely rely on the particular features of the network structure and/or on certain diffusion models that predict the pattern of information or diseases spreads through network ties.

We provide simulation studies to demonstrate when popular centrality measures can agree with our causal measure of influence. As an illustrative example, we apply several popular centrality measures to the HIV risk network in the Transmission Reduction Intervention Project and demonstrate the assumptions under which each centrality can represent the causal influence of each participant in the study.

## Data 

We use the Transmission Reduction Intervention Project (TRIP) data for the application study. The study design is detailed at [here](https://www.nature.com/articles/srep38100).
The data are available on request from Dr. Georgios Nikolopoulos (nikolopoulos.georgios@ucy.ac.cy).  


We provide the sample data `Data/sample_net.RData`, which has the same data structure (e.g., the same number of nodes (277) and edges (542)) as the TRIP data we analyzed in the manuscript. This data was generated at random and is only meant to be used as an example, not to be used to replicate the results.


## Code for Reproducibility

* `aux_functions.R`
This file contains auxiliary functions to generate different diffusion processes.

* `sim.R`  
This `R` file generates the main simulation study in Section 5. We consider three data generating models: (i) homogeneous direct interference, (ii) traffic-dependent process, and (iii) homogeneous diffusion process. For each model, we calculate out-degree, betweenness, and diffusion centrality and compare them to our influence measure $\tau$. 

* `read_sim.R`
This file reads the simulation results from `sim.R` and creates the tables and figures that demonstrate the congruence of each centrality to the causal measure of influence. This code can be used to reproduce Table 2 and Figure 3 in the main manuscript and Figures S2 and S3 in the Supporting Information.


* `sim_multi.R` 
This file explores the difference between identifying a single most influential subjects and multiple influential subjects with $N=10$ nodes. The result is presented in the Supporting Information. This code can be used to reproduce Figure S4 and Tables S1, S2, and S3 in the Supporting Information. 


* `sim_missing.R`
This `R` file runs the same simulation as in `sim.R` with random missing edges.The simulation results are presented in the Supporting Information. This code can be used to reproduce Table S4 in the Supporting Information.


* `sim_sensitivity.R`
This `R` file runs the simulation with unit-varying $\alpha_{i}$ and $\beta_{ij}$ and random errors $\epsilon_{i} \overset{i.i.d}{\sim} N(0, \sigma^{2})$ from the following structural model (Equation (3) in the main manuscript) 

$$Y_{i}(\mathbf{a}_{\{k\}}) = \delta_{i} + \alpha_{i} I(a_{i} = a_{k}) + \sum_{j=1, j \neq i}^{N} \beta_{ji} I(a_{j} = a_{k}) + \epsilon_{i}, ~ k \in \{1,2,\ldots, N\}$$

This code can be used to reproduce Tables S6, S7, and S8 in the Supporting Information.

* `read_TRIP.R`
This `R` file is used to analyze the TRIP data, calculating the descriptive statistics and the centrality measures. For the purpose of illustration, we use `Data/sample_net.RData`. This code can be used to reproduce Figures 1 and 4 and Table 3 based on the pseudo data. 

