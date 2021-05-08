# Personalized-Disease-Diagnosis-Assistants
LOADDx and SCADDx can produce a short personalized ranked list of the most likely diseases for each patient at a requested time-point, from a knowledge base with thousands of diseases. They enable personalized disease diagnosis by combining patients’ temporal gene expression data with a Knowledge Graph. The proposed algorithms can predict the most likely diseases for each patient along with their most affected genes, which can support health care professionals in their decision-making.

## How to Install and Run the Software
LOADDx and SCADDx is written in R progrmming language, and can run on Windows, MacOS or Linux. However, the R programming language is required to be installed before you run the software/code.

### Required Programming Language:
R version 3.6.2 or above

You can download the latest version of R from here:
* [Install R](https://www.r-project.org/)


### Steps to run the code:
1. Download the provided R code, gene expression datasets and knowledge graph, and keep them in the same folder. 
2. Open the terminal
3. Go to the folder where you downloaded all the code and datasets
4. Run SCADDx using this command: 
```
R CMD BATCH SCADDx.R
```
5. Run LOADDx using this command: 
```
R CMD BATCH LOADDx.R
```
These command will create a .Rout (output) file in the same folder. This .Rout file will contain all the results. 
