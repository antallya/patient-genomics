# Description
Hi! I am an applied mathematician/statistician wanting to learn more about the world of DS and CS.
This is a project I did during one of my master's classes where we analyzed a single-cell RNA-seq dataset.
My part in this project was to do the dimensionality reduction and data processing.

## Dataset Overview
The data set has 78 rows (patients) 4949 columns with 4948 gene expression measurements of cancer tissue, each column representing a ‘gene’ with column 4949 having the information of a class variable with two values: 1 and 2.

## Missing values
We use Bootstrap to replace the missing values assuming NA are data *missing at random*, which makes sense as this is a dataset of genes from patients. 
In Bootstrap we take the variable column as our sampling set (all non-NA variables) and then replace the NAs with these values. 
Bootstrap is more stable than mean imputation when having a small amount of NAs because the mean is sensitive to outliers and influential points.

## Data Standardization
Due to the presence of genes with extremely high magnitudes compared to each other, we standardize the data (center and scale)
* center: subtracting the mean of each variable from all its values
* scale: dividing each variable value by its standard deviation

## Dimensionality reduction

### Methods that were not optimal
  
  ### PCA
  PCA might not be a great way of reducing dimensionality as the resulting variables (PCs) are linear combinations of the variables, so information from the original variables is lost.
  
  ### Feature Selection
  We tried reducing the variable set using feature selection, specifically backward selection because of its efficiency, however, due to the high dimensionality/colinearity and it being more adept for continuous class
  variables, I didn't choose it.

  ### GLM

### Chosen Method


