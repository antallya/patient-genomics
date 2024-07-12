# Description
Hi! I am an applied mathematician/statistician wanting to learn more about the world of DS and CS.
This is a project I did during one of my master's classes where we analyzed a **single-cell RNA-seq dataset**.
My part in this project was to do the dimensionality reduction and data processing.

## Dataset Overview
The dataset used is a subset of a larger single-cell RNA-seq dataset compiled by the **Allen Institute**, I retrieved it from Kaggle
This data frame consists of 78 rows (patient's cell) 4949 columns with 4948 gene expression measurements of cancer tissue, each column representing a ‘gene’ with column 4949 having the information of a class variable with two values: 1 and 2 (invasive or non-invasive cancer).

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
  
  #### PCA
  PCA might not be a great way of reducing dimensionality as the resulting variables (PCs) are linear combinations of the variables, so information from the original variables is lost.

  #### LDA
  Didn't work well due to the high dimensionality and colinearity between each gene.
  
  #### Feature Selection using regsubsets()
  We tried reducing the variable set using feature selection, specifically using the regsubsets() function backward selection because of its efficiency, however, since it is based on linear regression and our variables
  are categorical, I didn't use it.

  #### GLM
  GLM was computationally expensive as it had to check pairwise each gene.

### Chosen Method: RFEs
Its recursive aspect is less computationally expensive.
We aim to have approximately the same number of variables as observations or less, if we have more variables than observations this will cause underfitting in our training set, however, if we have more observations than variables this might cause overfitting on our training set. (note that because we were required to use QDA, we aimed for fewer variables in our training set)

In recursive random forests, we do these steps:

1. Rank the importance of all features using backward selection (start with all features and then remove the least important ones iteratively based on their significance yielded after doing t test on each variable).
2. Eliminate the least important features. (removes a small amount of variables per loop that have collinearity)
3. Build a model using the remaining features.
4. Repeat steps 1-3 until the desired number of features is reached.

*RFE trains the model during each iteration and evaluates its performance using a cross-validation set.*




