
# Reading data
getwd();
InitialData <- read.csv(file="gene-expression-invasive-vs-noninvasive-cancer.csv")

# Random index subset of of 2000 genes (columns)
set.seed(2315782)
team.gene.subset <- rank(runif(1:4948))[1:2000] 

genes <- InitialData[,team.gene.subset]
head(genes)


## Part 1: dimensionality reduction

#note 1: (missing values)
#We do bootstrap to replace the missing values assuming NA are data missing at random, which makes sense as this is a dataset of genes from patients.In bootstrap we take the variable column as our sampling set (all non NA variables) and then replace the NAs with this values. Bootstrap is more stable than mean imputation when having a small amount of NAs because mean is sensitive to outliers and influential points.

#note 2: PCA (unsupervised dim red)
#princomp uses jus the covariance matrix and assumes that the **number of observations is greater than the number of variables** and that we have a fairly square matrix. 
#prcomp function works here because it uses a different algorithm to princomp, prcomp is based on singular value decomposition (SVD) of our covariance matrix, which is a more general method for computing principal components.

# MISSING DATA

# Function that replaces the missing values with bootstrap samples (also checks if the number of items to replace is not a muiltiple of replacement length, i.e. columns without NA's)

replace_missing <- function(x) {
  non_missing <- x[!is.na(x)]
  if (length(non_missing) > 0) {
    bootstrap_sample <- sample(non_missing, replace = TRUE, size = sum(is.na(x)))
    x[is.na(x)] <- bootstrap_sample
  }
  return(x)
}

x = is.na(genes)
sum(x)
anyNA(genes)


bootstrap_genes <- as.data.frame(lapply(genes, replace_missing))
dim(bootstrap_genes)

x = is.na(bootstrap_genes)
sum(x)
anyNA(bootstrap_genes)


# we center and scale the values:
# center: subtracting the mean of each variable to all its values
# scale: dividing each variable value by its standard deviation 

# NOTE: PCA might not be a great way of reducing dimensionality as the resulting variables, i.e. PC's, are linear combinations of the variables, so information of the orginal variables is lost.

# DIMENSIONALITY REDUCTION (unsupervised)

# PCA
pca_result <- prcomp(bootstrap_genes, center = TRUE, scale = TRUE)
# Visualize the variance explained by each principal component
par(mfrow=c(1,1))
plot(pca_result)

pcs <- pca_result$x

summary(pca_result)
# Plot of cumulative proportion of variance explained

plot(cumsum(pca_result$sdev^2) / sum(pca_result$sdev^2), 
     xlab = "Number of Principal Components", 
     ylab = "Cumulative Proportion of Variance Explained",
     type = "b")

abline(v = 60, col = "red", lty = 2)
abline(h = 0.95, col = "blue", lty = 2)
text(65, 0.8, "60 PC's", col = "red", pos = 4)
text(10, 0.96, "95%", col = "blue", pos = 2)

# Principal components that explain 95% of the variance 
num_components <- which(cumsum(pca_result$sdev^2) / sum(pca_result$sdev^2) >= 0.95)[1]
# 60

selected_pcs <- pcs[, 1:num_components]

dim(selected_pcs)

new_genes_pca <- as.data.frame(selected_pcs)



# We tried reducing the variable set using feature selection, specifically backwards selection because of its efficiency, however due to the high dimensionality and colinearity of our dataset it did not yield a great dimensional reduction and it was computationally expensive as well.



# WARNING JUST RUN ONCE  !!!!!!!!
# We are retrieve our bootstrapped dataframe of variables and we create another dataframe called genes_full that includes the class response variable on the first column
# adding back class labels (new column)

dim(bootstrap_genes)
genes_full <- bootstrap_genes

genes_full$new_column <- InitialData[,4949] #here we are creating a dataset that includes the class variable like [class response variable, all other variables]

#1 if it's invasive, 2 otherwise
genes_full$new_column <- as.factor(genes_full$new_column)

genes_full <- genes_full[-54, ] #obs 54 is deleted because it has more than 70% missing data

dim(genes_full)




# We will use RFEs in order to reduce the dimensionality of our database and keep significant variables without loosing information like PCA does. Our aim is to have approximately the same number of variables as observations to perform ML agorithms in questions 2 and 3, if we have more variables than observations this will cause underfitting in our training set, however if we have more observations than variables this might cause overfitting on our training set.
# As we have 77 we are aiming for a reduction of about 77 variables (or a number close to it, preferably less). This method is computationally expensive but not as glm.  

# In recursive random forests we do this steps

# 1. Rank the importance of all features using backwards selection (start with all features and then remove the least important ones iteratively based on their significance yielded after doing t test on each variable).
# 2. Eliminate the least important features. (removes small amount of varibales per loop that have collinearity)
# 3. Build a model using the remaining features.
# 4. Repeat steps 1-3 until the desired number of features is reached.

# RFE trains the model during each iteration and evaluates its performance using a cross-validation set.


# DIMENSIONALITY REDUCTION SUPERVISED

# RANDOM FOREST RECURSIVE SELECTION
library("randomForest")
library("caret")
set.seed(2315782)

#Splitting Data
#use 70% of dataset as training set and 30% as test set
sample_rfe <- sample(c(TRUE, FALSE), nrow(genes_full), replace=T, prob=c(0.7,0.3))

xtrain <- genes_full[sample_rfe, -ncol(genes_full)]
dim(xtrain)
#(genes_full[, 2001])

xtest <- genes_full[!sample_rfe, -ncol(genes_full)]
dim(xtest)

ytrain  <- genes_full$new_column[sample_rfe]
ytrain

ytest   <- genes_full$new_column[!sample_rfe]
ytest
#RFE Supervised Dimension Reduction
# Define the control using a random forest selection function
control <- rfeControl(functions = rfFuncs, # random forest
                      method = "repeatedcv", # cv to test performance
                      repeats = 5, # number of iterations
                      number = 10) # number of folds in cv set

# Run RFE
result_rfe1 <- rfe(x = xtrain, 
                   y = ytrain, 
                   sizes = c(1:40),
                   rfeControl = control) # we aim for size of 77 variables

# Features selected by rfe
predictors(result_rfe1)

# Print the results
g1 <- result_rfe1$results #preparation for better visualization
dim(g1)

selected_variables <- result_rfe1$optVariables
print(selected_variables)

summary(selected_variables)
new_genes_subset <- genes_full[, c(selected_variables, "new_column")]
new_genes_subset
dim(new_genes_subset)


## Part 2: unsupervised 

#### K-means clustering

# To choose the hyperparameter k, the number of clusters, we will compute the sum of squared distances between each data point and each assigned cluster center. Where sum of the square (euclidean) distances can be computed like:

# $\sum_{i=1}^{k} \sum_{j=1}^{n_i} \left\| x_{ij} - c_i \right\|^2$

# We checkif there is an "elbow" where the sum of squares starts to decrease in size and we take that number of clusters. In this case 4.

# note: If the data were to follow a different distribution than normal (qqplot checking) we need to check it's empirical distribution to check if they follow another known distribution, and then assess which distance to choose depending on MLE of their variance.


# K-MEANS CLUSTERING

# K-MEANS CLUSTERING

new_genes_subset_wo <- new_genes_subset[, -ncol(new_genes_subset)] # our reduced dataset without the class variable "invasive"

# K-MEANS FOR GENES

# 1. Choosing  hyperparameter k:
library(ggplot2)


# Computes sum of squared distances between each  datapoint and each assigned cluster center
k_values <- 1:10
ssq <- numeric(length(k_values))

for (k in k_values) {
  kmeans_model <- kmeans(new_genes_subset_wo, centers = k)
  ssq[k] <- kmeans_model$tot.withinss
}

# Elbow curve plot
ggplot(data.frame(k = k_values, ssq = ssq), aes(x = k, y = ssq)) +
  geom_line() +
  geom_point() +
  labs(x = "Number of Clusters (k)", y = "Sum of squared distances from a point to centroid") +
  theme_minimal()

#K-means for genes
library(fpc)
clusters <- kmeans(new_genes_subset_wo, centers = 3, nstart = 10)
plotcluster(new_genes_subset_wo, clusters$cluster, xlab = "Clusters of genes")

# 3. Hierarchical clustering
clust <- hclust(dist(new_genes_subset_wo))
plot(clust, xlab = "Clusters", y = "Height of Clusters")


# K-MEANS FOR PATIENTS

# 1. Choosing  hyperparameter k:
library(ggplot2)


# Compute sum of squared distances between each  datapoint and each assigned cluster center
k_values <- 1:10
ssq <- numeric(length(k_values))

for (k in k_values) {
  kmeans_model <- kmeans(t(new_genes_subset_wo), centers = k)
  ssq[k] <- kmeans_model$tot.withinss
}

# Elbow curve plot
ggplot(data.frame(k = k_values, ssq = ssq), aes(x = k, y = ssq)) +
  geom_line() +
  geom_point() +
  labs(x = "Number of Clusters (k)", y = "Sum of squared distances from a point to centroid") +
  theme_minimal()


# K-means for patients
new_genes_subset_pat <- t(new_genes_subset_wo)
clusters <- kmeans(new_genes_subset_pat, centers = 3, nstart = 10)
plotcluster(new_genes_subset_pat, clusters$cluster, xlab = "Clusters of patients")

# 3. Hierarchical clustering
clust <- hclust(dist(t(new_genes_subset_wo)))
plot(clust, xlab = "Clusters", y = "Height of Clusters")


# K-MEANS for genes
# 
# We notice that the elbow curve has an edge at 3. So we choose our hyperparameter k to be 3..Each cluster represents a group of patients who have similar patterns of gene expression across the 55 variables.
# 
# We see that in the dendogram that objects located at the same height share similar characteristics and visually the height when there is three defined clusters starts getting bigger and less patients start joining to the cluster. Observation 53 is joins very late to them and it is dissimilar to the other ones because of the high number of NAs the observation contained.
# 
# K-MEANS for patients
# 
# We notice that the elbow curve for patients has no defined edge, however we will chose 3 centers because that seems like a moderate edge. Each gene is now treated as a data point, and the clustering algorithm groups genes together based on similarities in their expression patterns across the 77 patients.
# 
# When we plot the kmeans clustering for patients we notice that almost all the datapoints overlap, this is a reasonable assumption for two reasons.
# 
# 1. Most data in a row might be similar or equal because gene expression levels might coincide between patients because there might be an interval with a very small variance of that gene level before it being considered abnormal.
# 
# 2. The patients are independent from each other, unlike the genes where they shared high levels of colinearity still after the dimensionality was reduced.
# 
# We also see that in the dendogram here there seems to be many similar genes that start joining fast to clusters and then it slows down. 
# ***Note that this can also be a good method of reducing dimensionality too, but we would end up keeping 8 variables and we run the risk of overfitting with very few genes variables and a lot of patients observations.



## Part 3: supervised


# RANDOM FORESTS
library(randomForest)

# We divide our labeled dataset into train (70% of data) and test data (30% of data)

ind <- sample(2, nrow(new_genes_subset), replace = TRUE, prob = c(0.7, 0.3))
#ind <- sample(1:nrow(new_genes_subset), 0.8 * nrow(new_genes_subset))
train_data <- new_genes_subset[ind==1,]
#train_data$Label<-factor(train_data$new_column,label=c("invasive","noninvasive"))
dim(train_data)

test_data <- new_genes_subset[ind==2,]
dim(test_data)

rf_model <- randomForest(new_column~., 
                         data = train_data, 
                         mtry = 3, importance = TRUE
                         #proximity=TRUE
) 
print(rf_model)

importance(rf_model)
varImpPlot(rf_model)

predicted_classes <- predict(rf_model, newdata = test_data)

## Performace 
confusion_matrix <- confusionMatrix(data = predicted_classes, reference = test_data$new_column)

# Extract performance metrics
accuracy <- confusion_matrix$overall['Accuracy']
kappa <- confusion_matrix$overall['Kappa']
precision <- confusion_matrix$byClass['Precision']
recall <- confusion_matrix$byClass['Recall']
f1_score <- confusion_matrix$byClass['F1']

# Print the performance metrics
cat("Accuracy:", accuracy, "\n")
cat("Kappa:", kappa, "\n")
cat("Precision:", precision, "\n")
cat("Recall:", recall, "\n")
cat("F1 Score:", f1_score, "\n")
##





## Logistic Model


dim(new_genes_subset)
train_indices <- sample(1:nrow(new_genes_subset), 0.7 * nrow(new_genes_subset))  # 80% for training
train_data <- new_genes_subset[train_indices, ]
test_data <- new_genes_subset[-train_indices, ]

train_labels <- new_genes_subset$new_column[train_indices] 
test_labels <- new_genes_subset$new_column[-train_indices] 

number_of_columns <- ncol(train_data) - 1

logistic_model <- glm(train_labels~. , data = train_data[, 1:number_of_columns], family = binomial)
summary(logistic_model)

predicted_probabilities <- predict(logistic_model, 
                                   newdata = test_data, type = "response")


predicted_classes <- ifelse(predicted_probabilities <  0.5, 1, 2)  # Threshold at 0.5
predicted_classes

# Evaluate Model Performance

conf_matrix <- table(Actual = test_labels, Predicted = predicted_classes)
print("Confusion Matrix:")
print(conf_matrix)

accuracy <- sum(diag(conf_matrix)) / sum(conf_matrix)
print(paste("Accuracy:", accuracy))

precision <- conf_matrix[2, 2] / sum(conf_matrix[, 2])
print(paste("Precision:", precision))

recall <- conf_matrix[2, 2] / sum(conf_matrix[2, ])
print(paste("Recall:", recall))

f1_score <- 2 * precision * recall / (precision + recall)
print(paste("F1-score:", f1_score))

# Create Data Frame for Metrics
metrics_df <- data.frame(Metric = c("Accuracy", "Precision", "Recall", "F1-Score"),
                         Value = c(accuracy, precision, recall, f1_score))

ggplot(metrics_df, aes(x = Metric, y = Value)) +
  geom_bar(stat = "identity", fill = "skyblue") +
  labs(title = "Performance Metrics of Logistic Model",
       y = "Value") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


## LDA

library(MASS)
linear_training <- lda(train_data$new_column~., train_data[, 1:number_of_columns])
#linear_training
attributes(linear_training)

p <- predict(linear_training, train_data[, 1:number_of_columns])
summary(p)
p$x
ldahist(data = p$x[,1], g = train_data$new_column)


## LDA Accuracy
train_accuracy <- mean(p$class == train_data$new_column)
train_precision <- sum(p$class == 1 & train_data$new_column == 1) / sum(p$class == 1)
train_recall <- sum(p$class == 1 & train_data$new_column == 1) / sum(train_data$new_column == 1)
train_f1_score <- 2 * (train_precision * train_recall) / (train_precision + train_recall)

cat("Training Accuracy:", train_accuracy, "\n")
cat("Training Precision:", train_precision, "\n")
cat("Training Recall:", train_recall, "\n")
cat("Training F1 Score:", train_f1_score, "\n")

###

linear_testing <- lda(test_data$new_column~., test_data[,1:number_of_columns])
#linear_testing
attributes(linear_testing)

# Histogram Testing
pt <- predict(linear_testing, test_data[,1:number_of_columns])
pt$x
ldahist(data = pt$x[,1], g = test_data$new_column)

### LDA Testing Accuracy
test_accuracy <- mean(pt$class == test_data$new_column)
test_precision <- sum(pt$class == 1 & test_data$new_column == 1) / sum(pt$class == 1)
test_recall <- sum(pt$class == 1 & test_data$new_column == 1) / sum(test_data$new_column == 1)
test_f1_score <- 2 * (test_precision * test_recall) / (test_precision + test_recall)

cat("Testing Accuracy:", test_accuracy, "\n")
cat("Testing Precision:", test_precision, "\n")
cat("Testing Recall:", test_recall, "\n")
cat("Testing F1 Score:", test_f1_score, "\n")
##




library(MASS)

# LDA
lda.invasive <- lda(new_column ~ . , data=new_genes_subset)

# Predicts the class label of each observation based on training data 
invasive.predicted <- predict(lda.invasive)$class 

# Contingency table 
table(Predicted = invasive.predicted, Actual = new_genes_subset$new_column) 


# QDA
qda.invasive <- qda(new_column ~ ., data = new_genes_subset)

qda.invasive <- qda(new_column ~ . , data=new_genes_subset)

# Predicts the class label of each observation based on training data 
invasive.predicted2 <- predict(qda.invasive)$class 

# Contingency table
table(Predicted = invasive.predicted2, Actual = new_genes_subset$new_column) 




## KNN

library(caret)
library(class)

dim(train_data)
summary(train_labels)

knn_model <- train(train_data[, 1:number_of_columns], train_labels, 
                   method = "knn")
summary(knn_model)

# Get predictions on training data
train_predictions <- predict(knn_model, newdata = train_data[, 1:number_of_columns])

# Evaluate performance
accuracy <- confusionMatrix(train_predictions, train_labels)$overall["Accuracy"]
precision <- confusionMatrix(train_predictions, train_labels)$byClass["Precision"]
recall <- confusionMatrix(train_predictions, train_labels)$byClass["Recall"]
f1_score <- confusionMatrix(train_predictions, train_labels)$byClass["F1"]

# Confusion matrix
confusion_matrix <- confusionMatrix(train_predictions, train_labels)
print(confusion_matrix)
# Apply KNN model End


library(e1071)
library(caret)

svm_model <- svm(train_labels ~ ., data = train_data[, 1:number_of_columns],
                 kernel = "radial")
summary(svm_model)

predicted_classes <- predict(svm_model, newdata = test_data)
summary(predicted_classes)

tbl <- table(Actual = test_labels, Predicted = predicted_classes)
print(tbl)




cf <- caret::confusionMatrix(data=predicted_classes,
                             reference=test_labels)
print(cf)
# Extract Metrics
accuracy <- cf$overall["Accuracy"]
precision <- cf$byClass["Precision"]
recall <- cf$byClass["Sensitivity"]
f1_score <- cf$byClass["F1"]

# Create Data Frame for Metrics
metrics_df <- data.frame(Metric = c("Accuracy", "Precision", "Recall", "F1-Score"),
                         Value = c(accuracy, precision, recall, f1_score))

# Plot Bar Graph
ggplot(metrics_df, aes(x = Metric, y = Value)) +
  geom_bar(stat = "identity", fill = "skyblue") +
  labs(title = "Performance Metrics of SVM after PCA",
       y = "Value") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


# Perform Cross-Validation
# Discussion
number_of_columns <- ncol(train_data) - 1 #
print(number_of_columns)

dim(train_data)
?train
train_data[, 26]

library(caret)
# Define your control parameters (optional)
#ctrl <- trainControl(method = "cv", number = 10)  # Example of 10-fold cross-validation
ctrl <- trainControl(method = "cv",
                     number = 10,
                     summaryFunction = twoClassSummary,
                     classProbs = TRUE, 
                     # Specify performance measures
                     # Here, 'Accuracy', 'Kappa', 'MAE', 'RMSE', 'Rsquared' will be calculated
                     savePredictions = TRUE,
                     verboseIter = TRUE)

# Train the SVM model using the train function
#svm_cv <- train(x = train_data[, 1:number_of_columns],  # Predictors (features)
#                y = as.numeric(train_labels),          # Outcome (class labels)
#                method = "svmLinear",                  # SVM with linear kernel
#                trControl = ctrl,                      # Cross-validation settings
#                preProcess = c("center", "scale"),      # Pre-processing steps
#                tuneGrid = expand.grid(C = seq(0, 2, length = 20))  # Hyperparameter tuning (optional)
#)
# Train the SVM model using the train function
train_labels_modified <- train_labels
levels(train_labels_modified) <- make.names(levels(train_labels_modified))
train_labels_modifiedf <- as.factor(train_labels_modified) # changing variable names to X1 and X2

svm_cv <- train(x = train_data[, 1:number_of_columns],  # Predictors (features)
                y = train_labels_modifiedf,          # Outcome (class labels)
                method = "svmLinear",                  # SVM with linear kernel
                trControl = ctrl,                      # Cross-validation settings
                preProcess = c("center", "scale")      # Pre-processing steps
)
summary(svm_cv)


svm_model <- svm(x = train_data[, 1:number_of_columns],  # Predictors (features)
                 y = train_labels_modifiedf,            # Outcome (class labels)
                 kernel = "linear",                      # Linear kernel
                 scale = TRUE)                           # Scale the data


# Make predictions on the scaled test data
predictions <- predict(svm_model, test_data)

# Calculate accuracy
accuracy <- sum(predictions == test_labels) / length(test_labels)

# Calculate confusion matrix and accuracy
conf_mat <- confusionMatrix(predictions, train_labels)
accuracy <- conf_mat$overall['Accuracy']

# Print accuracy
print(accuracy)
summary(svm_cv)

#svm_model <- svm(as.numeric(train_labels), data = train_data[, 1:number_of_columns], kernel = "linear")
#summary(svm_model);
train_labels_modified <- make.names(train_labels)

rf_cv <- train(train_data[,1:number_of_columns], train_labels_modified, method = "rf", trControl = ctrl)
summary(rf_cv);

knn_cv <- train(train_data[,1:number_of_columns], train_labels_modified, method = "kknn", trControl = ctrl)
summary(knn_cv);

# Compare Performance
compare_models <- resamples(list(
  SVM = svm_cv, 
  Random_Forest = rf_cv, 
  KNN = knn_cv))
summary(compare_models)


