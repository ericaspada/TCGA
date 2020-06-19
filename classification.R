# create train and test sets using the caret package
# split 65% and 35%
set.seed(17)
in_train = createDataPartition(y=data_fcbf$Class, p=.65, list=FALSE) # randomly samples the data (default) and returns the indexes 
train = data_fcbf[in_train,]
test = data_fcbf[-in_train,]
test[1:25,1:25]
train$Class = as.factor(train$Class)

# experiment with cross-validation
# createFolds(y=merged_r$Class, k=5, list=FALSE)

# Implement SVM using the e1071 package
library(e1071)

# by default, the svm function scales the data to mean 0 and variance 1
svm_fit_1 = svm(Class~., data = train, kernel = "linear", cost=10) # linear decision boundary
svm_fit_2 = svm(Class~., data = train, kernel = "radial", cost=10) # default gamma=1/ncol
svm_fit_3 = svm(Class~., data = train, kernel = "radial", cost=1) # smaller C
svm_fit_4 = svm(Class~., data = train, kernel = "radial", cost=10, gamma=.1) # larger gamma
svm_fit_5 = svm(Class~., data = train, kernel = "radial", cost=10, gamma=.001) # smaller gamma

# define a function to evaluate the fitted models
eval_svm = function(x){
  predictions = predict(x, newdata=test[,2:99])
  df1 = as.data.frame(predictions)
  df2 = as.data.frame(test[,1])
  compare = cbind(df1, df2)
  compare$correct = (compare[,1]==compare[,2])
  precision = count(filter(compare,correct==TRUE))/count(compare)
  confusion_matrix = table(predicted=predictions,actual=test$Class)
  print(precision)
  print(confusion_matrix)
}

# evaluate the model fit with linear kernel function
eval_svm(svm_fit_1)
# evaluate the model fit with RBF kernel function, C=10, gamma=1/ncol
eval_svm(svm_fit_2)
# evaluate the model fit with RBF kernel function, C=1, gamma=1/ncol
eval_svm(svm_fit_3)
# evaluate the model fit with RBF kernel function, C=10, gamma=.1
eval_svm(svm_fit_4)
# evaluate the model fit with RBF kernel function, C=10, gamma=.001
eval_svm(svm_fit_5)
# models 1 and 5 have the highest precision (1 error)

  