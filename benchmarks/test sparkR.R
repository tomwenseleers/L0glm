# Set the system environment variables
Sys.setenv(SPARK_HOME = "C:/spark/spark-3.4.1-bin-hadoop3")
.libPaths(c(file.path(Sys.getenv("SPARK_HOME"), "R", "lib"), .libPaths()))

# for windows installation instructions of spark see
# https://phoenixnap.com/kb/install-spark-on-windows-10

#load the Sparkr library
library(SparkR)

# Create a spark context and a SQL context
sc <- sparkR.session(master = "local")

# demo linear regression ####
path <- file.path(Sys.getenv("SPARK_HOME"), "data/mllib/sample_linear_regression_data.txt")
df <- read.df(path, source="libsvm")

# fit Linear Regression Model
model <- spark.lm(df, label ~ features, regParam = 0.01, maxIter = 1,
                  solver="auto") # auto, l-bfgs or normal
# argument weightCol for weights

# get the summary of the model
summary(model)

# coefficients of the model
summary(model)$coefficients

# make predictions
predictions <- predict(model, df)

# save and load the model
path <- "path/to/model"
write.ml(model, path)
savedModel <- read.ml(path)
summary(savedModel)


# demo glm ####
# Create the DataFrame
df <- createDataFrame(iris)

# Fit a linear model on the dataset
system.time(model <- spark.glm(Sepal_Length ~ Sepal_Width + Species, 
                         data = df, 
                         family = gaussian(identity),
                         stringIndexerOrderType = "alphabetDesc")) 
# 2.9s

library(microbenchmark)
microbenchmark(glm.fit(y = iris$Sepal.Length,
                    x = model.matrix(~Sepal.Width + Species, data=iris), 
                    family = gaussian(identity)), unit="s") 
# 0.0004s, 7500x faster


# Model coefficients are returned in a similar format to R's native glm().
summary(model)

# maybe check out H20 instead? or tensorflow? or torch?
# https://docs.h2o.ai/sparkling-water/2.3/latest-stable/doc/rsparkling.html
# https://docs.h2o.ai/sparkling-water/2.3/latest-stable/doc/ml/sw_glm.html?highlight=linear%20regression
# https://docs.h2o.ai/h2o/latest-stable/h2o-docs/data-science/glm.html?highlight=linear%20regression






# some more tests

#create a sparkR DataFrame
DF <- createDataFrame(faithful)
head(DF)

# Create a simple local data.frame
localDF <- data.frame(name=c("John", "Smith", "Sarah"), age=c(19, 23, 18))

# Convert local data frame to a SparkR DataFrame
df <- createDataFrame(localDF)

# Print its schema
printSchema(df)
# root
#  |-- name: string (nullable = true)
#  |-- age: double (nullable = true)

# Create a DataFrame from a JSON file
path <- file.path(Sys.getenv("SPARK_HOME"), "examples/src/main/resources/people.json")
peopleDF <- read.json(path)
printSchema(peopleDF)

# Register this DataFrame as a table.
createOrReplaceTempView(peopleDF, "people")

# SQL statements can be run by using the sql methods provided by sqlContext
teenagers <- sql("SELECT name FROM people WHERE age >= 13 AND age <= 19")

# Call collect to get a local data.frame
teenagersLocalDF <- collect(teenagers)

# Print the teenagers in our dataset 
print(teenagersLocalDF)

# Stop the SparkContext now
sparkR.stop()