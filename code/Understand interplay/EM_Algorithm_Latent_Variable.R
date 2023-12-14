library(bnlearn)
setwd("/Users/barajascervantes.a/Documents/research/r_bayes_nets/")
data <- read.csv("df_work_binary.csv", header = TRUE, stringsAsFactors = TRUE)
df <- data[,c("Y","REVEL","SAI")]
latdf <- data.frame(df, LAT = factor(rep(NA, nrow(df)), levels = c("A", "B")))

# Randomly impute values
imputed = latdf 
imputed$LAT <- sample(factor(c("A", "B")), nrow(df), replace = TRUE)
table(imputed$LAT)


evidence.net <- empty.graph(names(latdf))
modelstring(evidence.net) = "[Y][LAT|Y][REVEL|LAT][SAI|LAT]"
graphviz.plot(evidence.net)
evidence.net

fitted <- bn.fit(evidence.net, imputed, method="bayes")
wl = data.frame(from=c("Y","LAT","LAT"), to=c("LAT","REVEL","SAI"))

r <- structural.em(latdf, fit = "bayes", impute="bayes-lw", start=fitted, maximize.args=list(whitelist = wl), return.all = TRUE)  

r

table(r$imputed$LAT)
