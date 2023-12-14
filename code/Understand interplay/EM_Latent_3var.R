library(bnlearn)
setwd("/Users/barajascervantes.a/Documents/research/r_bayes_nets/")
data <- read.csv("df_work_binary.csv", header = TRUE, stringsAsFactors = TRUE)
df <- data[,c("Y","REVEL","SAI")]
latdf <- data.frame(df, LAT = factor(rep(NA, nrow(df)), levels = c("A", "B","C")))

# Randomly impute values
imputed = latdf 
imputed$LAT <- sample(factor(c("A", "B", "C")), nrow(df), replace = TRUE)

evidence.net <- empty.graph(c("Y","REVEL", "SpliceAI"))
modelstring(evidence.net) = "[Y][REVEL|Y][SpliceAI|Y]"
graphviz.plot(evidence.net)
evidence.net

fitted <- bn.fit(evidence.net, imputed, method="mle")
wl = data.frame(from=c("Y","LAT","LAT"), to=c("LAT","REVEL","SAI"))

r <- structural.em(latdf, fit = "bayes", impute="bayes-lw", start=fitted, maximize.args=list(whitelist = wl), return.all = TRUE)  

r

table(r$imputed$LAT)
