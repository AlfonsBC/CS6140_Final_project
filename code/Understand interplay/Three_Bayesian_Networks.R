library(bnlearn)
setwd("/Users/barajascervantes.a/Documents/research/r_bayes_nets/")

# Modelling the Dataset obtained from ClinVar
data <- read.csv("finalboot.csv", header = TRUE, stringsAsFactors = TRUE)
df <- data[,c("Y","REVEL","SAI")]
nrow(df)
colnames(df) <- c("Y", "R", "S")
grant_graph <- empty.graph(names(df))
fc_graph <- empty.graph(names(df))


modelstring(grant_graph) = "[Y][S|Y][R|Y]"
modelstring(fc_graph) = "[Y][S|Y][R|Y][R|S]"
df$Y <- replace(df$Y, df$Y == 0, "B")
df$Y <- replace(df$Y, df$Y == 1, "P")
df$Y <- as.factor(df$Y)

fitem_grant <- bn.fit(grant_graph, df, method = "mle")
fitem_fc <- bn.fit(fc_graph, df, method = "mle")
graphviz.plot(fitem_fc)
graphviz.plot(fitem_grant)


fitem_grant$Y
fitem_grant$R
fitem_grant$S
#####################
# Counting each cell of BN
sum(df$R == "-4" & df$S == "I")
levels_sai <- c("I", "+1", "+2")
levels_rev <- c("-4", "-3", "-2","-1", "I", "+1", "+2", "+3","+4")

mat_cell <- matrix(0, 27, 2)
row_idx <- 1
for (j in levels_rev){
  for (i in levels_sai){
    mat_cell[row_idx, 1] = sum(df$R == j & df$S == i)   
    mat_cell[row_idx, 2] = round((sum(df$R == j & df$S == i) / 490153),5)
    row_idx <- row_idx + 1
  }
}
mat_cell

##############################################################################################################################
# Naive Bayes
library(gRain)
# BN into a crafted tree 
junction <- compile(as.grain(fitem_grant))
# Revel has 8 levels 
# Splice AI has 3 levels
# 24 LR+ in total if jointly and  LR
## R = BI and S = B 
# P(Y=1 | R="-4" , S="I") 
jrevsai <- setEvidence(junction, nodes = c("R","S"), states = c("-4", "I"))
posterior_path <- querygrain(jrevsai, nodes = "Y",type="conditional")
alpha_prior <- fitem_grant$Y$prob[2]
LR_plus_joint <- unname((posterior_path /(1 - posterior_path)) / (alpha_prior / (1 - alpha_prior)))

# P(Y=1 | R="")
jrev <- setEvidence(junction, nodes = "R", states = "BI")
post_path_rev <- querygrain(jrev, nodes = "Y")$Y[2][1]
LR_plus_rev <- unname((post_path_rev /(1 - post_path_rev)) / (alpha_prior / (1 - alpha_prior)))


# P(Y=1 | S="")
jsai <- setEvidence(junction, nodes = "S", states = "B")
post_path_sai <- querygrain(jsai, nodes = "Y")$Y[2][1]
LR_plus_sai <- unname((post_path_sai /(1 - post_path_sai)) / (alpha_prior / (1 - alpha_prior)))

# Comparing the product
LR_plus_sai*LR_plus_rev
# with the joint LR+
LR_plus_joint


get_LR_revel_sai <- function(R_state, S_state){
  c <- log(1124)
  # P(Y=1 | R=" , S="") 
  jrevsai <- setEvidence(junction, nodes = c("R","S"), states = c(R_state, S_state))
  posterior_path <- querygrain(jrevsai, nodes = "Y",type="conditional")[2]
  alpha_prior <- fitem_grant$Y$prob[2]
  LR_plus_joint <- unname( (posterior_path /(1 - posterior_path)) / (alpha_prior / (1 - alpha_prior)) )
  ptsjoint <- round( (8 * log(LR_plus_joint) / c) , 1)
  # P(Y=1 | R="")
  jrev <- setEvidence(junction, nodes = "R", states = R_state)
  post_path_rev <- querygrain(jrev, nodes = "Y",type="conditional")[2]
  LR_plus_rev <- unname((post_path_rev /(1 - post_path_rev)) / (alpha_prior / (1 - alpha_prior)))
  ptsrev <- round( (8*log(LR_plus_rev) / c)  ,1)
  
  # P(Y=1 | S="")
  jsai <- setEvidence(junction, nodes = "S", states = S_state)
  post_path_sai <- querygrain(jsai, nodes = "Y",type="conditional")[2]
  LR_plus_sai <- unname((post_path_sai /(1 - post_path_sai)) / (alpha_prior / (1 - alpha_prior)))
  ptssai <- round( (8*log(LR_plus_sai) / c)  ,1)
  # Comparing the product
  ptssum <- round(ptsrev + ptssai, 1)
  numbers <- c(ptsrev, ptssai, ptssum, max(ptsrev,ptssai), ptsjoint )
  return(numbers)
}
levels_sai <- c("I", "+1", "+2")
levels_rev <- c("-4", "-3", "-2","-1", "I", "+1", "+2", "+3","+4")

mat <- matrix(0, 27, 7)
row_idx <- 1
for (j in levels_rev){
  for (i in levels_sai){
    mat[row_idx , 1] = j
    mat[row_idx , 2] = i
    mat[row_idx , 3:7] = get_LR_revel_sai(j, i)
    row_idx <- row_idx + 1
  }
}
colnames(mat) <- c("REV", "SAI","PtsREVEL", "PtsSAI", "SUM", "MAX" ,"JOINT")

strtop <-"""
\begin{table}[h!]
\centering
\begin{tabular}{||c c | c c || c c c ||} 
 \hline
"""
# -4& 0& 456456 & 456 & 456& 456 & 456 & 456 & 456 \\   
names <- c("\textbf{REV", "SAI","PtsREVEL", "PtsSAI", "SUM", "MAX" ,"JOINT}")
names_str <- paste(names, collapse= "}&\textbf{")
final <- paste(names_str, "\\")

# Latex Print Table
row_idx <- 1
strfinal <- ""
for (i in 1:27){
  vectora <- as.vector(mat[i,1:7])
  strvector <- paste(vectora, collapse="&")
  final <- paste(strvector, "\\")
  strfinal <- paste(strfinal, final, collapse = "")
  row_idx <- row_idx + 1
}

strfinal

plot(df$S)
round(0.9853567456, 2)
################################################################################################################
### Fully-Connected Graph
junction <- compile(as.grain(fitem_fc))
get_LR_revel_sai <- function(R_state, S_state){
  c <- log(1124)
  # P(Y=1 | R=" , S="") 
  jrevsai <- setEvidence(junction, nodes = c("R","S"), states = c(R_state, S_state))
  posterior_path <- querygrain(jrevsai, nodes = "Y", type="conditional")[2]
  alpha_prior <- fitem_grant$Y$prob[2]
  LR_plus_joint <- unname( (posterior_path /(1 - posterior_path)) / (alpha_prior / (1 - alpha_prior)) )
  ptsjoint <- round( (8 * log(LR_plus_joint) / c) , 1)
  # P(Y=1 | R="")
  jrev <- setEvidence(junction, nodes = "R", states = R_state)
  post_path_rev <- querygrain(jrev, nodes = "Y", type="conditional")[2]
  LR_plus_rev <- unname((post_path_rev /(1 - post_path_rev)) / (alpha_prior / (1 - alpha_prior)))
  ptsrev <- round( (8*log(LR_plus_rev) / c)  ,1)
  
  # P(Y=1 | S="")
  jsai <- setEvidence(junction, nodes = "S", states = S_state)
  post_path_sai <- querygrain(jsai, nodes = "Y", type="conditional")[2]
  LR_plus_sai <- unname((post_path_sai /(1 - post_path_sai)) / (alpha_prior / (1 - alpha_prior)))
  ptssai <- round( (8*log(LR_plus_sai) / c)  ,1)
  # Comparing the product
  ptssum <- round(ptsrev + ptssai, 1)
  numbers <- c(ptsrev, ptssai, ptssum, max(ptsrev,ptssai), ptsjoint )
  return(numbers)
}
levels_sai <- c("I", "+1", "+2")
levels_rev <- c("-4", "-3", "-2","-1", "I", "+1", "+2", "+3","+4")

mat <- matrix(0, 27, 7)
row_idx <- 1
for (j in levels_rev){
  for (i in levels_sai){
    mat[row_idx , 1] = j
    mat[row_idx , 2] = i
    mat[row_idx , 3:7] = get_LR_revel_sai(j, i)
    row_idx <- row_idx + 1
  }
}
colnames(mat) <- c("REV", "SAI","PtsREVEL", "PtsSAI", "SUM", "MAX" ,"JOINT")

strtop <-"""
\begin{table}[h!]
\centering
\begin{tabular}{||c c | c c || c c c ||} 
 \hline
"""
# -4& 0& 456456 & 456 & 456& 456 & 456 & 456 & 456 \\   
names <- c("\textbf{REV", "SAI","PtsREVEL", "PtsSAI", "SUM", "MAX" ,"JOINT}")
names_str <- paste(names, collapse= "}&\textbf{")
final <- paste(names_str, "\\")

# Latex Print Table
row_idx <- 1
strfinal <- ""
for (i in 1:27){
  vectora <- as.vector(mat[i,1:7])
  strvector <- paste(vectora, collapse="&")
  final <- paste(strvector, "&\\")
  strfinal <- paste(strfinal, final, collapse = "")
  row_idx <- row_idx + 1
}

strfinal





##############################################################################################################################
### Latent Variable in Python adding the CPTs obtained
#### Latent --> 2 values
# Variables and their levels
Y.lv <- c("B", "P")
L.lv <- c(0, 1)
R.lv <- c("+1", "+2", "+3","+4", "-1", "-2", "-3","-4", "I")
S.lv <-  c("+1", "+2", "I")


Y.prob <- array(c(0.957106,0.042894), dim=2, dimnames = list(Y=Y.lv))
Y.prob
L.prob <- array(c(0.066065,0.933935, 0.999999, 0.000001), dim = c(2,2), dimnames = list(L = L.lv, Y = Y.lv))
L.prob
S.prob <- array(c(0.016114, 0.002346, 0.981540, 0.002527, 0.000600, 0.996873 ), dim = c(3,2), dimnames = list(S= S.lv, L = L.lv))
S.prob
R.prob <- array(c(0.158307, 0.231500, 0.141220, 0.147479,0.041892,0.030209, 0.013849, 0.032354, 0.203190, 1.556179e-02, 1.619213e-04, 1.420845e-09,  2.229892e-18, 1.402624e-01, 2.793804e-01, 2.717417e-01, 1.311810e-01, 1.617109e-01  ), dim = c(9,2), dimnames = list(R = R.lv, L = L.lv))
R.prob
dag2 <- model2network("[Y][L|Y][S|L][R|L]")
dag2

cpt <- list(Y = Y.prob, L = L.prob,  S = S.prob , R = R.prob)
bn <- custom.fit(dag2,  cpt)
bn

graphviz.plot(bn)

junction <- compile(as.grain(bn))

levels_sai <- c("I", "+1", "+2")
levels_rev <- c("-4", "-3", "-2","-1", "I", "+1", "+2", "+3","+4")

get_LR_revel_sai <- function(R_state, S_state){
  c <- log(1124)
  # P(Y=1 | R=" , S="") 
  jrevsai <- setEvidence(junction, nodes = c("R","S"), states = c(R_state, S_state))
  posterior_path <- querygrain(jrevsai, nodes = "Y", type="conditional")[2]
  alpha_prior <- fitem_grant$Y$prob[2]
  LR_plus_joint <- unname( (posterior_path /(1 - posterior_path)) / (alpha_prior / (1 - alpha_prior)) )
  ptsjoint <- round( (8 * log(LR_plus_joint) / c) , 1)
  # P(Y=1 | R="")
  jrev <- setEvidence(junction, nodes = "R", states = R_state)
  post_path_rev <- querygrain(jrev, nodes = "Y", type="conditional")[2]
  LR_plus_rev <- unname((post_path_rev /(1 - post_path_rev)) / (alpha_prior / (1 - alpha_prior)))
  ptsrev <- round( (8*log(LR_plus_rev) / c)  ,1)
  
  # P(Y=1 | S="")
  jsai <- setEvidence(junction, nodes = "S", states = S_state)
  post_path_sai <- querygrain(jsai, nodes = "Y", type="conditional")[2]
  LR_plus_sai <- unname((post_path_sai /(1 - post_path_sai)) / (alpha_prior / (1 - alpha_prior)))
  ptssai <- round( (8*log(LR_plus_sai) / c)  ,1)
  # Comparing the product
  ptssum <- round(ptsrev + ptssai, 1)
  numbers <- c(ptsrev, ptssai, ptssum, max(ptsrev,ptssai), ptsjoint )
  return(numbers)
}


mat2 <- matrix(0, 27, 7)
row_idx <- 1
for (j in levels_rev){
  for (i in levels_sai){
    mat2[row_idx , 1] = j
    mat2[row_idx , 2] = i
    mat2[row_idx , 3:7] = get_LR_revel_sai(j, i)
    row_idx <- row_idx + 1
  }
}
colnames(mat2) <- c("REV", "SAI","PtsREVEL", "PtsSAI", "SUM", "MAX" ,"JOINT")

mat2
# Latex Print Table
row_idx <- 1
strfinal <- ""
for (i in 1:27){
  vectora <- as.vector(mat2[i,1:7])
  strvector <- paste(vectora, collapse="&")
  final <- paste(strvector, "&\\")
  strfinal <- paste(strfinal, final, collapse = "")
  row_idx <- row_idx + 1
}

strfinal


### Latent Variable in Python adding the CPTs obtained
#### Latent --> 3 values
# Variables and their levels
Y.lv <- c("B", "P")
L.lv <- c(0, 1)
R.lv <- c("+1", "+2", "+3","+4", "-1", "-2", "-3","-4", "I")
S.lv <-  c("+1", "+2", "I")


Y.prob <- array(c(0.957106,0.042894), dim=2, dimnames = list(Y=Y.lv))
Y.prob
L.prob <- array(c(0.024404,0.975596, 9.999995e-01,5.062106e-07), dim = c(2,2), dimnames = list(L = L.lv, Y = Y.lv))
L.prob
S.prob <- array(c(0.020235, 0.003326, 0.976439, 0.002656, 0.000581, 0.996763 ), dim = c(3,2), dimnames = list(S= S.lv, L = L.lv))
S.prob
R.prob <- array(c(0.126825, 0.188583, 0.228961, 0.277469, 0.017278, 0.009967, 0.000853, 0.000047, 0.150017, 4.497046e-02, 1.815782e-02, 3.786084e-04, 1.116209e-13, 2.094805e-01, 2.908890e-01 , 8.307532e-02, 1.489809e-02 , 3.381502e-01 ), dim = c(9,2), dimnames = list(R = R.lv, L = L.lv))
R.prob
dag3 <- model2network("[Y][L|Y][S|L][R|L]")
dag3

cpt <- list(Y = Y.prob, L = L.prob,  S = S.prob , R = R.prob)
bn2 <- custom.fit(dag3,  cpt)
bn2

graphviz.plot(bn2)

junction <- compile(as.grain(bn2))

levels_sai <- c("I", "+1", "+2")
levels_rev <- c("-4", "-3", "-2","-1", "I", "+1", "+2", "+3","+4")

mat2 <- matrix(0, 27, 7)
row_idx <- 1
for (j in levels_rev){
  for (i in levels_sai){
    mat2[row_idx , 1] = j
    mat2[row_idx , 2] = i
    mat2[row_idx , 3:7] = get_LR_revel_sai(j, i)
    row_idx <- row_idx + 1
  }
}
colnames(mat2) <- c("REV", "SAI","PtsREVEL", "PtsSAI", "SUM", "MAX" ,"JOINT")

mat2
# Latex Print Table
row_idx <- 1
strfinal <- ""
for (i in 1:27){
  vectora <- as.vector(mat2[i,1:7])
  strvector <- paste(vectora, collapse="&")
  final <- paste(strvector, "&\\")
  strfinal <- paste(strfinal, final, collapse = "")
  row_idx <- row_idx + 1
}

strfinal

### Naive Bayes
## grant_graph, fitem_grant
score(grant_graph, data=df, type="bic")
### Fully-connected
## fg_graph, fitem_fc
score(fc_graph, data=df, type="bic")
### Lat with two var
## dag2
### Lat with three var
## dag3


