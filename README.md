# CS6140 Final project
The following repository contains code and final project report to the project carried out in the course CS6140: Machine Learning in the Fall of 2023 . Primarily, this project uses Probabilistic Graphical Models to predict pathogenicity in Non-synonymous SNPs.

## Introduction.
---
In 2015, the American College of Medical Genetics and Genomics (**ACMG**) and the Association for Molecular Pathology (**AMP**) published a guideline that provides a framework for sequence variant interpretation ([Richards et al. 2015](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4544753/)).  The guidelines defined $28$ criteria with criterion code that each was assigned the strength of evidence: _stand-alone (A), very strong (VS), strong (S), moderate (M) or supporting (P)_. We can note in the table the overview of the ACMG/AMP for classifying sequence variants organized by data type and strength ([Harrison et al. 2019](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4544753/)) there is a category "_Functional and Computational data_" where **BP4** and **PP3** provide of supporting evidence of Genign and Pathogenic criteria respectively. **PP4** Multiple lines of computational (_in silico_) tools have been produced to predict the pathogenicity of a variant, e.g. whether a variant will disrupt the function of a gene product such as RNA or protein. It can be based on conservation, evolutionary, splicing impact, etc. Whereas **BP4** is where many computational tools are used to predict the begnitity of a variant. 

Following the 2015 version of the ACMG/AMP guidelines ([Richards et al. 2015](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4544753/)) there is a caveat in that each computational algorithm should not be counted as an independent criterion. Given that fact for the code **PP3**, **BP4**, the current project tries to understand the interplay between the evidence of two in-silico tools: 
-  [SpliceAI](https://doi.org/10.1016/j.cell.2018.12.015) a deep-learning tool to identify splice variants.
-  [REVEL](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5065685/) an ensemble method for predicting the pathogenicity of missense variants on the basis of individual in-silico tools: MutPred, FATHMM, VEST, PolyPhen, SIFT, PROVEAN, MutationAssessor, MutationTaster, LRT, GERP, SiPhy, phyloP, and phastCons.

There is a myriad of in-silico tools which design is motivated by the discovery of  novel variants ad hypothesis generation for experimental purposes. The purpose of this project is to consider the two well known computational tools to predict pathogenicity on variants depending on the _splice effect_ SpliceAI and REVEL which is a _mix of protein and splice effect_ on missene variants and meaasure the dependence across strengths of evidence using a probabilistic approach.

In short, the motivation of this project is to find the dependence between two lines of evidence and test the hypothesis that two or more lines of evidence in the code **PP3, BP$** cannot be considered to be the sum of individuals contributions.

## Setup.
---
### Modelling
The aim of this research project is to measure to what extent their predictions can be considered independent pieces of evidence or not, using a Bayesian approach with positive likelihood ratio. We will model the $LR^{+}$ using a two different approaches but in essence relaying on Prabiblistic principles. A Probabilistic Graphical Model approach, in particular Bayesian Belief Networks and the second is a Bayes Rule Analysis, which basically is applying the Bayes Rule. Moreover, this project aims to understand the role it plays to the pathogenicity and the in-silico tools predictions in terms of a Latent Variable Model. According to the compact representation for the set of _conditional independencies_ in Bayesian Networks denoted by $\mathcal{I}_{\ell}(\mathcal{G})$, where for each variable $X_i$

$$(X_i \perp NonDescendants_{X_i} | Pa^{\mathcal{G}}_{X_i})$$

Thus, adding a hidden variable acording to expert knowledge will make the contributions from each in-silico predictors independent given the hidden variable, but they will be dependent from both the Latent Variable and the Pathogenicity $Y$ variable. This approach we are going to use using a Latent Variable Model via _Expectation Maximization Algorithm_. The parameters to be learned are the _Conditional Probability Tables_ for each random variable. We are learning a model wth four discrete random variables $Y, L, SAI, REVEL$ as we see in the image bellow.



