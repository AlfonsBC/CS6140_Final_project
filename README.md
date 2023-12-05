# CS6140 Final project
The following repository contains code and final project report to the project carried out in the course CS6140: Machine Learning in the Fall of 2023 . Primarily, this project uses Probabilistic Graphical Models to predict pathogenicity in Non-synonymous SNPs.

## Introduction.
---
In 2015, the American College of Medical Genetics and Genomics (**ACMG**) and the Association for Molecular Pathology (**AMP**) published a guideline that provides a framework for sequence variant interpretation ([Richards et al. 2015](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4544753/)).  The guidelines defined $28$ criteria with criterion code that each was assigned the strength of evidence: _stand-alone (A), very strong (VS), strong (S), moderate (M) or supporting (P)_. We can note in the table the overview of the ACMG/AMP for classifying sequence variants organized by data type and strength ([Harrison et al. 2019](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4544753/)) there is a category "_Functional and Computational data_" where **BP4** and **PP3** provide of supporting evidence of Genign and Pathogenic criteria respectively. **PP4** Multiple lines of computational (_in silico_) tools have been produced to predict the pathogenicity of a variant, e.g. whether a variant will disrupt the function of a gene product such as RNA or protein. It can be based on conservation, evolutionary, splicing impact, etc. Whereas **BP4** is where many computational tools are used to predict the begnitity of a variant. 

Following the 2015 version of the ACMG/AMP guidelines ([Richards et al. 2015](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4544753/)) there is a caveat in that each computational algorithm should not be counted as an independent criterion. Given that fact for the code **PP3**, **BP4**, the current project tries to understand the interplay between the evidence of two in-silico tools: 
-  [SpliceAI](https://doi.org/10.1016/j.cell.2018.12.015) a deep-learning tool to identify splice variants.
-  [REVEL](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5065685/) an ensemble method for predicting the pathogenicity of missense variants on the basis of individual in-silico tools: MutPred, FATHMM, VEST, PolyPhen, SIFT, PROVEAN, MutationAssessor, MutationTaster, LRT, GERP, SiPhy, phyloP, and phastCons.
  
Then, measure to what extent their predictions can be considered independent pieces of evidence or not, using a Bayesian approach with positive likelihood ratio and with a Probabilistic Graphical Model approach, in particular Bayesian Belief Networks. Moreover, this project aims to understand the role it plays the pathogenicity and the in-silico tools predictions in terms of a Latent Variable Model. According to the compact representation for the set of _conditional independencies_ in Bayesian Networks denoted by $\mathcal{I}_{\ell}(\mathcal{G})$, where for each variable $X_i$

$$(X_i \perp NonDescendants_{X_i} | Pa^{\mathcal{G}}_{X_i})$$

Thus, adding a hidden variable acording to expert knowledge will make the contributions from each in-silico predictors independent given the hidden variable, but they will be dependent from both the Latent Variable and the Pathogenicity $Y$ variable.

