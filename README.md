# CS6140 Final project
The following repository contains code and final project report to the project carried out in the course CS6140: Machine Learning in the Fall of 2023 . Primarily, this project uses Probabilistic Graphical Models to predict pathogenicity in Non-synonymous SNPs.

## Introduction.
---
In 2015, the American College of Medical Genetics and Genomics (**ACMG**) and the Association for Molecular Pathology (**AMP**) published a guideline that provides a framework for sequence variant interpretation ([Richards et al. 2015](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4544753/)).  The guidelines defined $28$ criteria with criterion code that each was assigned the strength of evidence: _stand-alone (A), very strong (VS), strong (S), moderate (M) or supporting (P)_. We can note in the table the overview of the ACMG/AMP for classifying sequence variants organized by data type and strength ([Harrison et al. 2019](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4544753/)) there is a category "_Functional and Computational data_" where **BP4** and **PP3** provide of supporting evidence of Genign and Pathogenic criteria respectively. **PP4** Multiple lines of computational (_in silico_) tools have been produced to predict the pathogenicity of a variant, e.g. whether a variant will disrupt the function of a gene product such as RNA or protein. It can be based on conservation, evolutionary, splicing impact, etc. Whereas **BP4** is where many computational tools are used to predict the begnitity of a variant. 

Following the 2015 version of the ACMG/AMP guidelines ([Richards et al. 2015](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4544753/)) there is a caveat in that each computational algorithm should not be counted as an independent criterion. Given that fact for the code **PP3**, **BP4**, the current project tries to understand the interplay between the evidence of two in-silico tools: 
-  [SpliceAI](https://doi.org/10.1016/j.cell.2018.12.015) a deep-learning tool to identify splice variants.
-  [REVEL](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5065685/) an ensemble method for predicting the pathogenicity of missense variants on the basis of individual in-silico tools: MutPred, FATHMM, VEST, PolyPhen, SIFT, PROVEAN, MutationAssessor, MutationTaster, LRT, GERP, SiPhy, phyloP, and phastCons.

There is a myriad of in-silico tools which design is motivated by the discovery of  novel variants ad hypothesis generation for experimental purposes. The purpose of this project is to consider the two well known computational tools to predict pathogenicity on variants depending on the _splice effect_ SpliceAI and REVEL which is a _mix of protein and splice effect_ on missene variants and meaasure the dependence across strengths of evidence using a probabilistic approach.

In short, the motivation of this project is to find the dependence between two lines of evidence and test the hypothesis that two or more lines of evidence in the code **PP3, BP4** cannot be considered to be the sum of individuals contributions.

## Experimentation.
---
### Modelling
The aim of this research project is to measure to what extent their predictions can be considered independent pieces of evidence or not, using a Bayesian approach with positive likelihood ratio. We will model the $LR^{+}$ using a two different approaches but in essence relaying on Prabiblistic principles. A Probabilistic Graphical Model approach, in particular Bayesian Belief Networks and the second is a Bayes Rule Analysis, which basically is applying the Bayes Rule. Moreover, this project aims to understand the role it plays to the pathogenicity and the in-silico tools predictions in terms of a Latent Variable Model. According to the compact representation for the set of _conditional independencies_ in Bayesian Networks denoted by $\mathcal{I}_{\ell}(\mathcal{G})$, where for each variable $X_i$

$$(X_i \perp NonDescendants_{X_i} | Pa^{\mathcal{G}}_{X_i})$$

Thus, adding a hidden variable acording to expert knowledge will make the contributions from each in-silico predictors independent given the hidden variable, but they will be dependent from both the Latent Variable and the Pathogenicity $Y$ variable. This approach we are going to use using a Latent Variable Model via _Expectation Maximization Algorithm_. The parameters to be learned are the _Conditional Probability Tables_ for each random variable. We are learning a model wth four discrete random variables $Y, L, SAI, REVEL$ as we see in the image bellow.

### Data set
--

To perform this project, we needed to obtain a data set to predict pathogenicity according to a strict protocol that clinicians and researchers had agreed to which are the following:

1. ClinVar variants from 2020 and 2023. There was some processing to be done to the data to have nearly $970,000$ variants in nearly $18,000$ genes. 
2. Then, we keep non-VUS variants with a star rating greater than 1. The number of variants are $81,000$ in $10,661$ genes
3. Then we remove variants in genes containning only benign variants. Because there were genes that had great density of benign variant rather than pathogenic. After this process we got around $44,629$ variants in $1,629$ genes
4. Then we remove variants with an _Allele Frequency_ greater or equal than $0.01$. For this, we use the gnomAD data set which provides us of those variants that are common in the population. Thus, we will be obtaining the variants that are rare in the population. That means we got nearly $43,361$ variants for $1,629$ genes
5. Obtain REVEL scores which is an Ensemble method that outputs the prediction combining multiple predictors.  This reduces to have nearly $42,572$ variants in $1,603$ genes
6. Finally, we obtain the probabili of Splicing using the Neural Network tool via SpliceAI tool. This in the end gives us close to $42,013$ variants and $1,584$ genes

In the end, we obtain nearly $21,311$ _Pathogenic/Likely pathogenic variants_ and $20,702$ _Benign/Likely benign_ variants. This leads to have a prevalence of pathogenicty close to half percent. Nevertheless, this is not the correct to calculate the prior, that is why we followed closely the literature review with Pejaver et. al. that calculated the prevalence of pathogenicty to be $0.0441$.


## Results
---

### Understanding the Problem via Machine Learning

To understand the importance of this project, we first provided of code in the folder `/code/Understand_problem` that builds 4 different Machine Learning Models with data that aligns to reality, where for every missense mutation that is pathogenic or disrups the function of the protein, there are nearly $22.3$ benign variants. We model the binary classification problem with the following 2 classical approaches, one ensemble method and a neural network. 

- Logistic Regression
- Naive Bayes 
- Gradient Boosting Decision Tree
- Multi Layer Perceptron

We fitted the models with bootstrap data according to reality, to understand the importance of this problem in the area of Machine Learning. First, we use nearly $496,831$ number of boostrap variants to the before mentioned models, and obtained the following results using a standard approach of train/test split of $70\%$ and $30\%$ :

- The model that had the worst accuracy was the Naive Bayes approach where it has the naive hypothesis that the two scores Revel and SpliceAI are independent of each other. The rest three models were marginally close in accuracy in both train/test set nearly to $97\%$.

![Train/Test Accuracy Model Comparison](https://github.com/AlfonsBC/CS6140_Final_project/blob/master/images/accuracy_comparison.png?raw=true)

- Nevertheless, when comparing their respectives confusion matrices, Naive Bayes model approach was the one that obtained the greatest True Positive in the test set with $54\%$, the rest of the models were close to the $40\%$ True Positive rates. The following image is the confusion matrix for Naive Bayes model in the Test set.

![Naive Bayes Conf Matrix](https://github.com/AlfonsBC/CS6140_Final_project/blob/master/images/conf_matrix_naive.png?raw=true)

We can add more in the discussion, but a quick intuition that this is one of the best that actually predicts what we want to predict, e.g. Pathogenic variants, is that is based in a Probabilistic Framework. Thus, in the following part of the project, we will discuss the results obtained in the Bayesian Networks following a Probabilistic Graphical Model.

### Understanding the interplay of lines of evidence

- Bayes Rule Analysis
- Bayesian Networks:
	- Fully Connected Network
	- Naive Bayes approach
- Bayesian Network Latent Variable Model
	- With Latent Variable that can take two values $A,B$
	- With Latent Variable that can take three values $A,B,C$


## Discussion
---


## Conclusion
---


## References
---


