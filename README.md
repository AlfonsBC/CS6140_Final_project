# CS6140 Final project
The following repository contains code and final project report to the project carried out in the course CS6140: Machine Learning in the Fall of 2023. Primarily, this project uses Probabilistic Graphical Models to predict pathogenicity in Non-synonymous SNPs. This project have a `code` folder where focuses on two parts of this project:

a)  `/code/Understand_problem` : Understanding the Problem via Machine Learning

b) `/code/Understand_interplay` : Understand the effect of lines of evidence to predict pathogenicity of variants


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
---

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

Given the literature review such as [Tavtigian et. al 2018](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6336098/)  and [Pejaver et. al 2022](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC9748256/), the probabilistic bayesian approach to quantify the strength of evidence is given by the formula:

![equation](https://latex.codecogs.com/svg.image?&space;LR^&plus;(\mathbf{X}=x)=\frac{\text{posterior&space;odds&space;of&space;pathogenicity}}{\textrm{prior&space;odds&space;of&space;pathogenicity}}=\frac{\displaystyle\frac{P(Y=1|\mathbf{X}=x)}{1-P(Y=1|\mathbf{X}=x)}}{\displaystyle\frac{P(Y=1)}{1-P(Y=1)})


- Bayesian Networks: 
	- Fully Connected Network: The Graphical model trained with the bootstraped data with a fully-connected network was the following: ![FC_network](https://github.com/AlfonsBC/CS6140_Final_project/blob/master/images/FC_model.png?raw=true)
	
	- Naive Bayes approach: The Graphical model trained with a naive assumption that the two predictors are independent from each other is the following:
![NB_network](https://github.com/AlfonsBC/CS6140_Final_project/blob/master/images/NB_model.png?raw=true)

- Bayesian Network Latent Variable Model
	- With Latent Variable that can take two values $A,B$ using Expectation-Maximization algorithm that solves a higly complex non-linear optimization for the likelihood function:
![Latent_newtork](https://github.com/AlfonsBC/CS6140_Final_project/blob/master/images/BN_proposed.png?raw=true)


The results that we had obtained from the Bayesian Network perspective following the Bayesian approach from Tavtigian et al. 2018, gives us understanding in each of the pieces of evidence given the REVEL score and provide with scores SpliceAI to predict Pathogenicity, where in different regions the sum of their contributions are independent, but in other regions their sum of their points is not equal to the contributed by both. Meaning that this regions are dependent. In the following image, we plot the main results from our Bayesian Network perspective:
![Main_results](https://github.com/AlfonsBC/CS6140_Final_project/blob/master/images/results.png?raw=true)

## Discussion
---

The interesting reults that we get from deriving four machine learning algorithms including: Logistic Regression, Naive Bayes, Gradient Boosting Decision Tree and Multi-layer Perceptron is that in reality they provide of a good approach to predict benign variant. Nevertheless, the rare deseases caused for missense mutations are rare, less than $8\%$ of the population. That means, giving a framework to understand the regions to look at in regards scores to predict pathogenicity and then start building machine learning models to predict based on data is of vital importance. We noted that the best True Positive rate was the Probbilistic approach, which heavily relies on understanding the uncertainty via Probaiblity rules. In addition to have a naive hypothesis that says that the to score predictors are independent given the pathogenicity. But this is not what the data tells as we can see in the table above.

- **Significance of Results** with a User Case:

Suppose we have a variant which gives a REVEL score of $-4$ which will indicate that is a Benign with Very Strong be the strength of evidence. But that same variant, gives a SpliceAI score of $+2$ which means that the variant is highly disruptive. Using Tavtigian framework we are able to obtain the Individual points that a Clinician will use of $-4.9$ for REVEL being $Benign Very Strong (-4)$ and for SpliceAI score will use $2.4$ points fiven the $Pathogenic Moderate (+2)$. Clinician with the current guidelines to predict pathogenicity, will add the points which will be $-2.5$ which will be giving a result of Benign Moderate. But in reality, given the framework from literature review and Bayesian Networks, we are able to understand that that variant actually is a Variant of Uncertain Significance a VUS, as the joint score according to the Probabilities rules would be of $0.5=0$ and $0$ is Inconclusive.

What we can note, is that in the area of a well-known REVEL score from $-1$ to $-4$ that should be benign according to clinicians, but if there is a variant which has a REVEL score of $-1$ and SpliceAI score of $+2$ then that variant will be Pathogenic Very strong with $\infity$ be the maximum of point possible to reach which are $8$ points. The table will help clinicians to predict better using well-known predictors such as REVEL and SpliceAI. We notice that there are independence in the second and third row every three rows as we see in the table, where the predictors give SpliceAI score of $+1,+2$ primarily. 


## Conclusion
---

 Primarily, the results from the Fully-connected graphical mode are of vital importance to understand the effect that different lines of evidence within the code $BP4, PP3$ have, and understanding it will play a significant role to predict Pathogenicity in Varaints or Missense Mutations. The results of this project have shown that is necessarry to continue understanding the immmediate effect that the in-silico predictors have between one another to predict pathogenicity of a variant.


 Further work will be obtaining multiple predictors of Pathogenicity such as CADD, MutPred2, AlphaMissense, etc to predict Pathogenicity based on the Probabilistic Bayesian Framework of Tavtigian et al. There are a lot of experimantation to be done in the Latent/Hidden Variable Models to understand better the role of dependency that play REVEL and SpliceAI scores to pathogenicity, implementing code to Continuous/Discrete Bayesian Networks with Latent Variables would be an area to further explore as there is difficult to find one programming language that have those required functionalities.




## References
---
1. Tavtigian SV, Greenblatt MS, Harrison SM, Nussbaum RL, Prabhu SA, Boucher KM, Biesecker LG; ClinGen Sequence Variant Interpretation Working Group (ClinGen SVI). [Modeling the ACMG/AMP variant classification guidelines as a Bayesian classification framework.](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6336098/) Genet Med. 2018 Sep;20(9):1054-1060. doi: 10.1038/gim.2017.210. Epub 2018 Jan 4. PMID: 29300386; PMCID: PMC6336098. 
2. Pejaver V, Byrne AB, Feng BJ, Pagel KA, Mooney SD, Karchin R, O'Donnell-Luria A, Harrison SM, Tavtigian SV, Greenblatt MS, Biesecker LG, Radivojac P, Brenner SE; [ClinGen Sequence Variant Interpretation Working Group. Calibration of computational tools for missense variant pathogenicity classification and ClinGen recommendations for PP3/BP4 criteria.](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC9748256/) Am J Hum Genet. 2022 Dec 1;109(12):2163-2177. doi: 10.1016/j.ajhg.2022.10.013. Epub 2022 Nov 21. PMID: 36413997; PMCID: PMC9748256.
3. Nir Friedman. 1997. [Learning Belief Networks in the Presence of Missing Values and Hidden Variables.](https://dl.acm.org/doi/10.5555/645526.657145) In Proceedings of the Fourteenth International Conference on Machine Learning (ICML '97). Morgan Kaufmann Publishers Inc., San Francisco, CA, USA, 125–133.
4. Lauritzen, S.L., 1995. [The EM algorithm for graphical association models with missing data.](https://www.stats.ox.ac.uk/~steffen/papers/em95.pdf)Computational statistics & data analysis, 19(2), pp.191-201.
5. Packages:
- [bnlearn](https://www.bnlearn.com):  a packages for Bayesian Network learning and inference in R 
- [sklearn](https://scikit-learn.org/stable/): Machine Learning in Python
- [causalnex](https://causalnex.readthedocs.io/en/latest/01_introduction/01_introduction.html): Combining Machine Learning and Bayesian networks for causal reasoning.


