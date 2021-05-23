# Mixture_models_simulation

Codes for my Master's final project.

$Abstract
Finite mixture models provide a convenient framework for model-based clustering. Traditionally, the model parameters are estimated by maximum likelihood estimation, fulfilled by the expectation-maximization (EM) algorithm. Such approach to clustering has many advantages but also several pitfalls. Some of those issues can be overcome by varying the EM algorithm. We describe two variants of the EM algorithm, namely the Classification EM (CEM) and the Stochastic EM (SEM). 

We study the performance of the standard EM, CEM, and SEM measured by the Adjusted Rand index in simulation studies for two different mixtures. First, we examine a finite Gaussian mixture model which is, by far, the most popular and widely studied mixture model. 
%Based on our study, the three procedures suffer more from higher overlap among clusters rather than from increasing the number of dimensions $p$, except for CEM which appeared to be more impacted by higher $p$. 
%The results obtained by the standard EM were the most accurate while CEM and SEM were faster especially for the first couple of steps. 
%When it comes to accuracy, SEM outperformed CEM in every simulation scenario. 
Then, we present our results for a finite mixture of Markov chains. We conducted a simulation study similar to the one with Gaussian mixture but additionally, we studied how frequently the three procedures identify the correct number of components $K$. 
%Based on the results, the accuracy was very similar for all the three procedures, yet the standard EM appeared to be identifying the correct $K$ earlier (for lower sample sizes and shorter sequences) and resulted in slightly more accurate clusterings. As for Gaussian mixtures, running CEM or SEM for the first couple of steps and then running the standard EM seems advisable to achieve faster progress while maintaining high accuracy.

**Keywords**: Finite mixture models, clustering, EM algorithm, classification EM, stochastic EM, adjusted Rand Index.
