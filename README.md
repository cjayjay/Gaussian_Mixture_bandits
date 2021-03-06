# Gaussian_Mixture_bandits
In multi-armed bandits problems, one central goal is to utilize structures in arms to 1) solve difficult problems or 2) further minimize regrets. Most work is focused on using known information. Lipshitz bandits assumed lipshitiz reward functions and discretized the arms space to solve uncountably many arms problem. Spectral bandits uses a given similarity graph to represent the pair-wise similarity information on arms, and assumed a smooth function on this graph to extract more information every time an arm is played. Gaussian bandits assumed a correlated gaussian prior on the expected rewards. Therefore, the rewards from one arm will give information about rewards on all arms.



There are some works that learns the structure between arms in an online manner. Gentile et al proposed to dynamically cluster the arms on a graph as samples are collected. 

In our project, we propose a two level mixture model which is capable of utilizing known information and learning certain structure between arms. We assume arms with similar rewards can be clustered together, the expected rewards on each cluster is generated from one common distribution. In our algorithm, we first select a promising cluster, then select an arm within that cluster.  We show that this will give better upper bound on the bayes regret, and performs better on simulated experiments compared to popular bandit algorithms.
