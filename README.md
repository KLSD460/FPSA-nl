# FPSA-nl

Some code for experiments done in the paper "A single-loop Proximal Subgradient Algorithm for A Class Structured Fractional Programs".

References:
Deren Han, Min Tao and Zihao Xia.
A Single-loop Proximal Subgradient Algorithm for A Class Structured Fractional Programs,
https://arxiv.org/abs/2503.12176

# FPSA-nl_L1oSK

This folder include the MATLAB code of FPSA-nl for solving $L_1/S_{\kappa}$ sparse signal recovery problem in section 5.1.

- testL1oSK.m: the test file for $L_1/S_{\kappa}$ sparse signal recovery problem.

- FPSA_nl_L1oSK.m: the code for FPSA-nl Algorithm.

- data_generator_A.m: the file for generating matrix $A$.

- norm_maxK.m: compute the max K-norm for $x$.

- randsample_separated.m: random sampling K integers from 1-N with spacing at least L.

- StatRes_L1oLK.m: compute the StatRes for L1oLk problem.
