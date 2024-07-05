# Nullspace-code

This repo contains the implementation of randomized small-block Lanczos for null space computation proposed in [1] and scripts to reproduce the numerical experiments.

The main function of randomized small-block Lanczos method is TRlanczos.m in the num_exp folder.

To reproduce Table 1 in Section 1, run main_intro.m in the exp_intro folder.

To reproduce tables and figures in Section 5, uncomment and run corresponding codes at main.m in the num_exp folder.

[1]: Daniel Kressner and Nian Shao. ``A randomized small-block Lanczos method for large-scale null space computations.''
