This directory contains the different algorithms applied in the paper.

The "our_schedulers" subfolder contains all of our scheduling algorithms, as well as the implementations for the Cilk, ETF and BL-EST baselines. More details in the README file within.

The separate "hdagg_baseline" subfolder contains the code necessary to run the external HDagg baseline by Zarebavani et al. This is essentially a pruned version of the same scripts as in the "our_schedulers" folder, only keeping the files and classes necessary for reading input files, invoking the compiled version of HDagg, and evaluating the cost of the returned partition. More details in the README file within.

