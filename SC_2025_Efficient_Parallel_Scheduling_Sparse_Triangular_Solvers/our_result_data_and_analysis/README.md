Our data files are processed through several Python Jupyter notebooks. The packages required for the specific notebooks are always imported and also explicitly specified at the beginning of the notebook. The first halves of the notebooks are essentially identical, just loading data from different folders. The end of the notebooks then produce the tables and figures from the tables.

We note that in the raw data in the notebooks, the data sets and algorithms have slightly different names than in the paper:
- for the data sets, florida (SuiteSparse), florida_metis (METIS), florida_Pchol (iChol), er (Erdos_Renyi), rb (Narrow Bandwidth)
- for the algorithms, HDAGG_BIN (HDagg), SpMP (SpMP), SMGreedyBspGrowLocalAutoCoresParallel (GrowLocal), SMFunOriGrowlv2 (Funnel+GL)


Data_Analysis_SC_paper.ipynb:
-----------------------------

Data_Analysis_SC_paper.ipynb processes most of our SpTRSV runs, and allows to reproduce most of the tables and figures in the paper.

Cells 38-42 process our main SpTRSV speed-ups for each of the 5 data sets. The "Geommean_serial" column provides the numbers in Table 7.1 for each data set. The "Geommean_supersteps_relative_to_Wavefront" column provides the numbers in Table 7.2 for each data set. The other two columns are only the relative speed-ups between the, i.e. the ratios of the values in the other two columns. We note that the numbers for Narrow Bandwidth are slightly different from those in the paper, since due to an error, a single instance was left out for the SpMP and HDagg runs; however, the difference is very small.

Cell 43 processes the amortization threshold on the florida (SparseSuite) data set, thus producing the numbers for Table 7.6 in the q25, median_profitability and q75 columns for each algorithm. Note that some of the values slightly differ from the paper, but this is only random noise: we ran these experiments again for Table 7.7 (see later), and used the same value in both tables to avoid confusion.

Cell 44 produces the violin plot that appears as Figure 1.2 in the paper.

Cell 46 produces a plot of the scheduling time that is currently included in the supplementary material.

Cell 50 produces the performance plot that appears as Figure 7.1 in the paper. We note that due to a small fix in the Funnel+GL algorithm, one of the curves slightly differs from the variant in the paper.


Data_Analysis_SC_paper_scheduling_thread.ipynb:
-----------------------------------------------

Data_Analysis_SC_paper_scheduling_thread.ipynb considers the block decomposition approach in our paper. The data here can be used to reproduce Table 7.7 in the paper. Cell 33 produces a table where the "median_profitability" column gives the data in "Amort. Threshold" column of Table 7.7. Cell 34 produces a table where the columns "GM_MulThr_sched_time_speedup" "GM_MulThr_FLOP_decrease" and "GM_MulThr_superstep_increase", respectively, give the values in the "Sched. Time", "Flops/s" and "Supersteps" columns of Table 7.7.


Data_Analysis_paper_permutation_clean.ipynb:
--------------------------------------------

Data_Analysis_paper_permutation_clean.ipynb considers the GrowLocal scheduler with and without the reordering technique, thus reproducing Table 7.3 in the paper. In particular, Cell 31 shows all the values for Table 7.3, in the appropriate order of the data sets. LOOP_PROCESSORS corresponds to the "Reordering" column, and NO_PERMUTE corresponds to the "No Reordering" column.


Data_Analysis_paper_arch_clean.ipynb:
-------------------------------------

Data_Analysis_paper_permutation_clean.ipynb is for Table 7.4 in the paper, the SpTRSV results on the SuiteSparse data set on different architectures. In Cell 30, the "Geomean" column of the table shows the speed-ups for each specific architecture and algorithm in Table 7.4.

Data_Analysis_paper_scaling_clean.ipynb:
-------------------------------------

Data_Analysis_paper_permutation_clean.ipynb shows the scaling of our algorithms: Table 7.5 and Figure 7.2 in the paper. Cell 34 shows the speed-ups of GrowLocal for different numbers of cores (note that Table 7.5 only shows some selected rows). Cell 35 generates the plot for Figure 7.2; this looks slightly different from the paper due to some last minute data fixes, but its essential characteristics are the same.
