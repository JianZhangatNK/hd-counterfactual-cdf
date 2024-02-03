# hd-counterfactual-cdf
Refer to our working paper "Jun Cai, Jian Zhang, Yahong Zhou(2024). Estimation and Inference of Counterfactual Cumulative Distribution Function in a High-dimension Framework: A Distributional Oaxaca–Blinder Decomposition Application". Jun Cai(Huazhong University of Science and Technology), Jian Zhang(Nankai University), Yahong Zhou(Shanghai University of Finance and Economics). 
In this paper, we propose two approaches to estimate counterfactual CDF, namely, "DML" and "PSDB". Our approaches attain semiparametric efficiency bounds, have double/de-bias properties and can be used to estimate high dimensional counterfactual CDF.
File "counterCDF.R" contains functions for estimating counterfactual CDF and multiplier bootstraping samples.
File "counterCDF_stats.R" contains functions for calculating average structrural effects(ASE), average treatment effects on the treated(ATT), quantile treatment effects on the treated(QTT), and p-values of Kolmogorov-Smirnov test.
File "counterCDF_example.R" is a simple example that illustrates how our code works in a simulated dataset.
File "counterCDF_graph.R" is used for visualizing our estimation results and produce graphs.
