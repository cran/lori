-   [LORI (LOw-Rank Interaction) model for count data with
    covariates](#lori-low-rank-interaction-model-for-count-data-with-covariates)
    -   [Example: Aravo data set](#example-aravo-data-set)

LORI (LOw-Rank Interaction) model for count data with covariates
================================================================

The LORI model is designed to analyse count data with covariates, using
a Poisson log-linear model. In particular, it can be used to assess the
effect of temporal and geographical covariates on species abundances.

Let *Y* ∈ ℕ<sup>*n* × *p*</sup> be a (incomplete) matrix of counts, and
*L* ∈ ℝ<sup>*n**p* × *K*</sup> a matrix of covariates about the rows and
columns of *Y*. For example if *Y* counts the abundance of species
across sites (rows) and time stamps (columns), *L* might contain
temporal, spatial, and spatio-temporal information.

    library(lori)
    library(glmnet)

    ## Warning: package 'glmnet' was built under R version 3.4.4

    ## Loading required package: Matrix

    ## Warning: package 'Matrix' was built under R version 3.4.4

    ## Loading required package: foreach

    ## Loaded glmnet 2.0-16

    library(gridExtra)
    data("aravo")

Example: Aravo data set
-----------------------

The {Aravo data set} measures the abundance of 82 species of alpine
plants in 75 sites in France. The data consist of a contingency table
collecting the abundance of species across sampling sites. Covariates
about the environments and species are also available.

    # Environment characteristics
    head(aravo$env)

    ##      Aspect Slope Form PhysD ZoogD Snow
    ## AR07      7     2    1    50    no  140
    ## AR71      1    35    3    40    no  140
    ## AR26      5     0    3    20    no  140
    ## AR54      9    30    3    80    no  140
    ## AR60      9     5    1    80    no  140
    ## AR70      1    30    3    40    no  140

    # Species traits
    head(aravo$traits)

    ##           Height Spread Angle  Area Thick  SLA N_mass Seed
    ## Agro.rupe      6     10    80  60.0  0.12  8.1 218.70 0.08
    ## Alop.alpi      5     20    20 190.9  0.20 15.1 203.85 0.21
    ## Anth.nipp     15      5    50 280.0  0.08 18.0 219.60 0.54
    ## Heli.sede      0     30    80 600.0  0.20 10.6 233.20 1.72
    ## Aven.vers     12     30    60 420.0  0.14 12.5 156.25 1.17
    ## Care.rosa     30     20    80 180.0  0.40  6.5 208.65 1.68

    d <- dim(aravo$spe)
    n <- d[1]
    p <- d[2]

    # Create covariate matrix, choose quantitative variables in Row and Column covariates
    # and use covmat function to replicate species/environments
    # center and scale covariate matrix
    cov <- scale(covmat(aravo$env[, c(1,2,4,6)], aravo$traits, n, p))

    lambda1 <- qut(aravo$spe, cov)

    ## 
     1 / 100
     2 / 100
     3 / 100
     4 / 100
     5 / 100
     6 / 100
     7 / 100
     8 / 100
     9 / 100
     10 / 100
     11 / 100
     12 / 100
     13 / 100
     14 / 100
     15 / 100
     16 / 100
     17 / 100
     18 / 100
     19 / 100
     20 / 100
     21 / 100
     22 / 100
     23 / 100
     24 / 100
     25 / 100
     26 / 100
     27 / 100
     28 / 100
     29 / 100
     30 / 100
     31 / 100
     32 / 100
     33 / 100
     34 / 100
     35 / 100
     36 / 100
     37 / 100
     38 / 100
     39 / 100
     40 / 100
     41 / 100
     42 / 100
     43 / 100
     44 / 100
     45 / 100
     46 / 100
     47 / 100
     48 / 100
     49 / 100
     50 / 100
     51 / 100
     52 / 100
     53 / 100
     54 / 100
     55 / 100
     56 / 100
     57 / 100
     58 / 100
     59 / 100
     60 / 100
     61 / 100
     62 / 100
     63 / 100
     64 / 100
     65 / 100
     66 / 100
     67 / 100
     68 / 100
     69 / 100
     70 / 100
     71 / 100
     72 / 100
     73 / 100
     74 / 100
     75 / 100
     76 / 100
     77 / 100
     78 / 100
     79 / 100
     80 / 100
     81 / 100
     82 / 100
     83 / 100
     84 / 100
     85 / 100
     86 / 100
     87 / 100
     88 / 100
     89 / 100
     90 / 100
     91 / 100
     92 / 100
     93 / 100
     94 / 100
     95 / 100
     96 / 100
     97 / 100
     98 / 100
     99 / 100
     100 / 100

    res <- lori(aravo$spe, cov, lambda1=lambda1, lambda2=0, trace.it = T)

    ## [1] "fitting model..."
    ## 
    1%
    2%
    3%
    4%
    5%
    6%
    7%
    8%
    9%
    10%
    11%
    12%
    13%
    14%
    15%
    16%
    17%
    18%
    19%
    20%
    21%
    22%
    23%
    24%
    25%
    26%
    27%
    28%
    29%
    30%
    31%
    32%
    33%
    34%
    35%
    36%
    37%
    38%
    39%
    40%
    41%
    42%
    43%
    44%
    45%
    46%
    47%
    48%
    49%
    50%
    51%
    52%
    53%
    54%
    55%
    56%
    57%
    58%
    59%
    60%
    61%
    62%
    63%
    64%
    65%
    66%
    67%
    68%
    69%
    70%
    71%
    72%
    73%
    74%
    75%
    76%
    77%
    78%
    79%
    80%
    81%
    82%
    83%
    84%
    85%
    86%
    87%
    88%
    89%
    90%
    91%
    92%
    93%
    94%
    95%
    96%
    96% - iter: 10 - error: 1.036751e-05 - objective: 0.6617773
    96% - iter: 20 - error: 1.035692e-05 - objective: 0.6617088
    96% - iter: 30 - error: 1.034633e-05 - objective: 0.6616403
    96% - iter: 40 - error: 1.033574e-05 - objective: 0.6615719
    96% - iter: 50 - error: 1.032516e-05 - objective: 0.6615035
    96% - iter: 60 - error: 1.031458e-05 - objective: 0.6614353
    96% - iter: 70 - error: 1.030399e-05 - objective: 0.6613671
    96% - iter: 80 - error: 1.029342e-05 - objective: 0.661299
    96% - iter: 90 - error: 1.028284e-05 - objective: 0.661231
    96% - iter: 100 - error: 1.027227e-05 - objective: 0.661163
    96% - iter: 110 - error: 1.026169e-05 - objective: 0.6610951
    96% - iter: 120 - error: 1.025112e-05 - objective: 0.6610273
    96% - iter: 130 - error: 1.024056e-05 - objective: 0.6609596
    96% - iter: 140 - error: 1.022999e-05 - objective: 0.660892
    96% - iter: 150 - error: 1.021943e-05 - objective: 0.6608244
    96% - iter: 160 - error: 1.020887e-05 - objective: 0.6607569
    96% - iter: 170 - error: 1.019831e-05 - objective: 0.6606895
    96% - iter: 180 - error: 1.018775e-05 - objective: 0.6606222
    96% - iter: 190 - error: 1.01772e-05 - objective: 0.6605549
    96% - iter: 200 - error: 1.016664e-05 - objective: 0.6604877
    96% - iter: 210 - error: 1.015609e-05 - objective: 0.6604206
    96% - iter: 220 - error: 1.014555e-05 - objective: 0.6603536
    96% - iter: 230 - error: 1.0135e-05 - objective: 0.6602866
    96% - iter: 240 - error: 1.012446e-05 - objective: 0.6602198
    96% - iter: 250 - error: 1.011392e-05 - objective: 0.6601529
    96% - iter: 260 - error: 1.010338e-05 - objective: 0.6600862
    96% - iter: 270 - error: 1.009284e-05 - objective: 0.6600196
    96% - iter: 280 - error: 1.008231e-05 - objective: 0.659953
    96% - iter: 290 - error: 1.007178e-05 - objective: 0.6598865
    96% - iter: 300 - error: 1.006125e-05 - objective: 0.6598201
    96% - iter: 310 - error: 1.005072e-05 - objective: 0.6597537
    96% - iter: 320 - error: 1.00402e-05 - objective: 0.6596875
    96% - iter: 330 - error: 1.002968e-05 - objective: 0.6596213
    96% - iter: 340 - error: 1.001916e-05 - objective: 0.6595552
    96% - iter: 350 - error: 1.000864e-05 - objective: 0.6594891
    97%
    97% - iter: 10 - error: 1.207383e-05 - objective: 0.6591226
    97% - iter: 20 - error: 1.206137e-05 - objective: 0.6590431
    97% - iter: 30 - error: 1.204891e-05 - objective: 0.6589637
    97% - iter: 40 - error: 1.203645e-05 - objective: 0.6588843
    97% - iter: 50 - error: 1.202399e-05 - objective: 0.6588051
    97% - iter: 60 - error: 1.201154e-05 - objective: 0.6587259
    97% - iter: 70 - error: 1.199908e-05 - objective: 0.6586468
    97% - iter: 80 - error: 1.198663e-05 - objective: 0.6585678
    97% - iter: 90 - error: 1.197418e-05 - objective: 0.6584889
    97% - iter: 100 - error: 1.196172e-05 - objective: 0.6584101
    97% - iter: 110 - error: 1.194927e-05 - objective: 0.6583314
    97% - iter: 120 - error: 1.193682e-05 - objective: 0.6582528
    97% - iter: 130 - error: 1.192438e-05 - objective: 0.6581743
    97% - iter: 140 - error: 1.191193e-05 - objective: 0.6580959
    97% - iter: 150 - error: 1.189948e-05 - objective: 0.6580175
    97% - iter: 160 - error: 1.188704e-05 - objective: 0.6579393
    97% - iter: 170 - error: 1.18746e-05 - objective: 0.6578611
    97% - iter: 180 - error: 1.186215e-05 - objective: 0.657783
    97% - iter: 190 - error: 1.184971e-05 - objective: 0.6577051
    97% - iter: 200 - error: 1.183728e-05 - objective: 0.6576272
    97% - iter: 210 - error: 1.182484e-05 - objective: 0.6575494
    97% - iter: 220 - error: 1.18124e-05 - objective: 0.6574717
    97% - iter: 230 - error: 1.179997e-05 - objective: 0.6573941
    97% - iter: 240 - error: 1.178753e-05 - objective: 0.6573166
    97% - iter: 250 - error: 1.17751e-05 - objective: 0.6572391
    97% - iter: 260 - error: 1.176267e-05 - objective: 0.6571618
    97% - iter: 270 - error: 1.175024e-05 - objective: 0.6570845
    97% - iter: 280 - error: 1.173782e-05 - objective: 0.6570074
    97% - iter: 290 - error: 1.172539e-05 - objective: 0.6569303
    97% - iter: 300 - error: 1.171297e-05 - objective: 0.6568533
    97% - iter: 310 - error: 1.170055e-05 - objective: 0.6567764
    97% - iter: 320 - error: 1.168813e-05 - objective: 0.6566996
    97% - iter: 330 - error: 1.167571e-05 - objective: 0.6566229
    97% - iter: 340 - error: 1.166329e-05 - objective: 0.6565463
    97% - iter: 350 - error: 1.165087e-05 - objective: 0.6564698
    97% - iter: 360 - error: 1.163846e-05 - objective: 0.6563934
    97% - iter: 370 - error: 1.162605e-05 - objective: 0.656317
    97% - iter: 380 - error: 1.161364e-05 - objective: 0.6562408
    97% - iter: 390 - error: 1.160123e-05 - objective: 0.6561646
    97% - iter: 400 - error: 1.158882e-05 - objective: 0.6560885
    97% - iter: 410 - error: 1.157642e-05 - objective: 0.6560125
    97% - iter: 420 - error: 1.156401e-05 - objective: 0.6559366
    97% - iter: 430 - error: 1.155161e-05 - objective: 0.6558608
    97% - iter: 440 - error: 1.153921e-05 - objective: 0.6557851
    97% - iter: 450 - error: 1.152682e-05 - objective: 0.6557095
    97% - iter: 460 - error: 1.151442e-05 - objective: 0.655634
    97% - iter: 470 - error: 1.150203e-05 - objective: 0.6555585
    97% - iter: 480 - error: 1.148964e-05 - objective: 0.6554832
    97% - iter: 490 - error: 1.147725e-05 - objective: 0.6554079
    97% - iter: 500 - error: 1.146486e-05 - objective: 0.6553327
    97% - iter: 510 - error: 1.145247e-05 - objective: 0.6552577
    97% - iter: 520 - error: 1.144009e-05 - objective: 0.6551827
    97% - iter: 530 - error: 1.142771e-05 - objective: 0.6551078
    97% - iter: 540 - error: 1.141533e-05 - objective: 0.6550329
    97% - iter: 550 - error: 1.140295e-05 - objective: 0.6549582
    97% - iter: 560 - error: 1.139058e-05 - objective: 0.6548836
    97% - iter: 570 - error: 1.13782e-05 - objective: 0.654809
    97% - iter: 580 - error: 1.136583e-05 - objective: 0.6547346
    97% - iter: 590 - error: 1.135346e-05 - objective: 0.6546602
    97% - iter: 600 - error: 1.13411e-05 - objective: 0.6545859
    97% - iter: 610 - error: 1.132873e-05 - objective: 0.6545118
    97% - iter: 620 - error: 1.131637e-05 - objective: 0.6544377
    97% - iter: 630 - error: 1.130401e-05 - objective: 0.6543636
    97% - iter: 640 - error: 1.129165e-05 - objective: 0.6542897
    97% - iter: 650 - error: 1.12793e-05 - objective: 0.6542159
    97% - iter: 660 - error: 1.126695e-05 - objective: 0.6541421
    97% - iter: 670 - error: 1.12546e-05 - objective: 0.6540685
    97% - iter: 680 - error: 1.124225e-05 - objective: 0.6539949
    97% - iter: 690 - error: 1.12299e-05 - objective: 0.6539215
    97% - iter: 700 - error: 1.121756e-05 - objective: 0.6538481
    97% - iter: 710 - error: 1.120522e-05 - objective: 0.6537748
    97% - iter: 720 - error: 1.119288e-05 - objective: 0.6537016
    97% - iter: 730 - error: 1.118055e-05 - objective: 0.6536284
    97% - iter: 740 - error: 1.116821e-05 - objective: 0.6535554
    97% - iter: 750 - error: 1.115588e-05 - objective: 0.6534825
    97% - iter: 760 - error: 1.114355e-05 - objective: 0.6534096
    97% - iter: 770 - error: 1.113123e-05 - objective: 0.6533369
    97% - iter: 780 - error: 1.11189e-05 - objective: 0.6532642
    97% - iter: 790 - error: 1.110658e-05 - objective: 0.6531916
    97% - iter: 800 - error: 1.109427e-05 - objective: 0.6531191
    97% - iter: 810 - error: 1.108195e-05 - objective: 0.6530467
    97% - iter: 820 - error: 1.106964e-05 - objective: 0.6529744
    97% - iter: 830 - error: 1.105733e-05 - objective: 0.6529021
    97% - iter: 840 - error: 1.104502e-05 - objective: 0.65283
    97% - iter: 850 - error: 1.103272e-05 - objective: 0.6527579
    97% - iter: 860 - error: 1.102042e-05 - objective: 0.652686
    97% - iter: 870 - error: 1.100812e-05 - objective: 0.6526141
    97% - iter: 880 - error: 1.099582e-05 - objective: 0.6525423
    97% - iter: 890 - error: 1.098353e-05 - objective: 0.6524706
    97% - iter: 900 - error: 1.097124e-05 - objective: 0.652399
    97% - iter: 910 - error: 1.095895e-05 - objective: 0.6523274
    97% - iter: 920 - error: 1.094667e-05 - objective: 0.652256
    97% - iter: 930 - error: 1.093439e-05 - objective: 0.6521846
    97% - iter: 940 - error: 1.092211e-05 - objective: 0.6521134
    97% - iter: 950 - error: 1.090983e-05 - objective: 0.6520422
    97% - iter: 960 - error: 1.089756e-05 - objective: 0.6519711
    97% - iter: 970 - error: 1.088529e-05 - objective: 0.6519001
    97% - iter: 980 - error: 1.087302e-05 - objective: 0.6518292
    97% - iter: 990 - error: 1.086076e-05 - objective: 0.6517584
    97% - iter: 1000 - error: 1.08485e-05 - objective: 0.6516876
    98%
    98% - iter: 10 - error: 1.29738e-05 - objective: 0.6507272
    98% - iter: 20 - error: 1.295938e-05 - objective: 0.6506428
    98% - iter: 30 - error: 1.294495e-05 - objective: 0.6505585
    98% - iter: 40 - error: 1.293053e-05 - objective: 0.6504744
    98% - iter: 50 - error: 1.291611e-05 - objective: 0.6503903
    98% - iter: 60 - error: 1.290169e-05 - objective: 0.6503064
    98% - iter: 70 - error: 1.288727e-05 - objective: 0.6502225
    98% - iter: 80 - error: 1.287285e-05 - objective: 0.6501388
    98% - iter: 90 - error: 1.285844e-05 - objective: 0.6500552
    98% - iter: 100 - error: 1.284403e-05 - objective: 0.6499716
    98% - iter: 110 - error: 1.282961e-05 - objective: 0.6498882
    98% - iter: 120 - error: 1.28152e-05 - objective: 0.6498049
    98% - iter: 130 - error: 1.28008e-05 - objective: 0.6497217
    98% - iter: 140 - error: 1.278639e-05 - objective: 0.6496386
    98% - iter: 150 - error: 1.277198e-05 - objective: 0.6495555
    98% - iter: 160 - error: 1.275758e-05 - objective: 0.6494726
    98% - iter: 170 - error: 1.274318e-05 - objective: 0.6493898
    98% - iter: 180 - error: 1.272878e-05 - objective: 0.6493071
    98% - iter: 190 - error: 1.271438e-05 - objective: 0.6492246
    98% - iter: 200 - error: 1.269999e-05 - objective: 0.6491421
    98% - iter: 210 - error: 1.268559e-05 - objective: 0.6490597
    98% - iter: 220 - error: 1.26712e-05 - objective: 0.6489774
    98% - iter: 230 - error: 1.265681e-05 - objective: 0.6488952
    98% - iter: 240 - error: 1.264242e-05 - objective: 0.6488132
    98% - iter: 250 - error: 1.262804e-05 - objective: 0.6487312
    98% - iter: 260 - error: 1.261365e-05 - objective: 0.6486493
    98% - iter: 270 - error: 1.259927e-05 - objective: 0.6485676
    98% - iter: 280 - error: 1.258489e-05 - objective: 0.6484859
    98% - iter: 290 - error: 1.257052e-05 - objective: 0.6484043
    98% - iter: 300 - error: 1.255614e-05 - objective: 0.6483229
    98% - iter: 310 - error: 1.254177e-05 - objective: 0.6482415
    98% - iter: 320 - error: 1.25274e-05 - objective: 0.6481603
    98% - iter: 330 - error: 1.251303e-05 - objective: 0.6480792
    98% - iter: 340 - error: 1.249866e-05 - objective: 0.6479981
    98% - iter: 350 - error: 1.24843e-05 - objective: 0.6479172
    98% - iter: 360 - error: 1.246994e-05 - objective: 0.6478364
    98% - iter: 370 - error: 1.245558e-05 - objective: 0.6477556
    98% - iter: 380 - error: 1.244122e-05 - objective: 0.647675
    98% - iter: 390 - error: 1.242687e-05 - objective: 0.6475945
    98% - iter: 400 - error: 1.241252e-05 - objective: 0.6475141
    98% - iter: 410 - error: 1.239817e-05 - objective: 0.6474338
    98% - iter: 420 - error: 1.238382e-05 - objective: 0.6473535
    98% - iter: 430 - error: 1.236948e-05 - objective: 0.6472734
    98% - iter: 440 - error: 1.235514e-05 - objective: 0.6471934
    98% - iter: 450 - error: 1.23408e-05 - objective: 0.6471135
    98% - iter: 460 - error: 1.232646e-05 - objective: 0.6470337
    98% - iter: 470 - error: 1.231213e-05 - objective: 0.646954
    98% - iter: 480 - error: 1.22978e-05 - objective: 0.6468744
    98% - iter: 490 - error: 1.228414e-05 - objective: 0.6467949
    98% - iter: 500 - error: 1.226914e-05 - objective: 0.6467155
    98% - iter: 510 - error: 1.225482e-05 - objective: 0.6466362
    98% - iter: 520 - error: 1.22405e-05 - objective: 0.646557
    98% - iter: 530 - error: 1.222618e-05 - objective: 0.646478
    98% - iter: 540 - error: 1.221187e-05 - objective: 0.646399
    98% - iter: 550 - error: 1.219756e-05 - objective: 0.6463201
    98% - iter: 560 - error: 1.218325e-05 - objective: 0.6462413
    98% - iter: 570 - error: 1.216895e-05 - objective: 0.6461626
    98% - iter: 580 - error: 1.215465e-05 - objective: 0.6460841
    98% - iter: 590 - error: 1.214035e-05 - objective: 0.6460056
    98% - iter: 600 - error: 1.212606e-05 - objective: 0.6459272
    98% - iter: 610 - error: 1.211176e-05 - objective: 0.645849
    98% - iter: 620 - error: 1.209748e-05 - objective: 0.6457708
    98% - iter: 630 - error: 1.208319e-05 - objective: 0.6456927
    98% - iter: 640 - error: 1.206891e-05 - objective: 0.6456148
    98% - iter: 650 - error: 1.205463e-05 - objective: 0.6455369
    98% - iter: 660 - error: 1.204036e-05 - objective: 0.6454591
    98% - iter: 670 - error: 1.202608e-05 - objective: 0.6453815
    98% - iter: 680 - error: 1.201182e-05 - objective: 0.6453039
    98% - iter: 690 - error: 1.199755e-05 - objective: 0.6452265
    98% - iter: 700 - error: 1.198329e-05 - objective: 0.6451491
    98% - iter: 710 - error: 1.196903e-05 - objective: 0.6450718
    98% - iter: 720 - error: 1.195478e-05 - objective: 0.6449947
    98% - iter: 730 - error: 1.194053e-05 - objective: 0.6449176
    98% - iter: 740 - error: 1.192628e-05 - objective: 0.6448407
    98% - iter: 750 - error: 1.191203e-05 - objective: 0.6447638
    98% - iter: 760 - error: 1.189779e-05 - objective: 0.6446871
    98% - iter: 770 - error: 1.188356e-05 - objective: 0.6446104
    98% - iter: 780 - error: 1.186932e-05 - objective: 0.6445339
    98% - iter: 790 - error: 1.185509e-05 - objective: 0.6444574
    98% - iter: 800 - error: 1.184087e-05 - objective: 0.6443811
    98% - iter: 810 - error: 1.182665e-05 - objective: 0.6443049
    98% - iter: 820 - error: 1.181243e-05 - objective: 0.6442287
    98% - iter: 830 - error: 1.179822e-05 - objective: 0.6441527
    98% - iter: 840 - error: 1.178401e-05 - objective: 0.6440767
    98% - iter: 850 - error: 1.17698e-05 - objective: 0.6440009
    98% - iter: 860 - error: 1.17556e-05 - objective: 0.6439251
    98% - iter: 870 - error: 1.17414e-05 - objective: 0.6438495
    98% - iter: 880 - error: 1.172721e-05 - objective: 0.643774
    98% - iter: 890 - error: 1.171301e-05 - objective: 0.6436985
    98% - iter: 900 - error: 1.169883e-05 - objective: 0.6436232
    98% - iter: 910 - error: 1.168465e-05 - objective: 0.6435479
    98% - iter: 920 - error: 1.167047e-05 - objective: 0.6434728
    98% - iter: 930 - error: 1.16563e-05 - objective: 0.6433977
    98% - iter: 940 - error: 1.164213e-05 - objective: 0.6433228
    98% - iter: 950 - error: 1.162796e-05 - objective: 0.643248
    98% - iter: 960 - error: 1.16138e-05 - objective: 0.6431732
    98% - iter: 970 - error: 1.159964e-05 - objective: 0.6430986
    98% - iter: 980 - error: 1.158549e-05 - objective: 0.643024
    98% - iter: 990 - error: 1.157134e-05 - objective: 0.6429496
    98% - iter: 1000 - error: 1.15572e-05 - objective: 0.6428753
    99%
    99% - iter: 10 - error: 1.373306e-05 - objective: 0.6413019
    99% - iter: 20 - error: 1.371648e-05 - objective: 0.6412139
    99% - iter: 30 - error: 1.369991e-05 - objective: 0.641126
    99% - iter: 40 - error: 1.368334e-05 - objective: 0.6410382
    99% - iter: 50 - error: 1.366677e-05 - objective: 0.6409506
    99% - iter: 60 - error: 1.365021e-05 - objective: 0.6408631
    99% - iter: 70 - error: 1.363364e-05 - objective: 0.6407757
    99% - iter: 80 - error: 1.361709e-05 - objective: 0.6406884
    99% - iter: 90 - error: 1.360053e-05 - objective: 0.6406012
    99% - iter: 100 - error: 1.358398e-05 - objective: 0.6405141
    99% - iter: 110 - error: 1.356743e-05 - objective: 0.6404272
    99% - iter: 120 - error: 1.355088e-05 - objective: 0.6403403
    99% - iter: 130 - error: 1.353434e-05 - objective: 0.6402536
    99% - iter: 140 - error: 1.35178e-05 - objective: 0.6401671
    99% - iter: 150 - error: 1.350126e-05 - objective: 0.6400806
    99% - iter: 160 - error: 1.348472e-05 - objective: 0.6399942
    99% - iter: 170 - error: 1.346819e-05 - objective: 0.639908
    99% - iter: 180 - error: 1.345166e-05 - objective: 0.6398219
    99% - iter: 190 - error: 1.343514e-05 - objective: 0.6397359
    99% - iter: 200 - error: 1.341862e-05 - objective: 0.63965
    99% - iter: 210 - error: 1.34021e-05 - objective: 0.6395642
    99% - iter: 220 - error: 1.338559e-05 - objective: 0.6394786
    99% - iter: 230 - error: 1.336908e-05 - objective: 0.639393
    99% - iter: 240 - error: 1.335257e-05 - objective: 0.6393076
    99% - iter: 250 - error: 1.333606e-05 - objective: 0.6392223
    99% - iter: 260 - error: 1.331956e-05 - objective: 0.6391371
    99% - iter: 270 - error: 1.330307e-05 - objective: 0.6390521
    99% - iter: 280 - error: 1.328658e-05 - objective: 0.6389671
    99% - iter: 290 - error: 1.327009e-05 - objective: 0.6388823
    99% - iter: 300 - error: 1.32536e-05 - objective: 0.6387976
    99% - iter: 310 - error: 1.323712e-05 - objective: 0.638713
    99% - iter: 320 - error: 1.322064e-05 - objective: 0.6386285
    99% - iter: 330 - error: 1.320417e-05 - objective: 0.6385441
    99% - iter: 340 - error: 1.31877e-05 - objective: 0.6384599
    99% - iter: 350 - error: 1.317124e-05 - objective: 0.6383757
    99% - iter: 360 - error: 1.315478e-05 - objective: 0.6382917
    99% - iter: 370 - error: 1.313832e-05 - objective: 0.6382078
    99% - iter: 380 - error: 1.312187e-05 - objective: 0.638124
    99% - iter: 390 - error: 1.310542e-05 - objective: 0.6380404
    99% - iter: 400 - error: 1.308897e-05 - objective: 0.6379568
    99% - iter: 410 - error: 1.307253e-05 - objective: 0.6378734
    99% - iter: 420 - error: 1.30561e-05 - objective: 0.63779
    99% - iter: 430 - error: 1.303967e-05 - objective: 0.6377068
    99% - iter: 440 - error: 1.302324e-05 - objective: 0.6376237
    99% - iter: 450 - error: 1.300682e-05 - objective: 0.6375408
    99% - iter: 460 - error: 1.29904e-05 - objective: 0.6374579
    99% - iter: 470 - error: 1.297399e-05 - objective: 0.6373752
    99% - iter: 480 - error: 1.295758e-05 - objective: 0.6372925
    99% - iter: 490 - error: 1.294117e-05 - objective: 0.63721
    99% - iter: 500 - error: 1.292477e-05 - objective: 0.6371276
    99% - iter: 510 - error: 1.290838e-05 - objective: 0.6370453
    99% - iter: 520 - error: 1.289199e-05 - objective: 0.6369632
    99% - iter: 530 - error: 1.28756e-05 - objective: 0.6368811
    99% - iter: 540 - error: 1.285922e-05 - objective: 0.6367992
    99% - iter: 550 - error: 1.284285e-05 - objective: 0.6367173
    99% - iter: 560 - error: 1.282648e-05 - objective: 0.6366356
    99% - iter: 570 - error: 1.281011e-05 - objective: 0.636554
    99% - iter: 580 - error: 1.279375e-05 - objective: 0.6364726
    99% - iter: 590 - error: 1.27774e-05 - objective: 0.6363912
    99% - iter: 600 - error: 1.276105e-05 - objective: 0.6363099
    99% - iter: 610 - error: 1.274471e-05 - objective: 0.6362288
    99% - iter: 620 - error: 1.272837e-05 - objective: 0.6361478
    99% - iter: 630 - error: 1.271203e-05 - objective: 0.6360669
    99% - iter: 640 - error: 1.26957e-05 - objective: 0.6359861
    99% - iter: 650 - error: 1.267938e-05 - objective: 0.6359054
    99% - iter: 660 - error: 1.266306e-05 - objective: 0.6358248
    99% - iter: 670 - error: 1.264675e-05 - objective: 0.6357444
    99% - iter: 680 - error: 1.263044e-05 - objective: 0.635664
    99% - iter: 690 - error: 1.261414e-05 - objective: 0.6355838
    99% - iter: 700 - error: 1.259784e-05 - objective: 0.6355037
    99% - iter: 710 - error: 1.258155e-05 - objective: 0.6354237
    99% - iter: 720 - error: 1.256527e-05 - objective: 0.6353438
    99% - iter: 730 - error: 1.254899e-05 - objective: 0.6352641
    99% - iter: 740 - error: 1.253272e-05 - objective: 0.6351844
    99% - iter: 750 - error: 1.251645e-05 - objective: 0.6351049
    99% - iter: 760 - error: 1.250019e-05 - objective: 0.6350254
    99% - iter: 770 - error: 1.248393e-05 - objective: 0.6349461
    99% - iter: 780 - error: 1.246768e-05 - objective: 0.6348669
    99% - iter: 790 - error: 1.245144e-05 - objective: 0.6347878
    99% - iter: 800 - error: 1.24352e-05 - objective: 0.6347088
    99% - iter: 810 - error: 1.241897e-05 - objective: 0.63463
    99% - iter: 820 - error: 1.240275e-05 - objective: 0.6345512
    99% - iter: 830 - error: 1.238653e-05 - objective: 0.6344726
    99% - iter: 840 - error: 1.237031e-05 - objective: 0.6343941
    99% - iter: 850 - error: 1.235411e-05 - objective: 0.6343156
    99% - iter: 860 - error: 1.233791e-05 - objective: 0.6342373
    99% - iter: 870 - error: 1.232171e-05 - objective: 0.6341591
    99% - iter: 880 - error: 1.230552e-05 - objective: 0.6340811
    99% - iter: 890 - error: 1.228934e-05 - objective: 0.6340031
    99% - iter: 900 - error: 1.227317e-05 - objective: 0.6339253
    99% - iter: 910 - error: 1.2257e-05 - objective: 0.6338475
    99% - iter: 920 - error: 1.224084e-05 - objective: 0.6337699
    99% - iter: 930 - error: 1.222468e-05 - objective: 0.6336924
    99% - iter: 940 - error: 1.220853e-05 - objective: 0.633615
    99% - iter: 950 - error: 1.219239e-05 - objective: 0.6335377
    99% - iter: 960 - error: 1.217626e-05 - objective: 0.6334605
    99% - iter: 970 - error: 1.216013e-05 - objective: 0.6333834
    99% - iter: 980 - error: 1.214401e-05 - objective: 0.6333065
    99% - iter: 990 - error: 1.212789e-05 - objective: 0.6332296
    99% - iter: 1000 - error: 1.211178e-05 - objective: 0.6331529
    100%
    100% - iter: 10 - error: 1.431635e-05 - objective: 0.6309947
    100% - iter: 20 - error: 1.429754e-05 - objective: 0.6309045
    100% - iter: 30 - error: 1.427873e-05 - objective: 0.6308143
    100% - iter: 40 - error: 1.425993e-05 - objective: 0.6307243
    100% - iter: 50 - error: 1.424114e-05 - objective: 0.6306345
    100% - iter: 60 - error: 1.422235e-05 - objective: 0.6305447
    100% - iter: 70 - error: 1.420357e-05 - objective: 0.6304551
    100% - iter: 80 - error: 1.418479e-05 - objective: 0.6303656
    100% - iter: 90 - error: 1.416602e-05 - objective: 0.6302763
    100% - iter: 100 - error: 1.414726e-05 - objective: 0.6301871
    100% - iter: 110 - error: 1.41285e-05 - objective: 0.630098
    100% - iter: 120 - error: 1.410974e-05 - objective: 0.6300091
    100% - iter: 130 - error: 1.4091e-05 - objective: 0.6299202
    100% - iter: 140 - error: 1.407226e-05 - objective: 0.6298315
    100% - iter: 150 - error: 1.405352e-05 - objective: 0.629743
    100% - iter: 160 - error: 1.403479e-05 - objective: 0.6296546
    100% - iter: 170 - error: 1.401607e-05 - objective: 0.6295663
    100% - iter: 180 - error: 1.399735e-05 - objective: 0.6294781
    100% - iter: 190 - error: 1.397864e-05 - objective: 0.6293901
    100% - iter: 200 - error: 1.395994e-05 - objective: 0.6293021
    100% - iter: 210 - error: 1.394124e-05 - objective: 0.6292144
    100% - iter: 220 - error: 1.392255e-05 - objective: 0.6291267
    100% - iter: 230 - error: 1.390387e-05 - objective: 0.6290392
    100% - iter: 240 - error: 1.388519e-05 - objective: 0.6289518
    100% - iter: 250 - error: 1.386652e-05 - objective: 0.6288645
    100% - iter: 260 - error: 1.384785e-05 - objective: 0.6287774
    100% - iter: 270 - error: 1.38292e-05 - objective: 0.6286904
    100% - iter: 280 - error: 1.381055e-05 - objective: 0.6286035
    100% - iter: 290 - error: 1.37919e-05 - objective: 0.6285168
    100% - iter: 300 - error: 1.377327e-05 - objective: 0.6284302
    100% - iter: 310 - error: 1.375464e-05 - objective: 0.6283437
    100% - iter: 320 - error: 1.373601e-05 - objective: 0.6282573
    100% - iter: 330 - error: 1.37174e-05 - objective: 0.6281711
    100% - iter: 340 - error: 1.369879e-05 - objective: 0.628085
    100% - iter: 350 - error: 1.368019e-05 - objective: 0.6279991
    100% - iter: 360 - error: 1.36616e-05 - objective: 0.6279132
    100% - iter: 370 - error: 1.364301e-05 - objective: 0.6278275
    100% - iter: 380 - error: 1.362443e-05 - objective: 0.6277419
    100% - iter: 390 - error: 1.360586e-05 - objective: 0.6276565
    100% - iter: 400 - error: 1.35873e-05 - objective: 0.6275711
    100% - iter: 410 - error: 1.356874e-05 - objective: 0.6274859
    100% - iter: 420 - error: 1.355019e-05 - objective: 0.6274009
    100% - iter: 430 - error: 1.353165e-05 - objective: 0.6273159
    100% - iter: 440 - error: 1.351312e-05 - objective: 0.6272311
    100% - iter: 450 - error: 1.349459e-05 - objective: 0.6271464
    100% - iter: 460 - error: 1.347608e-05 - objective: 0.6270619
    100% - iter: 470 - error: 1.345757e-05 - objective: 0.6269774
    100% - iter: 480 - error: 1.343907e-05 - objective: 0.6268931
    100% - iter: 490 - error: 1.342058e-05 - objective: 0.6268089
    100% - iter: 500 - error: 1.340209e-05 - objective: 0.6267249
    100% - iter: 510 - error: 1.338361e-05 - objective: 0.626641
    100% - iter: 520 - error: 1.336515e-05 - objective: 0.6265572
    100% - iter: 530 - error: 1.334669e-05 - objective: 0.6264735
    100% - iter: 540 - error: 1.332823e-05 - objective: 0.6263899
    100% - iter: 550 - error: 1.330979e-05 - objective: 0.6263065
    100% - iter: 560 - error: 1.329136e-05 - objective: 0.6262232
    100% - iter: 570 - error: 1.327293e-05 - objective: 0.6261401
    100% - iter: 580 - error: 1.325451e-05 - objective: 0.626057
    100% - iter: 590 - error: 1.32361e-05 - objective: 0.6259741
    100% - iter: 600 - error: 1.32177e-05 - objective: 0.6258913
    100% - iter: 610 - error: 1.319931e-05 - objective: 0.6258087
    100% - iter: 620 - error: 1.318093e-05 - objective: 0.6257261
    100% - iter: 630 - error: 1.316256e-05 - objective: 0.6256437
    100% - iter: 640 - error: 1.314419e-05 - objective: 0.6255615
    100% - iter: 650 - error: 1.312584e-05 - objective: 0.6254793
    100% - iter: 660 - error: 1.310749e-05 - objective: 0.6253973
    100% - iter: 670 - error: 1.308915e-05 - objective: 0.6253154
    100% - iter: 680 - error: 1.307082e-05 - objective: 0.6252336
    100% - iter: 690 - error: 1.305251e-05 - objective: 0.6251519
    100% - iter: 700 - error: 1.30342e-05 - objective: 0.6250704
    100% - iter: 710 - error: 1.30159e-05 - objective: 0.624989
    100% - iter: 720 - error: 1.29976e-05 - objective: 0.6249077
    100% - iter: 730 - error: 1.297932e-05 - objective: 0.6248266
    100% - iter: 740 - error: 1.296105e-05 - objective: 0.6247455
    100% - iter: 750 - error: 1.294279e-05 - objective: 0.6246646
    100% - iter: 760 - error: 1.292454e-05 - objective: 0.6245839
    100% - iter: 770 - error: 1.290629e-05 - objective: 0.6245032
    100% - iter: 780 - error: 1.288806e-05 - objective: 0.6244227
    100% - iter: 790 - error: 1.286984e-05 - objective: 0.6243423
    100% - iter: 800 - error: 1.285162e-05 - objective: 0.624262
    100% - iter: 810 - error: 1.283342e-05 - objective: 0.6241818
    100% - iter: 820 - error: 1.281522e-05 - objective: 0.6241018
    100% - iter: 830 - error: 1.279704e-05 - objective: 0.6240219
    100% - iter: 840 - error: 1.277887e-05 - objective: 0.6239421
    100% - iter: 850 - error: 1.27607e-05 - objective: 0.6238624
    100% - iter: 860 - error: 1.274255e-05 - objective: 0.6237829
    100% - iter: 870 - error: 1.272441e-05 - objective: 0.6237035
    100% - iter: 880 - error: 1.270627e-05 - objective: 0.6236242
    100% - iter: 890 - error: 1.268815e-05 - objective: 0.623545
    100% - iter: 900 - error: 1.267004e-05 - objective: 0.6234659
    100% - iter: 910 - error: 1.265193e-05 - objective: 0.623387
    100% - iter: 920 - error: 1.263384e-05 - objective: 0.6233082
    100% - iter: 930 - error: 1.261576e-05 - objective: 0.6232295
    100% - iter: 940 - error: 1.259769e-05 - objective: 0.623151
    100% - iter: 950 - error: 1.257963e-05 - objective: 0.6230725
    100% - iter: 960 - error: 1.256158e-05 - objective: 0.6229942
    100% - iter: 970 - error: 1.254355e-05 - objective: 0.622916
    100% - iter: 980 - error: 1.252552e-05 - objective: 0.622838
    100% - iter: 990 - error: 1.25075e-05 - objective: 0.62276
    100% - iter: 1000 - error: 1.24895e-05 - objective: 0.6226822

    plot_cov(res)

![](example_files/figure-markdown_strict/apply-lori-1.png)

    plot_counts(res)

![](example_files/figure-markdown_strict/apply-lori-2.png)
