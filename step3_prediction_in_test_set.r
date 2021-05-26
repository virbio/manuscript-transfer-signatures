# code inspired from https://github.com/NikolayOskolkov/ClusteringHighDimensions/blob/master/easy_scrnaseq_tsne_cluster.R
# see great blog at https://towardsdatascience.com/how-to-cluster-in-high-dimensions-4ef693bacc6

# function that runs the unsupervised clustering and cluster inference; should be run for each training dataset - test set pair
# where the inputs are:
# OUTPUT_DIRECTORY - the directory where the outputs will be written (full path should be provided and the parent directory should already exist);
                   # the name should be specific to both the training and testing datasets;
# MATCOUNTS_TEST - the normalized gene count matrix from the test dataset
                 # where there is 1 row PER GENE (IN THE TRANSFER SIGNATURE) and 1 column per sample;
                 # see end of script step1_signature_evaluation_in_training.r (commented out section) for how to generate MATCOUNTS_TEST object
# MATCOUNTS_TRAIN - the normalized gene count matrix from the training dataset
                 # where there is 1 row PER GENE (IN THE TRANSFER SIGNATURE) and 1 column per sample;
                 # see end of script step1_signature_evaluation_in_training.r (commented out section) for how to generate MATCOUNTS_TRAIN object
# LABELS_TEST - the categorical phenotypes of the samples in the test dataset (if known),
              # must be provided in the same order than the test samples or column names of MATCOUNTS_TEST
# LABELS_TRAIN - the categorical phenotypes of the samples in the training dataset (2 categories),
              # must be provided in the same order than the training samples or column names of MATCOUNTS_TRAIN
# FILE_NAME - common part of the name of the output files
# PHENOTYPE_OF_INTEREST - a character string indicating the phenotype of interest in the training dataset
                        # (or the most likely to resemble the phenotype of interest in test dataset)
                        # has to match either one of the 2 categories provided in LABELS_TRAIN.
                        # the inferred cluster in the test set that resembles the most this category in the training set will be named "case",
                        # while the other cluster(s) will be named ctrl (for control) and interm (in the presence of 3 clusters).
# n_epochs, ntree, nn, learning_rate, metric, init, nn_method are all parameters for the umap functino that can be tuned if desired

# where the output is:
# <OUTPUT_DIRECTORY>/logs/<FILE_NAME>.log  # containing the standard output of the script and informing on the steps during the process
# <OUTPUT_DIRECTORY>/models/UMAP_model_<FILE_NAME>  # containing the UMAP model
# <OUTPUT_DIRECTORY>/clusters/UMAP_CLUSTER_ASSIGNMENT_<FILE_NAME>.txt  # table containing 8 columns
                                                                       # column 1: contains the test dataset sample ids
                                                                       # column 2: contains the X coordinates of the UMAP projection
                                                                       # column 3: contains the Y coordinates of the UMAP projection
                                                                       # column 4: contains the cluster number output by HDBSCAN and knn
                                                                       # column 5: contains the membership probability for the sample to belong in the cluster
                                                                       # column 6: contains the outlier score of each sample
                                                                       # column 7: contains the true labels of the test dataset samples (if known)
                                                                       # column 8: contains the inferred label category for each sample


umap_cluster=function(OUTPUT_DIRECTORY, MATCOUNTS_TEST, MATCOUNTS_TRAIN, LABELS_TEST='', LABELS_TRAIN, FILE_NAME, PHENOTYPE_OF_INTEREST,
                      n_epochs = 1000, ntree = 500, nn=25, learning_rate = .5, metric = 'manhattan', init='laplacian', nn_method = "annoy")
{
    library("dbscan")
    library(uwot)
    library("foreach")
    library("doParallel")
    library("parallel")
    library('class')
    numCores = detectCores()
    registerDoParallel(numCores)

    cross.entropy = function(p, phat){
      x = 0
      for (i in 1:length(p)){
        x = x + (p[i] * log(phat[i]))
      }
      return(-x)
    }

    get = get0

    setwd(OUTPUT_DIRECTORY)
    dir.create('logs')
    dir.create('clusters')
    dir.create('models')

    ################################################## READING DATA ##################################################
    sink(paste0("logs/",FILE_NAME,'.log'), split=TRUE)
    print(paste0("START WORKING WITH FILE ", FILE_NAME))
    print("READING AND FILTERING DATA")
    expr = suppressWarnings(as.matrix(MATCOUNTS_TEST))
    expr = na.omit(expr)
    N_smpl = dim(expr)[2]
    if (LABELS_TEST==''){LABELS_TEST = rep('Unknown', N_smpl)}
    print(paste0("DATASET CONTAINS ",dim(expr)[1]," GENES AND ",N_smpl," SAMPLES"))

    ################################### SELECTING SIGNIFICANT PRINCIPAL COMPONENTS ###################################
    # will be used at the end to pick the UMAP that has the closest distances to the distances displayed on the PC
    print("PERFORMING PRINCIPAL COMPONENT ANALYSIS")
    PC = prcomp(t(expr), center=TRUE, scale=FALSE)
    expl_var = PC$sdev^2/sum(PC$sdev^2)

    print("SELECTING SIGNIFICANT PRINCIPAL COMPONENTS")
    if(N_smpl <= 100){N_perm = 50}else{N_perm = 20}
    expl_var_perm = matrix(NA, ncol = length(PC$sdev), nrow = N_perm)
    for(k in 1:N_perm)
    {
      expr_perm = apply(expr,2,sample)
      PC_perm = prcomp(t(expr_perm), center=TRUE, scale=FALSE)
      expl_var_perm[k,] = PC_perm$sdev^2/sum(PC_perm$sdev^2)
      print(paste0("FINISHED ",k," PERMUTATIONS"))
    }
    pval = apply(t(expl_var_perm) >= expl_var,1,sum) / N_perm
    optPC=head(which(pval>=0.05),1)-1 # optPC stands for optimal number of PC
    print(paste0("OPTIMAL NUMBER OF PRINCIPAL COMPONENTS = ",optPC, ", THEY TOGETHER EXPLAIN ",round(sum(expl_var[1:optPC])*100,0),"% OF VARIANCE"))

    #################################### PERFORMING HDBSCAN HYPERPARAMETER TUNING on UMAP #####################################
    print("PERFORMING HDBSCAN HYPERPARAMETER TUNING")
    if(N_smpl <= 50 ){N_iter = 100}else{N_iter = 50}
    if(N_smpl <= 200){N_pt = min(20,floor(N_smpl/2))}else{N_pt = 50}
    score_umap = vector(length = length(max(2,ceiling(N_smpl/15)): max(N_pt, (ceiling(N_smpl/15)+1))))

    for(i in max(2,ceiling(N_smpl/15)): max(N_pt, (ceiling(N_smpl/15)+1)))
    {
        print(paste0("MIN POINTS PER CLUSTER = ",i))
        score_iter_umap = foreach (1:N_iter, .combine = c) %dopar%
        {
            umap_iter = umap(t(expr),
                      init=init,
                      nn_method = nn_method,
                      metric = metric,
                      n_neighbors=nn,
                      n_trees=ntree,
                      n_epochs = n_epochs,
                      learning_rate = learning_rate,
                      ret_model = FALSE, ret_nn = FALSE,
                      verbose = FALSE,
                      n_threads = numCores/2,
                      n_sgd_threads=1, tmpdir=".", n_components=2, approx_pow=T)

            res_umap = hdbscan(umap_iter, minPts = i)
            if (max(res_umap$cluster)==0){penalty=25}else{penalty=max(res_umap$cluster)^2}
            score_iter_umap_temp = ((sum(res_umap$cluster==0))/ N_smpl) * penalty
            return(score_iter_umap_temp)
        }
        score_umap[i - max(2, ceiling(N_smpl/15)) +1] = mean(score_iter_umap, na.rm = TRUE)
    }

    names(score_umap) = seq(from = max(2,ceiling(N_smpl/15)), to=max(N_pt, (ceiling(N_smpl/15)+1)), by=1)
    opt_score_umap    = as.numeric(names(score_umap)[score_umap==min(score_umap)])[1]

    print(paste0("OPTIMAL MIN SIZE OF CLUSTERS WITH UMAP = ", opt_score_umap))

    #################################### PERFORMING FINAL clustering with tuned parameter #####################################
    print("PERFORMING FINAL UMAP DIMENSIONALITY REDUCTION AND HDBSCAN CLUSTERING")

    if(N_smpl <= 50){N_umap = 50}else{N_umap = 20}
    umap_out = list(length = N_umap)
    res_iter_out = list(length = N_umap)
    CE = rep(NA, N_umap)
    for(k in 1:N_umap)
    {
        print(k)
        res_iter_tmp_clusterNb = 4 # allow a maximum of 3 clusters at this step
        index = 0
        opt_score_umap_tmp = opt_score_umap
        while (res_iter_tmp_clusterNb > 3 | res_iter_tmp_clusterNb == 0)
        {
            index = index +1
            # repeat until the number of clusters is lower than 4 to facilitate interpretation
            # if it forces multiple clusters together, it is okay as we can always rerun the algorithm sequentially
            umap_iter_tmp = umap(t(expr),
                            init=init,
                            nn_method = nn_method,
                            metric = metric,
                            n_neighbors=nn,
                            n_trees=ntree,
                            n_epochs = n_epochs,
                            learning_rate = learning_rate,
                            ret_model = TRUE, ret_nn = TRUE,
                            verbose = FALSE,
                            n_threads = numCores/2,
                            n_sgd_threads=1, tmpdir=".", n_components=2, approx_pow=T)
            if (index %% 10 ==0) # increase the minimal number of samples per cluster it the number of cluster after 10 iteration is still above 3
            {
                opt_score_umap_tmp = opt_score_umap_tmp + 1
            }
            res_iter_tmp = hdbscan(umap_iter_tmp$embedding[,1:2], minPts = opt_score_umap_tmp)
            res_iter_tmp_clusterNb = max(res_iter_tmp$cluster)
        }

      umap_out[[k]] = umap_iter_tmp
      res_iter_out[[k]] = res_iter_tmp
      CE[k] = cross.entropy(dist(PC$x[,1:optPC]), dist(umap_out[[k]]$embedding[,1:2]))
    }
    names(CE)=c(1:N_umap)
    CE = CE[!is.na(CE)]


    opt_umap = umap_out[[as.numeric(names(CE)[CE==max(CE, na.rm=T)][1])]]$embedding[,1:2]
    res_opt_umap = res_iter_out[[as.numeric(names(CE)[CE==max(CE, na.rm=T)][1])]]
    setwd('./models')
    save_uwot(umap_out[[as.numeric(names(CE)[CE==max(CE, na.rm=T)][1])]], file=paste0("UMAP_model_",FILE_NAME))
    setwd(OUTPUT_DIRECTORY)

    print(res_opt_umap)

    print("USING KNN FOR CLASSIFYING HDBSCAN OUTLIERS")
    if(length(res_opt_umap$cluster[res_opt_umap$cluster==0]) > 0 & length(res_opt_umap$cluster[res_opt_umap$cluster==1]) > 0)
    {
      res_opt_umap$cluster[res_opt_umap$cluster==0] = class::knn(opt_umap[which(res_opt_umap$cluster!=0),],
                                                      opt_umap[which(res_opt_umap$cluster==0),],
                                                      res_opt_umap$cluster[res_opt_umap$cluster!=0], k=3)
    }

    N_clust_umap = length(sort(unique(res_opt_umap$cluster)))
    dfout_umap = data.frame(SAMPLE=colnames(expr),
                            COORD_X=opt_umap[,1],
                            COORD_Y=opt_umap[,2],
                            CLUSTER=res_opt_umap$cluster,
                            MEMBERSHIP_PROB=res_opt_umap$membership_prob,
                            OUTLIER_SCORE=res_opt_umap$outlier_scores,
                            TRUE_LABELS=LABELS_TEST)

    rownames(dfout_umap) = colnames(expr)

    ################################### CLUSTER INFERENCE ###################################
    print("INFERRING CLUSTER PHENOTYPE CATEGORY BY COMPARISON WITH TRAINING DATASET")

    dfTrain_case = MATCOUNTS_TRAIN[, which(LABELS_TRAIN == PHENOTYPE_OF_INTEREST)]
    dfTrain_ctrl = MATCOUNTS_TRAIN[, which(LABELS_TRAIN != PHENOTYPE_OF_INTEREST)]

    rm(list=ls(pattern="dfTest_grp"))
    for (idx in sort(unique(dfout_umap$CLUSTER)))
    {
        assign(paste0('dfTest_grp',idx), MATCOUNTS_TEST[,rownames(subset(dfout_umap, CLUSTER==idx))]) # select only sample from a cluster==idx
    }
    dfTestGroups = ls(pattern="dfTest_grp")
    #create dataframe where we will monitor up and down direction of expression (updn stands for up and down)
    df_updn = matrix(NA, ncol=1+length(dfTestGroups), nrow=length(intersect(rownames(MATCOUNTS_TEST), rownames(MATCOUNTS_TRAIN))))
    rownames(df_updn) = intersect(rownames(MATCOUNTS_TEST), rownames(MATCOUNTS_TRAIN))
    colnames(df_updn) = c('Train_case', dfTestGroups)
    # create a dataframe where we will monitor the absolute difference in median expression (in case the direction is not sufficient to infer the cluster)
    df_absdiff = df_updn
    for (gene in intersect(rownames(MATCOUNTS_TEST), rownames(MATCOUNTS_TRAIN)))
    {
        # in the df_updn dataframe,
        # 1 stands for higher expression
        # 0 stands for lower expression
        # 0.5 stands for intermediate expression
        if (median(unlist(dfTrain_case[gene,])) > median(unlist(dfTrain_ctrl[gene,]))){tmpTrain=1}else{tmpTrain=0}
        df_updn[gene,'Train_case'] = tmpTrain
        df_absdiff[gene,'Train_case'] = abs(median(unlist(dfTrain_case[gene,])) - median(unlist(dfTrain_ctrl[gene,])))
        tmpvec = c()
        for (dfTestGrp in dfTestGroups)
        {
            tmpvec = c(tmpvec, median(unlist(get(paste(dfTestGrp))[gene,])))
        }
        tmpvec2 = replace(tmpvec, which(tmpvec==max(tmpvec))[1],1)
        tmpvec2 = replace(tmpvec2, which(tmpvec==min(tmpvec))[1],0)
        tmpvec2 = replace(tmpvec2, which(tmpvec2 !=0 & tmpvec2!=1), 0.5) # in the situation where there are 3 clusters
        df_updn[gene, c(2:ncol(df_updn))] = tmpvec2

        idxMax = which(tmpvec2==1)
        idxMin = which(tmpvec2==0)
        tmpvec3 = rep(NA, length(dfTestGroups))

        if (length(dfTestGroups)==3)
        {
            idxMed = which(tmpvec2==0.5)
            MaxDiff = abs(tmpvec[idxMax] - tmpvec[idxMed])
            MinDiff = abs(tmpvec[idxMin] - tmpvec[idxMed])
            MeanDiff = mean(c(MaxDiff, MinDiff))
            tmpvec3[idxMed] = MeanDiff

        }else{
            MaxDiff = MinDiff = abs(diff(tmpvec))
        }
        tmpvec3[idxMax] = MaxDiff
        tmpvec3[idxMin] = MinDiff
        df_absdiff[gene, c(2:ncol(df_absdiff))] = tmpvec3
    }
    ## infer which group is the most similar to phenotype of interest in the test sets, and
    # which group is the most similar to the other phenotype in the training dataset (as the training dataset has only 2 phenotypes)
    # (and which group is intermediate, if there are 3 groups/clusters)
    distMatch = c()
    fracMatch = c()
    for (dfTestGrp in dfTestGroups)
    {
        distMatch = c(distMatch, sum(abs(df_updn[,1]-df_updn[,dfTestGrp])))
        fracMatch = c(fracMatch, sum(df_updn[,1]==df_updn[,dfTestGrp])/nrow(df_updn))
    }

    caseIdx = which(fracMatch == max(fracMatch))
    if (length(caseIdx) >1 ) # in the situation where there are 2 cluster with the same proportion of matches
    {
        absDiffInFracMatch = rep(NA, length(caseIdx))
        for (i in 1:length(caseIdx))
        {
            idxMatch = which(df_updn[,1]==df_updn[,dfTestGroups[caseIdx[i]]])
            absDiffInFracMatch[i] = abs(sum(df_absdiff[idxMatch, dfTestGroups[caseIdx[i]]]))
        }
        caseIdx = caseIdx[which(absDiffInFracMatch == max(absDiffInFracMatch))]
    }
    ctrlIdx = which(distMatch == max(distMatch[-c(caseIdx)]) )
    if (length(ctrlIdx) >1 )
    {
        ctrlIdx = ctrlIdx[which(ctrlIdx != caseIdx)]
    }
    dfout_umap$inferred_cluster_category = replace(dfout_umap$CLUSTER, dfout_umap$CLUSTER==caseIdx, 'case')
    dfout_umap$inferred_cluster_category = replace(dfout_umap$inferred_cluster_category, dfout_umap$CLUSTER==ctrlIdx, 'ctrl')
    if (length(dfTestGroups)==3)
    {
        intermIdx = c(1:3)[-c(caseIdx,ctrlIdx)]
        dfout_umap$inferred_cluster_category = replace(dfout_umap$inferred_cluster_category, dfout_umap$CLUSTER==intermIdx, 'interm')
    }
    write.table(dfout_umap, paste0("clusters/UMAP_CLUSTER_ASSIGNMENT_",FILE_NAME,'.txt'),
                sep='\t', col.names=T, row.names=F, quote=F)
    sink()
}
