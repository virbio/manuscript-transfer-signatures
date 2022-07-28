# function that runs the random forest model should be run for each training dataset
# where the inputs are:
# TRAINING_DIRECTORY - the directory where the outputs will be written (full path should be provided and the parent directory should already exist);
                       # should be specific to the training dataset;
# SAMPLES - the list of sample ids - must be present in the colum names of the gene count matrix
# LABELS - the categorical phenotypes of the samples (2 categories), must be provided in the same order than the corresponding samples
# MATCOUNTS - the normalized gene count matrix - where there is one row per gene and one column per sample; see end of script (commented out section)
              #how to generate MATCOUNTS object
# SIGNATURES - a character vector of signature names (where each signature name object contains a list of genes)

# where the output is:
# RF_feat_importance_<signatureName><SUFFIX>.txt  # matrix containing the feature importance of each gene for each left out sample and for the combined model
                                                  # where each row is a gene from the signature (<signatureName>) and each column a sample +
                                                  # the "combined_model" column
# RF_combined_model_<signatureName><SUFFIX>.RData # robject containing the random forest model; not used in the rest of the pipeline but in case of
                                                  # interest to dig into
# <label1>/RF_output_df_prob_<signatureName><SUFFIX>.txt # matrix with first column being the sample name,
                                                         # 2nd column the true label of the sample and
                                                         # 3rd column being the RF output score (from 0 to 1) of the sample being a <label1>;
                                                         # <label1> can be for example a case. The closer to 1 the more likely the sample is <label1>.
# <label2>/RF_output_df_prob_<signatureName><SUFFIX>.txt # matrix with first column being the sample name, 2nd column the true label of the sample and
                                                         # 3rd column being the RF output score (from 0 to 1) of the sample being a <label2>;
                                                         # <label2> can be for example a control. The closer to 1 the more likely the sample is <label2>.
                                                         # the sum of the 3rd columns of <label1>/RF_output_df_prob_<signatureName>.txt and
                                                         # <label2>/RF_output_df_prob_<signatureName>.txt is equal to 1.
                                                         # the <case>/RF_output_df_prob_<signatureName>.txt output can be used to extract the ROC and PR values.

# when running the function on random set of genes replace SUFFIX with: "_random1", for the first set of genes, "_random2", for the second set of genes, etc

RF_LOOCV = function(TRAINING_DIRECTORY, SAMPLES, LABELS, MATCOUNTS, SIGNATURES, SUFFIX='', NTREE=1000)
{
    library(rfUtilities)
    library(randomForest, quietly = TRUE )
    get = get0

    dir.create(paste0(TRAINING_DIRECTORY))
    setwd(paste0(TRAINING_DIRECTORY))

    dataFrameLabel = data.frame(true_labels=LABELS)
    rownames(dataFrameLabel) = SAMPLES

    for (lbl in sort(names(table(LABELS)))) # lbl stands fro label
    {
        assign(paste('dataFrameLabel_',lbl, sep=''), dataFrameLabel)
        dir.create(paste0(lbl))
    }

    str_eval=function(x) {return(eval(parse(text=x)))} # will be used to combine RF

    for (signatureName in SIGNATURES)
    {
        print (signatureName)
        genes = get(paste(signatureName)) # genes is a vector of ENSG gene names present in the signature <signatureName>
        dataFrameLabel[,paste(signatureName)]   = rep(NA, length(SAMPLES))
        for (lbl in sort(names(table(LABELS))))
        {
            tmp = get(paste('dataFrameLabel_',lbl, sep=''))
            tmp[,paste(signatureName)] = rep(NA, length(SAMPLES))
            assign(paste('dataFrameLabel_',lbl, sep=''), tmp)
        }

        dataFrame  = MATCOUNTS[genes,SAMPLES]

        RF_models = c()
        imp_features = as.data.frame(replace(dataFrame, !is.na(dataFrame), NA)) # create a dataframe for importance feature with same dimensions than dataFrame
                                                                                # but filled with NA
        for (k in 1:ncol(dataFrame))
        {
            #results for each run of LOOCV
            # code modified from https://github.com/jasonzhao0307/R_lib_jason/blob/master/RF_output.R
            #train and test data
            data_rf_train = as.data.frame(t(as.matrix(dataFrame[,-k])))
            data_rf_test  = as.data.frame(t(as.matrix(dataFrame[,k])))
            colnames(data_rf_test) = genes

            samples_train = SAMPLES[-k]
            samples_test  = SAMPLES[k]
            labels_train  = factor(LABELS[-k], levels=sort(names(table(LABELS))))
            labels_test   = factor(LABELS[k], levels=sort(names(table(LABELS))))

            data_rf_train$pheno = labels_train
            data_rf_test$pheno  = labels_test

            varNames = paste(genes, collapse = "+")
            rf.form  = as.formula(paste("pheno", varNames, sep = " ~ "))
            my_forest= randomForest(rf.form, data=data_rf_train, importance=TRUE, ntree=NTREE, keep.forest=TRUE)
            assign(paste('my_forest_', k, sep=''), my_forest)
            RF_models = c(RF_models, paste('my_forest_', k, sep=''))
            imp_features[,k] = importance(my_forest, type=2)[rownames(imp_features),]

            for (lbl in sort(names(table(LABELS))))
            {
                if (lbl %in% rownames(my_forest$confusion))
                {
                    tmp = get(paste('dataFrameLabel_',lbl, sep=''))
                    tmp[k,paste(signatureName)] = predict(my_forest, type="prob", data_rf_test)[,lbl]
                    assign(paste('dataFrameLabel_',lbl, sep=''), tmp)
                }
            }
        }
        final_RF_model = str_eval(paste("rf.combine(",paste(lapply(RF_models, function(x){noquote(x)}), collapse=','),")", sep=''))
        imp_features$combined_model = importance(final_RF_model, type=2)[rownames(imp_features),]

        write.table(cbind(data.frame(samples = rownames(imp_features)), imp_features),
                file=paste('RF_feat_importance_',signatureName, SUFFIX, '.txt', sep=''),
                sep='\t', row.names=F, col.names=T, quote=F, eol='\n')

        save(final_RF_model, file=paste('RF_combined_model_', signatureName, SUFFIX, '.RData',sep=''))

        for (lbl in sort(names(table(LABELS))))
        {
            write.table(cbind(data.frame(samples = rownames(get(paste('dataFrameLabel_',lbl, sep='')))),
                              get(paste('dataFrameLabel_',lbl, sep=''))[,c('true_labels',signatureName)]),
                file=paste(lbl, '/RF_output_df_prob_',signatureName, SUFFIX, '.txt', sep=''),
                sep='\t', row.names=F, col.names=T, quote=F, eol='\n')
        }
    }
}


##########################################################################################################################################
## generate MATCOUNTS input needed for RF_LOOCV()
#
##########################################################################################################################################
## IF GENE EXPRESSION MATRIX COMES FROM RNASEQ EXPERIMENT:
## load the raw counts gene expression matrix as count_matrix (where rows are genes/ENSG and columns are samples)
#
## keep only the genes that had 20 reads in at least 10 percent of the samples
#count_matrix = count_matrix[apply(gcmatrix, 1, function(x){quantile(x, 0.9)}) > 20,]
## remove the the dot and anything coming after the dot in the ENSG gene nomenclature
#rownames(count_matrix) = gsub("\\..*", "", rownames(count_matrix))
## compute the read per million
#rpm = t(t(count_matrix) / (colSums(count_matrix)/1e6))
## alternatively the million reads could have been computed before removing the noisy genes...
#MATCOUNTS = log10(rpm + 1e-7) # add a pseudocount and use the log10
#
#
##########################################################################################################################################
## IF GENE EXPRESSION MATRIX COMES FROM MICROARRAY EXPERIMENT (PRENORMALIZED):
## load the prenormalized signal per probe matrix as count_matrix (where rows are probes and columns are samples)
## load the probe to gene (ENSG) dictionary dataframe as probeDico (where the rownames of the dataframe are the probes id and the one and only column (names "ENSG") are the corresponding gene id)
#
#count_matrix_meanPerGene = aggregate(count_matrix[rownames(probeDico),], list(ENSG = probeDico$ENSG), FUN=mean, na.rm=TRUE)
#count_matrix_meanPerGene = data.frame(count_matrix_meanPerGene[,-1], row.names=count_matrix_meanPerGene[,1])
## keep only the genes that had signal in at least 1 of the samples
#count_matrix_meanPerGene = count_matrix_meanPerGene[rowSums(count_matrix_meanPerGene, na.rm=T)>0,]
#MATCOUNTS = log10(count_matrix_meanPerGene + 1e-7)


