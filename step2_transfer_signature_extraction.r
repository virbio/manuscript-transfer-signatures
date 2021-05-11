# once the "step1_signature_evaluation_in_training_set.r" has been run with each literature signatures as well as 100 random signature of the same size,
# we can select the literature signatures that are above a given percentile threshold in terms of ROC AUC percentile and proceed to the following step

# function that extracts the transfer signatures should be run for each training dataset
# where the inputs are:
# TRAINING_DIRECTORY - the directory where the outputs from step1 have been written (full path should be provided and the parent directory should already exist)
                       # should be specific to the training dataset;
# PERCENTILE_THRESHOLD - the percentile threshold used (compared to the 100 random list of genes) to be selected to be part of the transfer signature
# SIGNATURES - a character vector of all literature signature names; importantly it should not contain the cognate signature for a give training dataset
# RANDOM_SUFFIX - a common string character suffix appended to the <signatureName> used for the output of the random list of genes, when running step1_signature_evaluation_in_training_set.r; f.ex. "_random" would be the common suffix of "_random1", "_random2", "_random3", etc
# COGNATE - a character vector of all the literature signature names that derived from the training dataset that is being assessed; is a subset of SIGNATURES
# PHENOTYPE_OF_INTEREST - a character string indicating the phenotype of interest, has to match either one of the categories provided as LABELS in step1.

# where the outputs are:
# transfer_signature_using_<PERCENTILE_THRESHOLD>'th_percentile.txt' # which is a list of genes in order of importance (most important on top)
                                                                     # making the transfer signature for the given training dataset
# literature_signatures_with_ROCAUCs_above_<PERCENTILE_THRESHOLD>th_percentile.txt # which contains a list of literature signature names that were above
                                                                                   # the percentile threshold


extract_transfer_signature = function(TRAINING_DIRECTORY, PERCENTILE_THRESHOLD=70, SIGNATURES, RANDOM_SUFFIX='_random', COGNATE, PHENOTYPE_OF_INTEREST)
{
    load("helper_functions.r")
    setwd(paste0(TRAINING_DIRECTORY))
    comb_gene_list = c() # comb_gene_list is a vector that contains the combined gene list
    comb_imp_scale = c() # comb_imp_scale is a vector that contains standardized gene importance feature from the signature that passed the PERCENTILE_THRESHOLD
    sig_list = c()       # sig_list is a vector that contains the list of signatures that are above the PERCENTILE_THRESHOLD

    for (signatureName in SIGNATURES)
    {
        if (!signatureName %in% COGNATE)
        {
            # extract a vector of ROC AUCs where at the first position we have the result obtained by the literature signature and
            # at position 2:101 we have the results obtained for the 100 random list of gene of the same size.
            rocaucs = extract_ROC_AUC(paste0(PHENOTYPE_OF_INTEREST,'/RF_output_df_prob_',signatureName,'.txt'))
            for (idx in c(1:100))
            {
                rocaucs = c(rocaucs, extract_ROC_AUC(paste0(PHENOTYPE_OF_INTEREST,'/RF_output_df_prob_',signatureName, RANDOM_SUFFIX, idx, '.txt')))
            }
            if (rocaucs[1] >= quantile(rocaucs,PERCENTILE_THRESHOLD/100))
            {
                featimp = read.table(paste0('RF_feat_importance_,'signatureName,'.txt'), header=T, rownames=1, sep='\t')
                comb_gene_list = c(comb_gene_list, rownames(featimp))
                comb_imp_scale = c(comb_imp_scale, scale(featimp$combined_model, center=T, scale=T))
                sig_list = c(sig_list, signatureName)
            }
        }
    }
    if (length(comb_gene_list) >0)
    {
        tmpordered = unique(comb_gene_list[order(comb_imp_scale, decreasing=T)])
        write.table(data.frame(tmpordered),
                    file=paste0('transfer_signature_using_',PERCENTILE_THRESHOLD, 'th_percentile.txt', sep=''),
                    quote=F, row.names=F, col.names=F)
        write.table(data.frame(sig_list),
                    file=paste0('literature_signatures_with_ROCAUCs_above_',PERCENTILE_THRESHOLD,'th_percentile.txt', sep=''),
                    quote=F, row.names=F, col.names=F)

    }
}



