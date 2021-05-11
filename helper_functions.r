# functions that extract ROC AUC and PR AUC
# where the input is:
# RF_output_score_matrix - character path to the <label>/RF_output_df_prob_<signatureName>.txt matrix output from "step1_signature_evaluation_in_training_set.r",
                         # where the first column is the sample name,
                         # 2nd column the true label of the samples and
                         # 3rd column is the RF output score (from 0 to 1) of the sample being a <label>;
                         # <label> can be for example a case. <label> has to be included in the path and has to be the first element before the "/".

# where the output is :
# either the ROC AUC or the PR AUC (as R object), it is not written to local directory

extract_ROC_AUC = function(RF_output_score_matrix)
    {
        RF_output_df = read.table(RF_output_score_matrix, header=T, rownames=1, sep='\t', stringsAsFactor=F)
        suppressWarnings(suppressMessages(library(AUC, quietly = TRUE)))
        detach("package:AUC")
        suppressWarnings(suppressMessages(library(pROC, quietly = TRUE)))
        # to make sure the auc() function will be the one from library pROC
        suppressWarnings(suppressMessages(library(ROCR, quietly = TRUE )))

        # LABELS is the true labels of the samples
        LABELS = RF_output_df[,1]
        # probs is a vector with the probabilities output per the RF (same order as LABELS
        probs = RF_output_df[,2]
        # pheno is the "category"  phenotype we are interested in; ex: 'case'
        pheno = gsub('\\/.*','', RF_output_df)

        # rename the labels into binary for LABELS (0 if it is not the correct label and 1 if it is the label of interest)
        binary_label = rep(NA, length(LABELS))
        binary_label = replace(binary_label, LABELS != pheno, 0)
        binary_label = replace(binary_label, LABELS == pheno, 1)
        binary_label = ordered(binary_label, levels=c(0,1))

        newpredicted = as.numeric(probs)

        roc_obj = roc(binary_label, newpredicted, quiet=TRUE, direction="<", algorithm=1)

        tmpdf = data.frame(sensitivity=roc_obj$sensitivities,
                           specificity=roc_obj$specificities)
        #tmpdf$AUC = auc(roc_obj)

        return(auc(roc_obj))
    }


extract_PR_AUC <- function(RF_output_score_matrix)
    {
        RF_output_df = read.table(RF_output_score_matrix, header=T, rownames=1, sep='\t', stringsAsFactor=F)
        suppressWarnings(suppressMessages(library(ROCR, quietly = TRUE )))
        suppressWarnings(suppressMessages(library(AUC, quietly = TRUE )))

        # LABELS is the true labels of the samples
        LABELS = RF_output_df[,1]
        # probs is a vector with the probabilities output per the RF (same order as LABELS
        probs = RF_output_df[,2]
        # pheno is the "category"  phenotype we are interested in; ex: 'case'
        pheno = gsub('\\/.*','', RF_output_df)

        # rename the labels into binary for LABELS (0 if it is not the correct label and 1 if it is the label of interest)
        binary_label = rep(NA, length(LABELS))
        binary_label = replace(binary_label, LABELS != pheno, 0)
        binary_label = replace(binary_label, LABELS == pheno, 1)

        binary_label = ordered(binary_label, levels=c(0,1))
        newpredicted = as.numeric(probs)

        pred = prediction(newpredicted, binary_label)
        PR = performance(pred, "prec", "rec" )
        # transform the PR object so that it is compatible with the method auc from the AUC package
        newPR = list(cutoffs=unlist(PR@x.values), measure=unlist(PR@y.values))
        # actually the y axis is precision... but it's just to be compatible with their normenclature
        # if there is NaN in the first position let's replace it appropriately
        if (is.na(newPR$measure[1]))
        {
            tri = cbind(as.numeric(as.vector(binary_label)), newpredicted)
            tri = tri[order(tri[,2], decreasing=T),]
            newPR$measure[1] = tri[1,1]
        }

        accu=accuracy(newpredicted, binary_label, perc.rank = TRUE)
        class(newPR)=class(accu)

        tmpdf = data.frame(precision=newPR$measure,
                           recall=newPR$cutoffs)
        #tmpdf$AUC = auc(newPR, min = 0, max = 1)

        return(auc(newPR, min = 0, max = 1))
    }


