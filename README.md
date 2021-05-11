## Date: 17AUG2020
## Author: Julia di Iulio

the scripts provided match the study design provided in Figure S1 of the manuscript.
the scripts provided in the folder are:

- step1_signature_evaluation_in_training_set.r (matching Figure S1A workflow)
- step2_transfer_signature_extraction.r (matching Figure S1B workflow)
- step3_prediction_in_test_set.r (matching Figure S1C workflow)
- helper_functions.r


step1_signature_evaluation_in_training_set.r contains a function (RF_LOOCV) to run random forest models.
helper_functions.r contains functions (extract_ROC_AUC and extract_PR_AUC) to extract ROC and PR AUCs.
step2_transfer_signature_extraction.r contains a function (extract_transfer_signature) to extract transfer signatures.
step3_prediction_in_test_set.r contains a function (umap_cluster) that performs, unsupervised clustering followed by cluster label inference.

inputs needed to run each function are described at the beginning of the respective scripts.
outputs generated by each function are described at the beginning of the respective scripts.