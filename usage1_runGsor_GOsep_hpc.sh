#!/bin/bash 

#PBS -l qos=biocore
#PBS -l nodes=1:ppn=1
#PBS -l walltime=10:00:00
#PBS -o /gpfs/data/biocore-analysis/yli/CRI-BIO-757-MED-ELengyel-yli-lszhu/pre-post_menopausal_analysis/FT_pre_post_menopausal_integration_results/functional_enrichment_results_preManuscript/run_DEGs_menstrual_proliferativeSecretoryComp_GOsep.sh.out
#PBS -e /gpfs/data/biocore-analysis/yli/CRI-BIO-757-MED-ELengyel-yli-lszhu/pre-post_menopausal_analysis/FT_pre_post_menopausal_integration_results/functional_enrichment_results_preManuscript/run_DEGs_menstrual_proliferativeSecretoryComp_GOsep.sh.err
#PBS -l mem=8gb 

module load biocore_pipelines/bioinfoDT-v1

echo hostname = $HOSTNAME

echo `date`
echo --------- 
workDir=/gpfs/data/biocore-analysis/yli/CRI-BIO-757-MED-ELengyel-yli-lszhu/pre-post_menopausal_analysis/FT_pre_post_menopausal_integration_results/functional_enrichment_results_preManuscript
repoDir=/group/bioinformatics/yli2/gsor_functional_enrichment_analysis/

##Rscript $repoDir/run_functional_analysis_main.R -h
##---
cd $workDir

Rscript $repoDir/run_functional_analysis_main.R --workingDirectory $workDir --repoDirectory $repoDir --inputFnames '/gpfs/data/biocore-analysis/yli/CRI-BIO-757-MED-ELengyel-yli-lszhu/pre-post_menopausal_analysis/preOnly_manuscript_plots/results_wNewAnnotation_DEGs/AFcomb_menstrual_proliferativeSecretoryComp_cellsWoSE1SE3preBP2EN2_UQnorm_bmiRaceCov/expCondCompDeMarkers_Proliferative-Secretory_adjSig_up_wNewAnnotation_15SelClusters_clusterST1.txt, /gpfs/data/biocore-analysis/yli/CRI-BIO-757-MED-ELengyel-yli-lszhu/pre-post_menopausal_analysis/preOnly_manuscript_plots/results_wNewAnnotation_DEGs/AFcomb_menstrual_proliferativeSecretoryComp_cellsWoSE1SE3preBP2EN2_UQnorm_bmiRaceCov/expCondCompDeMarkers_Proliferative-Secretory_adjSig_up_wNewAnnotation_15SelClusters_clusterST2.txt, /gpfs/data/biocore-analysis/yli/CRI-BIO-757-MED-ELengyel-yli-lszhu/pre-post_menopausal_analysis/preOnly_manuscript_plots/results_wNewAnnotation_DEGs/AFcomb_menstrual_proliferativeSecretoryComp_cellsWoSE1SE3preBP2EN2_UQnorm_bmiRaceCov/expCondCompDeMarkers_Proliferative-Secretory_adjSig_up_wNewAnnotation_15SelClusters_clusterST3.txt, /gpfs/data/biocore-analysis/yli/CRI-BIO-757-MED-ELengyel-yli-lszhu/pre-post_menopausal_analysis/preOnly_manuscript_plots/results_wNewAnnotation_DEGs/AFcomb_menstrual_proliferativeSecretoryComp_cellsWoSE1SE3preBP2EN2_UQnorm_bmiRaceCov/expCondCompDeMarkers_Proliferative-Secretory_adjSig_down_wNewAnnotation_15SelClusters_clusterST1.txt, /gpfs/data/biocore-analysis/yli/CRI-BIO-757-MED-ELengyel-yli-lszhu/pre-post_menopausal_analysis/preOnly_manuscript_plots/results_wNewAnnotation_DEGs/AFcomb_menstrual_proliferativeSecretoryComp_cellsWoSE1SE3preBP2EN2_UQnorm_bmiRaceCov/expCondCompDeMarkers_Proliferative-Secretory_adjSig_down_wNewAnnotation_15SelClusters_clusterST2.txt, /gpfs/data/biocore-analysis/yli/CRI-BIO-757-MED-ELengyel-yli-lszhu/pre-post_menopausal_analysis/preOnly_manuscript_plots/results_wNewAnnotation_DEGs/AFcomb_menstrual_proliferativeSecretoryComp_cellsWoSE1SE3preBP2EN2_UQnorm_bmiRaceCov/expCondCompDeMarkers_Proliferative-Secretory_adjSig_down_wNewAnnotation_15SelClusters_clusterST3.txt' --inputNames 'ST1_proliferative, ST2_proliferative, ST3_proliferative, ST1_secretory, ST2_secretory, ST3_secretory' --inputType geneListsInput --outputFname 'menstrual_proliferativeSecComp_DEGs_STs_gsor' --genome hs --genelistType GeneSymbol --resFnameSuffix 'menstrual_proSecComp_STs_DEGs' --runGO T --goPool F --runKegg F --runReactome F --runHumanDisease F --runNetworkCancerGenes F --runBroadMSigDB F --broadMSigDBCategory C7  > run_functional_analysis_menstrual_proliferativeSecretoryComp_STs_DEGs_GOsep.log 2>&1


echo ---------
echo `date` 

