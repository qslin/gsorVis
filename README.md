# genesets over representation functional enrichment analysis and visualizations #

The purpose of this repository is used to conduct functional enrichment analysis with GSOR and visualizations.

### repository structure ###

* 'R' folder includes all R scripts to conduct analysis with GSOR.
* 'run_functional_analysis_main.R' & 'run_summarizeEnrichGoRes_gsor.R' files: main files to call to run the GSOR and visualizations seperately.
* 'usage1_runGsor_GOsep_hpc.sh' & 'usage1_runGsor_kegg_local.sh': shell script examples calling 'run_functional_analysis_main.R' to conduct GSOR analysis on HPC and local PC respectively.
* 'usage2_gsor_summary_plot.R': dot-plot visualization of GSOR results.

### Analysis and visualization steps ###

* Step 1: run GSOR analysis 
  + `Rscript run_functional_analysis_main.R --help` instructs you what options are used for this GSOR analysis. In general, it includes below options:
      1. `--workingDirectory`: where the analysis should be conducted
      2. `--repoDirectory`: full path of this repository downloading locations
      3. `--inputFnames`: list of input file names seperated by ', ' (see examples), where includes the input gene lists
      4. `--inputNames`: corresponding names of above input gene lists saved in each input files with `--inputFnames`.
      5. `--inputType`: define types of input genes inside `--inputFnames`, such as 'geneListsInput'
      6. `--outputFname`: define where to save the output results
      7. `--genome`: which genome the input gene lists belong to, currently support hs for human, mm for mouse, rn for rat.
      8. `--genelistType`: define input gene list type, support 'GeneSymbol', 'ensembleID'.
      9. `--resFnameSuffix`: prefix of GSOR results, in case the same `--outputFname` used, results can still be saved with different file names.
      10. `--runGO` & `--goPool`: whether to run GSOR on GOs database with runing on BP, CC, and MF in combine `--goPool T` or seperately `--goPool F`.
      11. `--runKegg`: whether to run GSOR on KEGG database
      12. `--runReactome`: whether to run GSOR on Reactome database
      13. `--runNetworkCancerGenes`: whether to run GSOR on network cancer genes database
      14. `--runBroadMSigDB` & `--broadMSigDBCategory`: whether to run on Broad MSigDB database, if yes, use '--broadMSigDBCategory' to define which category to run.
  + Examples of running GSOR on HPC or local PC via shell script `usage1_runGsor_GOsep_hpc.sh` and `usage1_runGsor_kegg_local.sh` respectively.
  + GSOR analysis results are saved into the defined `--outputFname` subfolder defined prefix name in `--resFnameSuffix`.
  
* Step 2: dot-plot visualizations of GSOR results
  + The visulization example code: `usage2_gsor_summary_plot.R`, where 
     1. `bitbucketGitPath`: update the path where this repository is download 
     2. `enrichResDirName`: update where GSOR are saved, it is defined via `--outputFname` from above running
     3. `resFnamePrefix`: define as `--resFnameSuffix` from above running
     4. `toi.selection`: logical option to define whether to select certion terms of interes (toi), if yes, provide 1). `toi.sel`, a list of selected terms; 2) `precise.match`: whether to make precise match for `toi.sel`.
     5. `cluster.selection`: define whether to select certain clusters for visualizations, if yes, provide 1). `cluster.sel`: clusters names to be selected; 2). `cluster.sel.fnameSuffix`: file name prefix for dotplot visualization to save.
     6. `fdrOn`: logical option, whether to select enriched functions on FDR corrected p-value, if `fdrOn=F`, enriched functions are selected based on nominal p-values.
     7. `padjValListInput`: can be length of one for a unique FDR/p-value enriched functions selection, or a verctor of corresponding length of `cluster.selection` or all clusters to define different FDR/p-value cut-offs.
     8. `fulldotplotValInput`: list items corresponding to plot c(width, length), this depends on how many GSOR results saved inside `enrichResDirName`.
     9. `plotTop20`: whether to make a dot-plot of top 20 selected enriched functions dot-plot, if yes, please provide corresponding plot size via `topEnrichDotplotValInput`.
     10. `xLevelsReorder`: logical option, whether to order clusters presented on x-axix, if yes, please define the order via `xLevelsReorderVals`.
     11. Below are options if you would like to change, otherwise do not define them:
          + `axis.y.size`: if appears, can define y-axis font size, otherwise, by default = 20.
          + `axis.x.size`: if appears, can define x-axis font size, otherwise, by default = 20.
          + `legend.position`: define where to put legend, such as `legend.position = 'bottom'` or `legend.position = 'right'`.
          + `legend.size`: if appears, can define legend font size, by default = 14.
          + `legend.title.size`: if appears, can define legend title font size, by default = 15.
          + `ylimit`: if appears, can define the maximum values of FDR/p-value used for dot-plot.
          + `xlimit`: if appears, can define the maximum values of number of counts (dot radius size) for dot-plot.
     
