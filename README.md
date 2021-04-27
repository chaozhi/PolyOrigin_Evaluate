# Test PolyOrigin

## Evaluate with real 3x3 half-diallel potato data

* In "run_polyorigin", "step1_run_polyorigin.jl" perform haplotype reconstruction using input genofile "TableS2_dose.csv" and pedfile "TableS3_ped.csv".

* In "run_mappoly"
  * The 3x3 half-diallel crosses data: "TableS2_dose.csv" is split by sub-populations and chromosomes "potato_pop?_chr??-dose_mappoly.csv", and "TableS3_ped.csv" is split into "potato_pop?_ped.csv"
  * "step1_run_mappoly.R" perform haplotype reconstruction using MAPpoly for each sub-population and each chromosome. 

## Compare with MAPpoly and TetraOrigin for simulated F1 data

Run the following script files in the folder "simF1":

* "step1_simtest_f1data.R" produce simulated data using PedigreeSimR and saved in the empty folder F1data. See F1data.zip for the example of simulated data

* "step2-1_run_by_polyorigin.jl" perform haplotype reconstruction using PolyOrigin. Results will be saved in the empty folder "res_polyorigin"

* "step2-2_run_by_tetraorigin.nb" perform haplotype reconstruction using TetraOrigin. Results will be saved in the empty folder "res_tetraorigin"

* "step2-3-1_prepare_input_mappoly.R" transform the simulated data files in the input files of mappoly.  Results will be saved in the empty folder "res_mappoly"

* "step2-3-1_prepare_input_mappoly.R" transform the simulated data files in the input files of mappoly.  Results will be saved in the empty folder "res_mappoly"

* "step2-3-2_run_by_mappoly_dr0.R" (for simulated data without double reduction) and "step2-3-2_run_by_mappoly_dr0.5.R" ( with double reduction) perform haplotype reconstruction using MAPpoly. Results will be saved in the empty folder "res_mappoly"

* "step3-2_calculate accuracy_tetraorigin.nb" and "step3-3_calculate accuracy_mappoly.jl" caculate accuracy for the results of tetraorigin and mappoly, based on ground true-value files in F1data. The accuracies for PolyOrigin results have been calculated in the step2-1.  
