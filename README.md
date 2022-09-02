# mrflow
*A package for expedited Mendelian Randomization analysis of population data.*

## Functions
* `mr_analyze` quickly imports and executes 
## Future Implementation
* A standard MR data prep workflow, including one-shot execution of:
	* Interconversion between gene symbols, synonyms, ENSGs, Entrez, etc
	* SNP retrieval from database
	* Automatic pruning for LD--perhaps worth exploring writing own LD algo since it won't require any sort of web interfacing
	* Conduct Two-Sample MR w/ summary statistic SNPs ([`mr_analyze()`](https://github.com/aszatrowski/mrflow/blob/master/R/mr_analyze.R))
		* Fix "unknown or uninitialized column `trait`" error
	* Quick and efficient power analysis 

## Dependencies
`mrflow` currently depends on the following packages:
* [`MendelianRandomization`](https://cran.r-project.org/web/packages/MendelianRandomization/) for statistical analysis ([ref](https://doi.org/10.1093/ije/dyx034))
* [`TwoSampleMR`](https://github.com/MRCIEU/TwoSampleMR) for retrieval for data from the [IEU openGWAS database](https://gwas.mrcieu.ac.uk/) ([ref](https://doi.org/10.1101/2020.08.10.244293))
* [`mRnd`](https://github.com/kn3in/mRnd) for power analysis ([ref](https://doi.org/10.1093/ije/dyt179))