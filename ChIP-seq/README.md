To analyze ChIP-Seq data the following pipelines were adapted from https://github.com/lldelisle/myNGSanalysis/ and used depending on the type of dataset

For datasets from GEO that we wanted to reanalyze (Andrey et al., 2017; Sheth et al., 2016) we used the sbach command line `1a_Narrow_peak_MACS2_GEOdowloaded_ChIP_SE.sh` providing the FASTQ files through an external table like the one provided in `1b_Narrow_peak_MACS2_GEOdowloaded_table_example.txt`
This pipeline uses a Narrow Peak approach from MACS2 and triming Nextseq adaptors

For datasets produced in this study from FACS sorted cells we used the sbach command line `2a_Broad_peak_MACS2_FACS_sorted_cells.sh` providing the FASTQ files through an external table like the one provided in `2b_Broad_peak_MACS2_FACS_sorted_cells_table_example.txt`
This pipeline uses a Broad Peak approach from MACS2 and triming TrueSeq adaptors

Early, common, late putative enhancers classification on FACS sorted maintaining forelimb datasets we run the pipeline described in `pipeline_putative_early_late_enhancers_FACS_sorted_maintaining_forelimb.sh`


