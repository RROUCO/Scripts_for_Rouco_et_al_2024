To analyze C-HiC sequencing data the following steps were performed:

1st step:Create config file
The configuration file for HiCup was set up for each sample as shown in the example provided here `1_HiCup_Config_file.conf` 
The correct path for each FASTQ pair was provided. All the other parameters were maintained across samples. 

2nd step:load modules
Once the configuration file was created and prior to running the sbatch the modules indicated in the file `2_Modules_to_load.txt` were loaded.

3rd step: Run sbatch HiCup
A sbatch file `3_HiCup_sbatch.sh` to run HiCup was created for each sample indicating the corresponding configuration file (`1_HiCup_Config_file.conf`) for each sample.

4th step: Run sbatch Juicer normalization
A sbatch file `4_Juicer_generate_normalized_maps_sbatch.sh` to run Juicer was created for each sample. 
For all the samples the same chromosome coordinates (mm9;chr3:65000001-68500000) and binsze (5000bp and 10000bp) were used.
Only hicup.bam file path was modified between samples to provide the matching one.
Fragmented mm9 genome file `DpnII_starts_mm9.txt` provided in the `4_Juicer_generate_normalized_maps_sbatch.sh` was generated in R using the code described in `create_fragment_file_mm9.R`

5th step: Virtual Capture-C 
To run the Virtual Capture-C using the sbatch file `5a_V-CaptureC_command_line.sh` a sbatch file `5b_Generate_V-CaptureC_profiles_example.sh` was provided for each sample.
The `5b_Generate_V-CaptureC_profiles_example.sh` was modified for each sample providing the path for the hicup.bam
For all the samples the same chromosome coordinates (mm9;chr3:65000001-68500000) were used.
Coordinates for the viewpoint(s) of interest were provided in the file `5c_Virtual_viewpoints_5kb.txt`. 
In particular for Supplementary Figure 7 a 5kb viewpoint region, covering the Shox2 promoter, was used (mm9;chr3:66781632-66786630)

