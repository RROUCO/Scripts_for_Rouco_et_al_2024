To analyze C-HiC sequencing data the following steps were performed:

1st step:Create config file
The configuration file for HiCup was set up for each sample as shown in the example provided here `1_Config_file_example.conf` 
The correct path for each FASTQ pair was provided. All the other parameters were maintained across samples. 

2nd step:load modules
Once the configuration file was created and prior to running the sbatch the modules indicated in the file `2_Modules_to_load.txt" were loaded`.

3rd step: Run sbatch HiCup
A sbatch file to run HiCup was created for each sample indicating the corresponding configuration file for each sample.

4th step: Run sbatch Juicer normalization
A sbatch file to run Juicer was created for each sample. 
For all the samples the same chromosome coordinates (mm9;chr3:65000001-68500000) and binsze (5000bp and 10000bp) were used.
Only bam file path was modified between samples to provide the matching one.


