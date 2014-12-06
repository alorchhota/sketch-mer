Scripts:
========
We can categorize the scripts inside the working directory into 3 groups. 
1) Method implementation
    - countLeastSquares.py
    - KmerSupplier.py
2) Test data generation
    - exact_kmer_count.py 
    - generate_test_data.py
3) Run tests.
    - readUtil.py
    - generate_countminsketch_output.py
    - evaluate_accuracies.py


Directory Structure:
====================
To correctly run this program, we need a specific folder structure and also file names should follow a specific naming pattern. The working directory should contain two folders - 1) data, and 2) results.

data: All data are saved in "data" folder in the working directory. Three types of data are stored here.
    - Sequencing reads: Stored directly under 'data' folder in fasta file format. Filenames follow the pattern: "<DATASET>.fa".
    - Exact counts: We saved all unique 22-mers and their counts in "22mer_exact_counts" folder. Here filenames follow the pattern: "<DATASET>_counts.txt". We counted this files using a script, "exact_kmer_count.py".
    - Test data: We saved test datasets into "test" folder. Filenames follow the pattern: "test<DATASET>_<BATCHSIZE>.txt". We generated these files using a script, "generate_test_dataset.py".

results: Unless otherwise mentioned, we store results of any script into 'results' folder in the working directory. It must have a sub-directory "tmp", which we will use for generating exact k-mer counts.

Note: We created necessary directories and put necessary files in appropriate directories.

How to run main scripts:
========================
Given, test data are generated and put in correct directory, you can run our codes following the steps below.

Step-0: Go to the working directory (inside "sketch-mer" foldr) in terminal.
Step-1: Run the following command to generate estimated counts and save them into files.

python generate_countminsketch_output.py

This script will produce several output files into the "results" folder. File name patterns and their contents are described below:

    - pattern "cmresult_<DATASET>_<BATCHSIZE>.txt": contains estimated counts produced by countminsketch.
    - pattern "lsresult_<DATASET>_<BATCHSIZE>.txt": contains estimated counts produced by sketch-mer.

Note: This script may take several minutes to finish, as it runs all seven datasets, two of which take relatively long time.

Step-2: Run the following command to calculate accuracies and produce plots.

python evaluate_accuracies.py

This script will also produce several output files into the "results" folder. File name patterns and their contents are described below:
    
    - accuracies.txt: contains errors produced by countminsketch and sketch-mer for different datasets (column: dataset) and different batch sizes (column: batchsize). Note: the last two columns titled "countmin_mean_err", and "lsquare_mean_err" contains average error produced by countminsketch and sketch-mer, respectively. You may skip other different columns, as we produced those columns just for explorign purpose.
    - pattern "bar_<DATASET>.png": plot to show effects of batch size in <DATASET>.
    - pattern "hist_<DATASET>_<BATCHSIZE>.png": plot to show error distribution using batch size of <BATCHSIZE> in <DATASET>.
    - bar_err_comp_<BATCHSIZE>.png: plot to compare accuracies of all datasets using batch size of <BATCHSIZE>.
    - lineErrVsNumKmer_<BATCHSIZE>.png: plot to show effects of the number of distinc k-mers when batch size is <BATCHSIZE>.
    - lineErrVsHashSize_200.png: plot to show effects of hash size (controlled by epsilong parameter in countminsketch). Note this can be calculated after running the whole procedure for different values of epsilon. Apology for using hard-coded values in this function due to time constraint.

Note: Plots sometimes may not look pretty. We can tweak some parameters to make it pretty.


How to run helper scripts:
==========================
We wrote some helper scripts to generate test data. We used those scripts before we actually run our program. We kept outputs produced by these scripts into correct folder. However, if you want to generate those files, you can do so by following the steps below.

Step-1: set appropriate value of the variable 'dataset' at the end of the file "exact_kmer_count.py", and then run the following command.

python exact_kmer_count.py

The ouput is given into the file: "results/<DATASET>_counts.txt". Note: you have to put this file into the 'data/22mer_exact_counts/' directory to actually use it in the program.

Step-2: Run the following command to generate test datasets. 

python generate_test_dataset.py

It will generate test files for different datasets and batch sizes in the "results" folder. Filenames will follow the pattern "test_<DATASET>_<BATCHSIZE>.txt". Note: you have to copy these files into the "data/test/" directory to actually use them in the program.


