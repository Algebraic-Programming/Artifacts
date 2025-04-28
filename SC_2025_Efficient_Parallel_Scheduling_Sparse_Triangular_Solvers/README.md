# SC25 Artefacts

Code and Artefacts for the SpTRSV submission. 

## License
See license files in respective folders

## How to run everything (on Linux)

Our experiments are conducted on three different architectures **Intel x86**, **ARM**, and **AMD x86** EPYC. See items 3 & 4 for differences. 

0. Prerequisites
   - Install OMP, see https://www.openmp.org/ for details.
      To install, e.g., on Ubuntu type
      ```bash
      sudo apt-get install libomp-dev
      ```
   - Install METIS, see http://github.com/KarypisLab/METIS for details.
      To install, e.g., on Ubuntu type
      ```bash
      sudo apt-get install libmetis-dev
      ```
   - Install Eigen, see https://eigen.tuxfamily.org/dox/ for details.
      To install, e.g., on Ubuntu type
      ```bash
      sudo apt-get install libeigen3-dev
      ```
   - Install Boost and its dependencies, see https://www.boost.org/ for details.
     To install, e.g., on Ubuntu type
     ```bash
     sudo apt-get install build-essential g++ python-dev autotools-dev libicu-dev libbz2-dev libboost-all-dev
     ```
1. Build OneStopParallel
   <!-- See README.md in one-stop-parallel folder -->
   ```bash
   bash build_one-stop-parallel.sh
   ```
2. Generate the data set
   1. Install and compile Sympiler aggregation for data set generation
      ```bash
      bash downloadSympilerForMatrixConversion.sh
      ```
   2. Generate the data set
      ```bash
      bash downloadDatasets.sh
      ```
   3. Analyse the data set
      ```bash
      bash generateGraphStatistics.sh
      ```
3. Build Sympiler aggregation for runtime evaluation. The scripts replace the Sympiler installation of Step 2.\
   If **Intel** or **AMD x86**:
   ```bash
   bash downloadSpMPSympilerForExperiments.sh
   ```
   If **ARM**:
   ```bash
   bash downloadARMSympilerForExperiments.sh
   ```
4.  Run the test suite on the respective architecture\
   (If the experiments take too long remove matrices from the data set and evaluate on the smaller data set) \
   On **Intel x86**
    ```bash
    cd one-stop-parallel
    bash run_sc25_experiment_Intel.sh
    cd ..
    bash run_sympiler_experiments_Intel.sh
    ```
    On **AMD x86**
    ```bash
    cd one-stop-parallel
    bash run_sc25_experiment_AMD.sh
    cd ..
    bash run_sympiler_experiments_AMD.sh
    ```
    On **ARM**
    ```bash
    cd one-stop-parallel
    bash run_sc25_experiment_ARM.sh
    cd ..
    bash run_sympiler_experiments_ARM.sh
    ```
5.  Run the evaluation scripts\
    To do the entire evaluation, the experiments need to be conducted on the three architectures, and results have to be put together in one place in the proposed folder structure.
    See README.md in the sc25-evaluation folder 

## How to run evaluation scripts on provided data

1.  Run the evaluation scripts\
    See README.md in the our_result_data_and_analysis folder
