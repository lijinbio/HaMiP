# CMSIP: Hydroxymethylation anlaysis of CMS-IP data

A scalable, accurate, and efficient solution for hydroxymethylation analysis of CMS-IP sequencing data.

![Workflow of CMSIP.](https://github.com/lijinbio/cmsip/raw/master/cmsip_flowchart.png)

## Installation

CMSIP has been deployed in Bioconda at https://anaconda.org/bioconda/cmsip. It is encouraged to install CMSIP from Bioconda due to most runtime dependencies will be installed automatically. The following channels should be added in Conda. Namely,

```
conda config --add channels defaults
conda config --add channels bioconda
conda config --add channels conda-forge
conda install cmsip
```

Alternatively, CMSIP has been also deployed in PyPI at https://pypi.org/project/cmsip, and it can be installed via `pip`.

```
pip3 install cmsip
```

In some cases, users want to build CMSIP manually from source code at https://github.com/lijinbio/cmsip. Below is an example installation steps.

```
git clone https://github.com/lijinbio/cmsip.git
cd cmsip
python3 setup.py install
```

In order to run CMSIP after a manual installation, the following dependent software are required.

| Software | URL |
|-------|-------|
| Python 3 | https://www.python.org |
| Matplotlib | https://matplotlib.org |
| PyYAML | https://pyyaml.org |
| bedtools | https://bedtools.readthedocs.io |
| R software | https://www.r-project.org |
| R package DESeq2 | https://bioconductor.org/packages/release/bioc/html/DESeq2.html |
| R package genefilter | https://bioconductor.org/packages/release/bioc/html/genefilter.html |
| R package RVAideMemoire | https://cran.r-project.org/web/packages/RVAideMemoire/index.html |
| Gawk | https://www.gnu.org/software/gawk |
| MOABS | https://github.com/sunnyisgalaxy/moabs |

## Documentation

CMSIP takes in a configuration file for input data and program parameters. CMSIP can be run end-to-end, starting from raw FASTQ files to peak calling and differential hydroxymethylation identification. One can also start the pipeline from intermediate steps. For example, using alignment files as input so that mapping steps will be skipped.

### Inspection of configuration

The configuration file is in a YAML format. Two example templates are `config_fastq.yaml` and `config_bam.yaml` under https://github.com/lijinbio/cmsip/blob/master/config. `config_fastq.yaml` is used as a full CMSIP running from FASTQ inputs, while `config_bam.yaml` is adapted to input existing BAM files so that CMSIP will skip the long-time alignment step. The inspection of configuration is explained below.

1. sampleinfo

The sampleinfo section defines metadata information in analysis. Below metadata information can be specified.

| Parameter | Description |
|-------|-------|
| sampleinfo.sampleid | the unique identifier to one sample |
| sampleinfo.group | the biological group of the sample, e.g., KO or WT |
| sampleinfo.filenames | the absolute path of raw FASTQ files |
| sampleinfo.reference | the absolute path of the reference BAM file when `aligninfo.inputbam` is True |
| sampleinfo.spikein | the absolute path of the spike-in BAM file when `aligninfo.inputbam` and `aligninfo.usespikein` is True |

 2. groupinfo

This section defines biological comparison `group1` - `group2`, e.g., KO - WT.

| Parameter | Description |
|-------|-------|
| groupinfo.group1 | the first group in biological comparison |
| groupinfo.group2 | the second group in biological comparison |

3. resultdir

This directory is default working directory storing intemediate result files, such as BAM and BED files.

4. aligninfo

This section specifies parameters used in raw reads alignment.

| Parameter | Description |
|-------|-------|
| aligninfo.inputbam | True for BAM inputs. Default: FASTQ inputs. |
| aligninfo.reference | FASTA file of the reference genome, e.g. `hg38.fa`. |
| aligninfo.usespikein | True for spike-in libraries, otherwise False. This option controls the normalization method, either a spike-in normalization using spike-in mapping, or reduced to WIG sum in reference genome. |
| aligninfo.spikein | FASTA file of the spike-in genome, e.g. `mm10.fa`. |
| aligninfo.statfile | the output statistics file. This file includes quality control statistics as well as estimated normalization factors. |
| aligninfo.barplotinfo | a barplot of normalized WIG sums of samples. |
| aligninfo.numthreads | number of threads in alignment program. |
| aligninfo.verbose | Print verbose message |

5. genomescaninfo

This section defines parameters for CMS measurement construction.

| Parameter | Description |
|-------|-------|
| genomescaninfo.readextension | True to extend reads length before CMS measurement construction. |
| genomescaninfo.fragsize | the fixed fragment size to extend when readextension is True. |
| genomescaninfo.windowfile | an intermediate window file with fixed-size genomic regions. |
| genomescaninfo.referencename | the UCSC genome name to fetch reference genome size. E.g., hg38 or mm10. |
| genomescaninfo.windowsize | the window size |
| genomescaninfo.readscount | CMS measurement using readcount (True) or mean WIG (False). |
| genomescaninfo.counttablefile | the result count table file. |
| genomescaninfo.verbose | Print verbose message |

6. dhmrinfo

Parameters in this section is for DMR detection.

| Parameter | Description |
|-------|-------|
| dhmrinfo.method | The statistical method used in DHMR detection. Available methods: `ttest`, `chisq`, `gtest`, `nbtest`, `nbtest_sf`. `ttest` is calling Student's t-test to examine the mean difference of CMS measurements between two biological groups. `chisq` and `gtest` are Pearson’s Chi-squared and G-test to test if sums of CMS measurements fit the numbers of replicates between two biological condtions. `nbtest` applies negative binomial generalized linear model to formulate CMS measurements, and Wald test evaluates the significance of logarithmic fold change. By default, CMS measurement are adjusted by size factors using spike-in normalization. In `nbtest_sf`, CMS measurements are normalized by the median-ratio algorithm (previously used in DESeq2 for transcriptome measurements). |
| dhmrinfo.meandepth | Average depth to filter out low-depth windows. This step is essential to save computing resources and increase power of downstream statistical inference |
| dhmrinfo.testfile | The result file with statistical outputs for whole genome windows |
| dhmrinfo.qthr | q-value threshod for DHMW. |
| dhmrinfo.maxdistance | Maximum distance to merge adjacent DHMWs into DHMRs |
| dhmrinfo.dhmrfile | The final DHMR result file after merging adjacent DHMWs. |
| dhmrinfo.numthreads | The number of threads. |
| dhmrinfo.nsplit | The number of split of windows. This option controls parallelization with `dhmrinfo.numthreads`. |
| dhmrinfo.verbose | Print verbose message. |
| dhmrinfo.keepNA | Keep genome windows ruled out by independent filtering. |

7. useinput

To indicate if the input data is used during CMS-IP sequencing.

8. inputinfo

If `useinput` is True, this section is required to specify input data. When input data is used, peak windows are identified first by comparing CMS measurements between group 1/2 and their input data. Then, the union of peak windows are tested for DHMR between group 1 and group 2.

| Parameter | Description |
|-------|-------|
| inputinfo.group1 | The label for the first group input data. |
| inputinfo.group2 | The label for the second group input data. Group 1 and group 2 can share same set of input data. |
| inputinfo.method | The statistical method used in peak calling. See `dhmrinfo.method`. |
| inputinfo.qthr | q-value threshold for peak calling. |
| inputinfo.testfile1 | Statistical test results for group 1 peaking calling. |
| inputinfo.dhmrfile1 | Peak regions for group 1. |
| inputinfo.testfile2 | Statistical test results for group 2 peaking calling. |
| inputinfo.dhmrfile2 | Peak regions for group 2. |
| inputinfo.inputfilterfile | Union of peak regions in group 1 and group 2. |
| inputinfo.verbose | Print verbose message. |

### A toy example using BAM inputs

To facilitate the running of CMSIP, a toy example is generated using existing BAM inputs. The example is accessible at https://github.com/lijinbio/cmsip/blob/master/example. The `example` directory consists of running scripts and example BAM files. Below commands will generate the configuration file and run the example.

```
$ ./config.sh ## Generate config.yaml
$ ./fasta.sh ## download the reference genome and the spike-in genome under ./fasta
$ ./run.sh ## run the example
```

1. `config.sh`

This script will generate the running configuration file. The inspection of configuration file has been explained above. This example includes small BAM files for 2 KO and 2 WT samples, together with 3 input samples. Spike-in BAM files are also included for spike-in normalization. These BAM files are under the `./bamfile` directory. The `gtest` is used for peaking calling and DHMR detection.

2. `fasta.sh`

This script is to download required FASTA file for reference genome and spike-in genome. These FASTA files are used in MCALL for bisulfite conversion ratio (BCR) estimation. FASTA files are downloaded into a local directory `./fasta`.

3. `run.sh`

The simple command to run CMSIP:

```
$ cmsip -c config.yaml
```

Intermediate and results files are stored under `./outdir`.

## Contact

Maintainer: Jin Li, lijin.abc@gmail.com.
PI: De-Qiang Sun, dsun@tamu.edu.

