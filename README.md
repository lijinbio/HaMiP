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

In order to run CMSIP after a manual installation, the following dependent softwares are required.

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

### Example configuration file

Two example configurations are `template_cms.yaml` and `template_cms_bam.yaml` under https://github.com/lijinbio/cmsip.

## Contact

Maintainer: Jin Li, lijin.abc@gmail.com.
PI: De-Qiang Sun, dsun@tamu.edu.

