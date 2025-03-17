# glycobif
**Glycobif** is a tool for predicting carbohydrate utilization pathways in bifidobacterial genomes.
Using amino acid FASTA files as input, it analyzes the representation of 589 curated metabolic functional roles (catabolic enzymes, transporters, and transcriptional regulators) and 68 curated carbohydrate utilization pathways.

A preprint detailing:
* Metabolic reconstruction used in Glycobif
* Methods employed in the analysis
is available at [Arzamasov et al., 2024](https://doi.org/10.1101/2024.07.06.602360)

[Semen Leyn](https://github.com/sleyn/) and Marat Kazanov created the initial pipeline.

## Dependencies and installation
Glycobif requires the following dependencies:
 * [DIAMOND](https://github.com/bbuchfink/diamond) (tested v2.1.4 and v2.1.11)
 * Python libraries: `pandas`, `numpy`, `tqdm`, `scikit-learn`, `scipy`
 * R libraries: `data.table`, `openxlsx`, `caret`, `yaml`, `ranger`, `tidyverse`

**To create a mamba/conda environment with the required dependencies, download the repository and run:**
```
mamba env create -f env/glycobif.yml
```

**To test the pipeline, run:**
```
bash run_glycobif.sh
```
This will run the pipeline for *Bifidobacterium longum* subsp. *infantis* ATCC 15697 and output the results to the `output/` directory.

**Note:** If your mamba/conda environment has a name different from glycobif, manually specify it in `run_glycobif.sh` before running the script.

**Note2:** If you experience R package version errors after installation, reinstall `caret` through R in the created environment

## Usage
Glycobif uses two inputs:
* Amino acid FAA files
* Tab-separated file (genome list) with two columns: (i) genome_IDs (names of FAA files without the extension) and (ii) taxonomy 
```
genome_ID	curated_taxonomy
ATCC15697	Bifidobacterium longum subsp. infantis ATCC 15697 = JCM 1222
```

**Basic usage**
```
bash run_glycobif.sh [-t num_threads] [-i input_directory_faa] [-l input_list] [-o output_directory]
```
Parameters:
* `-t` - number of CPUs to use; default `1`
* `-i` - path to input directory; default `input/faa`
* `-l` - path to genome list file; default `input/genome_list.txt`
* `-o` - output directory;  default `output`

## Output
1. `BPM.xlsx` - a Binary Phenotype Matrix (BPM) that contains the binary representation of 68 carbohydrate utilization pathways (1 = complete pathway is present, 0 = pathway is absent or incomplete)
     * Full pathway descriptions are available in `dictionary\phenotype_metadata_carbs.txt`
2. `PredictionsFull.xlsx` - a table representing the presence of functional roles in each pathway (1 = role is present, 0 = role is absent)
     * Full descriptions of functional roles are available in`dictionary\functional_roles.txt`
3.`tmp` - a folder containing temporary files, such as annotation files with locus tags of genes.

## Best practices and limitations
1. Glycobif is designed for use with the genomes of human-colonizing *Bifidobacterium* species. While it will produce meaningful output for non-human *Bifidobacterium* species, they may contain pathways not represented in the curated set of 68 pathways (e.g., pectin degradation). Glycobif may also work with genomes from closely related genera like *Alloscardovia*, but it is not designed for use with other bacterial taxa
2. Glycobif works best with the genomes of cultured isolates as input. In this case, the presence of carbohydrate utilization pathways can be directly interpreted as predicted carbohydrate utilization capabilities (phenotypes)
3. Glycobif also works with metagenome-assembled genomes (MAGs), but there are several caveats:
     * Using MAGs with low completeness and high contamination may lead to false negative and false positive pathway predictions, respectively. We recommend using MAGs with completeness >97% and contamination <3% for more reliable results
     * MAGs may not correspond to individual bacterial strains. As a result, the predicted carbohydrate utilization pathways represent the metabolic potential of the community rather than a single organism
4. While the overall prediction accuracy is high (>94%), comparison with *in vitro* growth data in [Arzamasov et al., 2024](https://doi.org/10.1101/2024.07.06.602360) revealed the presence of false negative and false positive predictions:
	- False negative predictions for monosaccharides (e.g., ribose) and certain human milk oligosaccharides (LNT) can occur due to incomplete pathway predictions, resulting from limited knowledge about respective glycan transporters (i.e., when a catabolic pathway is present, but known transporter is absent, glycobif assigns 0 instead of 1). If you are particularly interested in these pathways, please check the output file `PredictionsFull.xlsx` for the genomes of interest and consider manual curation
	- False positive predictions may arise due to mutations in genes encoding functional roles. Glycobif, in its current implementation, does not account for such mutations

## Citation
If you used glycobif, please cite the following manuscript:
[Arzamasov et al., 2024](https://doi.org/10.1101/2024.07.06.602360)
