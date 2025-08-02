# glycobif
**Glycobif** is a tool for predicting carbohydrate utilization pathways encoded in bifidobacterial genomes.
Using amino acid FASTA files as input, it analyzes the representation of 589 curated metabolic functional roles (catabolic enzymes, transporters, and transcriptional regulators) and 68 curated carbohydrate utilization pathways.

A manuscript detailing:
* Metabolic reconstruction used in Glycobif
* Methods employed in the analysis
is available at [Arzamasov et al., 2025](https://doi.org/10.1038/s41564-025-02056-x)

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
     * Full pathway descriptions are available in `dictionary/phenotype_metadata_carbs.txt`
2. `PredictionsFull.xlsx` - a table representing the presence of functional roles in each pathway (1 = role is present, 0 = role is absent)
     * Full descriptions of functional roles are available in`dictionary/functional_roles.txt`
3. `tmp` - a folder containing temporary files, such as annotation files with locus tags of genes.

## Best practices and limitations
1. **Intended use**: glycobif is optimized for genomes of human-colonizing *Bifidobacterium* species. It can produce meaningful predictions for non-human *Bifidobacterium* genomes, though these may contain pathways not included in the curated set of 68 (e.g., pectin degradation). Predictions may also be generated for closely related genera such as *Alloscardovia*, but glycobif is not designed for broader application beyond these taxa
2. **Input genome type**: glycobif performs best with the genomes of cultured isolates. In these cases, predicted carbohydrate utilization pathways can be interpreted directly as predicted phenotypes (i.e., carbohydrate utilization capabilities)
3. Use with MAGs: glycobif can predict the pathway representation encoded in  metagenome-assembled genomes (MAGs), but several caveats apply:
     * MAGs with low completeness or high contamination (especially assembled from short-read data) can lead to false negatives and false positives, respectively. For reliable predictions, we recommend using MAGs with >97% completeness and <3% contamination
     * Short-read MAGs are particularly susceptible to assembly fragmentation, which can split gene clusters across multiple contigs and lead to loss of certain genes. A notable example is the HMO cluster I (H1) in *Bifidobacterium longum* subsp. *infantis*, where genes may appear to be missing due to assembly artifacts rather than true genomic variation. We recommend manually reviewing the presence and continuity of such clusters in MAGs of interest
     * MAGs may not represent individual strains. Consequently, pathway predictions may reflect the collective metabolic potential of the bin rather than a single organism
4. **Prediction accuracy and caveats**: While glycobif achieves high (>94%) overall phenotype prediction accuracy for isolates, benchmarking against *in vitro* growth data in [Arzamasov et al., 2025](https://doi.org/10.1038/s41564-025-02056-x) identified two sources of errors:
	- **False negatives** (growth despite predicted non-utilization), for certain monosaccharides (e.g., ribose) and human milk oligosaccharides (e.g., LNT), can occur when a complete catabolic pathway is present but the associated transporter is not detected, often due to gaps in current scientific knowledge. In cases where a perceived pathway incompleteness is observed, Glycobif assigns a binary phenotype of 0. If these pathways are of interest, users are encouraged to inspect the `PredictionsFull.xlsx` output to examine the representation of individual functional roles
	- **False positives** (predicted utilization but no growth) may result from disruptive mutations in genes encoding key pathway components. The current version of glycobif does not evaluate gene integrity or functionality at the nucleotide level

## Citation
If you used glycobif, please cite the following manuscript:
[Arzamasov et al., 2025](https://doi.org/10.1038/s41564-025-02056-x)
