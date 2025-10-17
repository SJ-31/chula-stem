# TCR Analysis Routine

Taking BD Rhapsody full-length TCR/BCR data as input, this analysis routine does the following:

+ Clonotype calling and clustering
    + The definitions follow that of [scirpy](https://scirpy.scverse.org/en/latest/generated/scirpy.tl.define_clonotypes.html)
        + clonotype/clone: a collection of T or B cells with identical CDR3 nucleotide sequences on both the VJ and VDJ chains
        + clonotype cluster: clonotypes with similar CDR3 amino acid sequences on both VJ and VDJ chains, with similarity defined using the `tcrdist` metric ([link to paper](https://www.nature.com/articles/nature22383))
+ Attempts to look up the top clones of each sample in external reference databases. Currently, two methods are included, which differ by the similarity metric used for comparison. In addition, not all of the databases can be used with each method.
    + `scirpy`: uses `tcrdist`
    + `tcrmatch`: uses a k-mer matching approach, described [here](https://www.frontiersin.org/journals/immunology/articles/10.3389/fimmu.2021.640725/full)
+ Predicts pMHC binding for the top clones in each sample with [MixTCRpred](https://github.com/GfellerLab/MixTCRpred). This is a deep learning-based approach which creates models for individual TCR-pMHC interactions. The routine uses a configurable selection of the authors' pre-trained models, which are described [here](https://github.com/GfellerLab/MixTCRpred/blob/main/pretrained_models/info_models.csv)
+ Generates various plots to aid interpretation, described in the next section
    + **Note:** all ranks in plots and preprocessing are sorted in decreasing size e.g., the rank 1 clonotype is the largest (has the most cells)

## Output details
+ `mixtcrpred_results.csv`: aggregated results of MixTCRpred's predictions of TCR-pMHC binding for the top clones in each sample. This file contains model metadata including the MHC class, peptide, and its origin.
    + Prediction quality is measured using two metrics:
        1. **score:** the higher the score, the more likely the sequence is to bind the pMHC complex
        2. **%rank:** "the fraction of TCR sequences with higher scores multiplied by 100. Low %rank score sequences are more likely to bind. Most true binders have been observed to have a % rank of less than 0.5–0.1"
+ `reference_comparison`: directory of sub-directories corresponding to each of the lookup methods listed previously
    + `scirpy`: results of querying [IEDB](https://www.iedb.org/), [VDJdb](https://vdjdb.cdr3.net/) and [McPAS](https://friedmanlab.weizmann.ac.il/McPAS-TCR/) with `scirpy`
    + `tcrmatch`: results of querying IEDB and [CEDAR](https://cedar.iedb.org/) databases with `tcrmatch`
        + `tcrmatch_invalid.csv` contains clones that could not be compared using the `tcrmatch` metric. This metric relies on the CDR3 beta chain and thus won't work if the chain is missing (`reason == no_cdr3_beta`) or if the receptor subtype is TRG+TRD (`reason == wrong_subtype`)
+ `sequences`: directory containing full-length TCR sequences for the top clones in each sample. Each clone should have at least two entries (VJ and VDJ), but there may be more if the clones possess different V, J, D sequences.
    + FASTA headers are annotated with the V, D, J, C calls (if available), chain subtype, rank in the sample, as well as clone size
    + Each sample has two fasta files, one for nucleotide sequences (`<sample_name>.fasta`), the other for amino acid sequences (`<sample_name>-aa.fasta`)
    + `ignored_sequences.csv` contains cells of top clones for which a sequence could not be extracted. The `reason` column takes one of four values
        1. **incomplete:** sequence alignment does not span the entire V(D)J region, the inverse of `complete_vdj` field in the [AIRR Rearrangement Schema](https://docs.airr-community.org/en/latest/datarep/rearrangements.html)
        2. **is_ambiguous:** cell's chain pairing is non-standard e.g. TRA+TRD, TRA+IGH
        3. **no_sequence:** sequence is shorter than a user-specified minimum length
        4. **no_ir:** cell has no immune receptor chain
+ `combined.h5mu`: MuData file aggregating all the AIRR data from each sample

### `plots` Directory
+ `clone_ranks.png`, `clone_cluster_ranks.png`: cohort-level plots showing the proportion of the sample taken up by the top-ranked clones and clone clusters, respectively
+ `public_private.png`, `public_private_clusters.png`: plots visualizing the top clones and clone-clusters in the cohort, and whether any are shared between samples
+ `clone_expansion`: directory containing sample-level plots showing the relationship between clone size and cell type
+ `cell_type_abundance`: directory containing sample-level plots showing the distribution of clonotypes across cell types.
+ `vdj`: directory containing annotated visualizations of full-length TCR sequences for the top clones, organized by sample and chain type.
+ `alpha_diversity.png`: cohort-level comparison of clonotype alpha diversity. The included metrics range from 0–1, and the higher they are, the more diverse the sample
    + `richness`: the number of unique clonotypes in the sample, normalized by the total number of unique clonotypes in the cohort
    + `pielou_e`: Pielou's Evenness index, higher when many clones are present in the sample and clone sizes are evenly distributed
    + `simpson`: Gini-Simpson index, which can be interpreted as the probability that two randomly selected cells are different clonotypes
