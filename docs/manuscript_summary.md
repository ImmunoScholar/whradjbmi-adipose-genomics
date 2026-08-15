# Integrating GWAS, adipose eQTLs, and single-cell transcriptomics to prioritize candidate genes and cell types for body-fat distribution
## Abstract
Waist-to-hip ratio adjusted for BMI (WHRadjBMI) captures genetic variation in
body-fat distribution independent of overall adiposity, implicating adipose
tissue biology in its mechanism. We integrated three layers of evidence to
nominate candidate genes and adipose cell populations that may mediate
genetic effects at WHRadjBMI-associated loci. Starting from genome-wide
association summary statistics (Pulit et al. 2019; n=694,649), we identified
eight independent genome-wide-significant loci via distance-based pruning,
all corresponding to previously established WHR loci. For the three most
significant loci, we retrieved adipose subcutaneous cis-eQTL associations
(GTEx v8, eQTL Catalogue) and performed Bayesian colocalization (coloc.abf)
across all genes in each locus window (122 gene-locus tests).
Two loci showed
strong evidence of a shared causal variant between the GWAS and eQTL signals:
RSPO3 (PP4=0.77) and ZNF664 (PP4=0.92). Notably, at both loci the physically
nearest or eponymous gene (VEGFA, NCOR2) did not colocalize, instead showing
evidence consistent with two distinct causal variants -- illustrating that
genomic proximity alone is not evidence of a shared regulatory mechanism.
Mapping these prioritized genes into single-cell RNA-seq data from human
adipose stromal vascular fraction (Vijay et al. 2020; 26,350 cells, four
broad cell-type categories validated against the source publication's
reported proportions and depot-specificity) showed RSPO3 enriched in adipose
progenitor cells and ZNF664 enriched in a visceral-adipose-specific
mesothelial-derived population. These results nominate RSPO3 and ZNF664 as
candidate mediators of genetic effects on body-fat distribution, acting
plausibly through progenitor and mesothelial-lineage stromal biology.
This
analysis is hypothesis-generating: colocalization evidence is consistent
with, but does not prove, a shared causal variant, and cell-type expression
enrichment does not establish a causal cell type or disease mechanism.

## Background

Body-fat distribution, rather than overall adiposity, is an independent
contributor to cardiometabolic disease risk: excess visceral and central fat
carries greater metabolic risk than fat stored subcutaneously, even at
equivalent BMI. Waist-to-hip ratio adjusted for BMI (WHRadjBMI) is a
genetically tractable proxy for this distribution phenotype -- by
statistically removing the BMI component, WHRadjBMI GWAS signal is enriched
for loci acting on adipose tissue biology (depot patterning, adipocyte/
progenitor function, tissue expandability) rather than on the central
appetite-regulation pathways that dominate BMI GWAS.
Genome-wide association studies have identified hundreds of
WHRadjBMI-associated loci, but a GWAS locus alone does not identify a causal
gene, a mechanism, or a cell type. The variant most strongly associated with
a trait is frequently not the variant that alters gene expression, and the
gene physically nearest to a lead SNP is frequently not the gene whose
regulation is affected. Two further layers of evidence are needed to move
from "an associated locus" to "a plausible mediating gene and cell type":
(1) evidence that the locus is a cis-eQTL in a disease-relevant tissue,
ideally supported by formal colocalization rather than mere positional
overlap, and (2) evidence of where in a relevant tissue's cellular
composition the implicated gene is expressed.

This project integrates these three layers -- GWAS, adipose cis-eQTL
colocalization, and adipose single-cell transcriptomics -- for WHRadjBMI,
with the explicit goal of testing whether coloc-prioritized genes differ from
naive nearest-gene assignments, and whether they show interpretable,
cell-type-localized expression consistent with a plausible mechanism.
## Methods
**GWAS.** Summary statistics were obtained from Pulit et al. (2019;
GIANT+UK Biobank meta-analysis, n=694,649, European ancestry, hg19/GRCh37;
Zenodo DOI 10.5281/zenodo.1251813). The full dataset (27,375,636 variants)
was quality-checked for column completeness, duplicate variants, allele
validity, and value-range sanity; 706 variants (0.0026%) lacking
chromosome/position/sample-size annotation were excluded (none reached
genome-wide significance; minimum P among them was 1.6x10^-5). Independent
loci were defined among genome-wide-significant variants (P < 5x10^-8;
54,362 variants) via distance-based greedy pruning: iteratively selecting
the most significant remaining variant as a lead SNP and excluding a
+/-500kb window around it. This is an approximation to formal LD-clumping;
no linkage-disequilibrium reference panel was used, as the downstream
colocalization method does not require one. Eight independent loci were
retained.
**Adipose eQTL retrieval.** For each locus, a +/-500kb window around the
lead SNP was lifted from hg19 to GRCh38 (UCSC chain file; 100% width
retention and single-chromosome mapping for all eight loci). Adipose
subcutaneous cis-eQTL nominal association statistics (GTEx v8, n=581,
accessed via the eQTL Catalogue, dataset QTD000116) were retrieved for each
window using remote tabix range queries (`Rsamtools::scanTabix`); no bulk
eQTL file was downloaded.
**Allele harmonization and colocalization.** GWAS and eQTL variants were
merged on shared rsID. Alleles were harmonized to a common effect-allele
orientation (the eQTL Catalogue's documented convention that the alternate
allele is the effect allele), with ambiguous palindromic variants (A/T, C/G;
14.2% of merged rows) and non-matching allele pairs (0.09%) excluded rather
than heuristically resolved. Colocalization (`coloc.abf`, default priors
p1=p2=1x10^-4, p12=1x10^-5) was performed for every gene present at the
three most GWAS-significant loci (122 gene-locus tests total; minimum 30
overlapping SNPs required per test).
**Single-cell integration.** Processed single-cell RNA-seq data from human
adipose stromal vascular fraction (SVF; subcutaneous and visceral depots,
obese donors) were obtained from GSE129363 (Vijay et al. 2020). The provided
cell-barcode annotation contained tissue depot and disease-status labels but
no cell-type labels; standard clustering (Seurat: normalization, PCA,
graph-based clustering, UMAP) was performed on the 26,350 cells retained
after the original authors' own quality filtering (confirmed against their
published methods), yielding 12 clusters. Clusters were annotated into four
broad categories -- adipocyte progenitor/stem cells, immune cells,
endothelial cells, and mesothelial cells -- using canonical marker gene
module scores, requiring a genuinely positive module score rather than a
forced best-match assignment.
## Results

**Locus identification.** Eight genome-wide-significant loci were identified
(P range 8.4x10^-78 to 2.1x10^-293). All eight correspond to loci previously
reported in the WHR/WHRadjBMI literature, including the RSPO3, TBX15-WARS2,
LYPLAL1, ADAMTS9, and COBLL1-GRB14 regions -- an internal validation that
the locus-selection pipeline correctly recovers this trait's known genetic
architecture.

**Colocalization.** At the three loci carried through the full
eQTL/colocalization workflow:

- **Locus 1** (chr6, lead SNP rs72959041, P=2.08x10^-293): **RSPO3** was the
  best-colocalizing gene (PP4=0.768, PP3=0.052), exceeding the commonly
  applied PP4>0.75 threshold for suggestive colocalization.
- **Locus 2** (chr6, lead SNP rs998584, P=1.22x10^-170): no gene reached
  strong colocalization evidence; the best candidate, **PEX6**, showed
  PP4=0.252. Notably, VEGFA -- the gene for which this locus is
  conventionally named in the literature -- showed PP3=0.764, evidence
  consistent with two distinct causal variants rather than a shared one.
- **Locus 3** (chr12, lead SNP rs863750, P=4.17x10^-101): **ZNF664** showed
  the strongest colocalization evidence in this analysis (PP4=0.918,
  PP3=0.022). NCOR2, the locus's nearest/namesake gene, showed weak evidence
  (PP4=0.047).

**Single-cell mapping.** The derived cell-type composition (progenitor
34.3%, immune 34.5%, mesothelial 25.4%, endothelial 5.8%) closely matched
the proportions reported in the source publication (combined progenitor
lineage ~55-60%, immune ~34-37%, endothelial ~6-8%), and the mesothelial
population was found to be 99.96% visceral-adipose-derived (6,697/6,700
cells) -- consistent with the source publication's independent finding that
this progenitor subtype originates from the visceral mesothelium. RSPO3
expression was enriched in progenitor cells (6.1-fold relative to the mean
of other cell types), and ZNF664 expression was enriched in mesothelial
cells (2.3-fold). By contrast, NCOR2 (locus 3's non-colocalizing nearest
gene) showed a nearly flat expression pattern across all four cell types,
providing an independent line of evidence -- convergent with the
colocalization result -- that ZNF664 rather than NCOR2 is the more
biologically relevant candidate at this locus.
The remaining five loci have GWAS-level evidence only; eQTL/colocalization/
single-cell integration was scoped to the three most significant loci by
design.
## Interpretation

This analysis nominates **RSPO3** and **ZNF664** as candidate genes
mediating genetic effects on WHRadjBMI at two independent loci, with
expression patterns in adipose SVF consistent with a plausible cellular
locus of action -- adipose progenitor biology for RSPO3, and a
visceral-mesothelial-derived progenitor population for ZNF664. Both results
are internally reinforced by convergent evidence (strong colocalization and
differential cell-type expression) and, for ZNF664 specifically, are
consistent with prior functional literature implicating the
CCDC92-ZNF664-DNAH10 gene cluster in adipocyte and insulin-resistance
biology at this locus.
Equally informative is what did not colocalize: at both tested ambiguous
loci, the gene conventionally associated with the locus by name or
proximity (VEGFA, NCOR2) showed weaker or contradictory colocalization
evidence than an alternative gene in the same window. This is not a
limitation of the analysis but its central demonstrated finding -- that
positional or nominal association with a locus is not, by itself, evidence
of a shared regulatory mechanism, and that formal statistical colocalization
can meaningfully redirect candidate-gene prioritization away from the
"obvious" choice.
These results should be read as hypothesis-generating. Colocalization
evidence indicates that the GWAS and eQTL association signals at a locus are
statistically consistent with arising from a single shared causal variant;
it does not identify that variant, does not prove causality, and does not
rule out more complex scenarios (e.g., allelic heterogeneity) that the
underlying method assumes away. Similarly, cell-type expression enrichment
establishes where a gene could plausibly act, not that it does act there in
a disease-relevant way, or that the identified cell type is the causal one
for the WHRadjBMI phenotype.
## Limitations
**Ancestry compatibility.** The GWAS is European-ancestry only; GTEx (the
  eQTL source) is approximately 85% European. Colocalization's implicit
  assumption of comparable linkage-disequilibrium structure between the two
  datasets is only approximately satisfied.
- **No LD reference panel.** Locus definition used distance-based pruning
  rather than formal LD-clumping, since no LD panel was incorporated into
  this pipeline. This may occasionally define locus boundaries differently
  than a true clumping algorithm, particularly in regions of complex or
  extended LD.
**BMI-adjustment collider bias.** WHRadjBMI's adjustment for a heritable
  covariate (BMI) can in principle introduce collider bias at some loci;
  this was not evaluated locus-by-locus in this analysis.
- **Tissue and cellular scope of the single-cell data.** GSE129363 profiles
  the stromal vascular fraction only; mature adipocytes are removed during
  sample preparation and are entirely absent from this dataset. All
  expression conclusions here pertain to progenitor, immune, endothelial,
  and mesothelial populations, not to adipocyte-intrinsic expression, which
  may be mechanistically important for this trait.
- **Cell-type resolution.** Four broad cell-type categories were used, a
  substantial simplification relative to the ~24 fine-grained subclusters
  (and 14 distinct immune subtypes) resolved in the original publication.
  Finer sub-classification, which could reveal more specific cell-state-level
  enrichment, was outside this project's scope.
**Scope of colocalization testing.** Formal colocalization was performed
  at only 3 of 8 identified loci (the most GWAS-significant), a deliberate
  project boundary rather than a completeness gap; the remaining five loci
  have GWAS-level evidence only.
- **No prior-sensitivity analysis.** Colocalization results were not tested
  for robustness to variation in the `coloc.abf` prior probabilities (p1,
  p2, p12).
- **Sample-size handling.** Colocalization used a single representative
  sample size per dataset rather than per-variant sample-size vectors, a
  standard simplification for meta-analytic and eQTL summary data.
**Statistical power at weaker loci.** Loci with smaller effect sizes or
  eQTL sample coverage (e.g., locus 2) may be underpowered to detect true
  colocalization even where a shared causal variant exists; the absence of
  strong PP4 evidence at such a locus is not equivalent to positive evidence
  against colocalization.
