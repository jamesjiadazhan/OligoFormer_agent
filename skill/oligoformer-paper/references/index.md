# Paper navigation

Use exact headings to locate current line numbers. Read relevant passages, not both complete documents. This index locates evidence; it does not summarize findings.

## Main paper — [paper.md](paper.md)

| Exact heading | Look here for |
| --- | --- |
| Abstract | One-paragraph summary of aims, method and headline results |
| 1 Introduction | RNAi/siRNA background and prior siRNA efficacy predictors |
| 2.1.1 Dataset collection | Nine source datasets, Huesken/Mixset/Takayuki grouping, normalization and redundancy removal; Table 1 |
| 2.1.2 Flanking region of target mRNAs | Choice of 19-nt flanking regions (57-nt mRNA input) |
| 2.1.3 Sequence normalization | Truncation/padding of siRNA and mRNA sequences |
| 2.1.4 RNA-FM model embedding | How RNA-FM pretrained embeddings are used |
| 2.1.5 Thermodynamic parameters | The 24 thermodynamic features |
| 2.2 The architecture of OligoFormer | Model architecture (Conv2D, BiLSTM, Oligo encoder, MLP); Figure 1 |
| 2.3 Cross-validation methods | Intra-dataset, inter-dataset and cross cell-line evaluation setups and metrics |
| 3.1 Performance on intra-dataset validation | Intra-dataset comparison with baseline methods |
| 3.2 Performance on inter-dataset validation | Huesken-to-Mixset results; Table 2 with AUC/PRC/F1/PCC of all methods |
| 3.3 Ablation study | Feature and module ablation results |
| 3.4 Model interpretation | Attention/feature interpretation and sequence motifs |
| 3.5 OligoFormer with off-target searching for siRNA design | Off-target (PITA, TargetScan), toxicity and functionality filtering pipeline |
| 4 Discussion | Limitations and future directions |
| Data availability | Code and data locations |

For long sections, search a narrower subsection or prompt. Headings inside fenced quotations are source content, not document section boundaries.

## Supplementary information — [supplement.md](supplement.md)

| Exact heading | Look here for |
| --- | --- |
| 1.1 Data preprocessing | Redundancy removal, efficacy normalization, flanking-length analysis |
| 1.2 Thermodynamic parameters | Equations for nucleotide indicators, ΔG and ΔH with worked example |
| 1.3 RNA-FM embedding | RNA-FM embedding details |
| 1.4 Oligo encoder embedding | Equations of the encoder layers and dimensions |
| 1.5 Prediction calculation of OligoFormer | How the final efficacy score is computed |
| 1.6 Details of model training and test | Training hyperparameters, evaluation metric definitions, cross-validation |
| 1.7 Comparison methods | Baselines: Monopoli-RF, OligoWalk, siRNAPred, i-Score, s-Biopredsi, DSIR |
| 1.8 Ablation study | Feature and module ablation setup |
| 1.9 Off-target | Functionality criteria, toxicity (seed viability) and off-target search procedure |
| 1.10 Mismatch Handling Module | Mismatch kon model |
| 1.11 User-friendly OligoFormer based on Docker | Docker usage commands |
| 2 Supplementary Figures | Supplementary Figs. 1-6 and captions |
| Supplementary Table 1.5. Cross validation | Cross cell-line and inter-dataset validation metrics |
| Supplementary Table 4. Ablation of OligoFormer components (check marks indicate included components per configuration column) | Ablation metrics per component configuration |

For long sections, search a narrower subsection or prompt. Headings inside fenced quotations are source content, not document section boundaries.

## Assets

Figure and table captions in the documents link to the files below. Open only the needed image; for a table, read its header and relevant rows first.

- `assets/figure/`: main figures (JPEG).
- `assets/supp_figs/`: supplementary and extended-data figures (JPEG).
- `assets/table/`: main tables (CSV, or JPEG when transcription is unreliable).
- `assets/supp_table/`: supplementary tables (CSV, or JPEG fallback).

Asset paths are relative to the skill root. CSVs retain internal blank rows; captions and merged headers may also occupy rows. Consult the document notes before treating every row as data.
