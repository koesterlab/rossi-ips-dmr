### Summary of GermLayerTracker Paper and Comparison to Our Analysis

**Goal of the Paper:**
The study aimed to identify, which germ layer human pluripotent stem cells (PSCs) are likely to differentiate into during early differentiation. They developed a tool called GermLayerTracker, which calculates a pluripotency score based on DNA methylation (DNAm) at specific CpGs, indicative of differentiation potential.

---

### Key Findings

1. **Comparison to classic marker genes:**

   * Traditional markers like SOX17, GATA6 (endoderm), GATA2, TBXT (mesoderm), SOX2, PAX6 (ectoderm) were not differentially expressed in our dataset.
   * This highlights the **added value of DNAm-based markers** over conventional transcriptional markers, especially for early differentiation.

2. **Scatterplots (Figure 1B and S2A/S2B):**

   * Show mean beta values (DNAm levels) for pairwise comparisons of PSCs vs germ layers.
   * Only CpGs with a **methylation difference ≥ 0.2** and **adjusted p-value ≤ 0.05** are highlighted.
   * Color coding appears as follows (interpretation based on your notes):

     * **Gray:** ΔDNAm < 0.2 or not significant
     * **Blue/Red/Black:** Significant hyper- or hypomethylation toward one germ layer
     * When both comparisons are significant, color assignment is unclear; currently, the analysis **excludes CpGs where both comparisons are non-significant**.
   * Observed trend: Endoderm and mesoderm share many DNAm changes, particularly in dataset 2.

3. **Correlation between DNAm and gene expression:**

   * Promoter hypomethylation generally correlates with **upregulated gene expression**, and hypermethylation with **downregulation** (Figure S2D).
   * Our data partially supports this, especially for **POU5F1**, although many genes differ between analyses.
   * Hypomethylation is defined as **lower DNAm in the germ layer vs PSC**.

---

### GermLayerTracker Scoring System

1. **PSC Score (Figure 2):**

   * Three top CpGs for PSCs:

     * cg00661673 (PALLD)
     * cg00933813 (intergenic)
     * cg21699252 (MYCNOS)
   * For each CpG, DNAm values are combined into a pluripotency score:

     * **Hypomethylated CpGs** in PSCs: score = 1 − DNAm
     * **Hypermethylated CpGs**: score = DNAm
   * This results in a **high score for PSCs** and **low score for differentiated cells**.
   * In our dataset, these three CpGs are **lowly methylated in PSCs** and highly methylated in differentiated layers.

     * If we sum DNAm directly, our score is inverted compared to the paper (PSC appears low).
     * Correct calculation: we should **invert the DNAm for hypomethylated PSC CpGs**, as in the paper.

2. **Germ-layer-specific CpGs:**

   * Endoderm (ENDO): cg20548013 (PHACTR1), cg14521421 (DENND2B), cg08913523
   * Mesoderm (MESO): cg14708360, cg08826152 (ADORA2B), cg11599718 (VPS37B)
   * Ectoderm (ECTO): cg01907071 (THSD4), cg18118164 (EFNA5), cg13075942 (RAD51B)
   * Combined endoderm-mesoderm (ENDOMESO): cg23385847 (CAMK4), cg24919344, cg11147278

3. **Germ-layer differentiation scores (Figure 3B):**

   * For **hypomethylated sites**, the score is **1 − DNAm**, to ensure it **increases with differentiation**.
   * Scores are calculated relative to the **mean DNAm of undifferentiated PSCs**.
   * This allows differentiation toward specific germ layers to be quantified quantitatively.

---

### Key Points for Our Analysis

* Focus on **individual CpGs** rather than DMRs.
* Use **PSC-specific hypomethylated CpGs** and apply the **1 − DNAm transformation** to ensure our pluripotency score aligns with the original paper.
* Scatterplots can be generated similarly, highlighting CpGs with **ΔDNAm ≥ 0.2** and **adjusted p-value ≤ 0.05**.
* Gene expression validation is challenging; DNAm changes may be more sensitive for **early differentiation events**.
* Endoderm and mesoderm share many changes—consider combining them for a more robust differentiation score.

---

I can also make a **visual table/diagram** showing how hypomethylation/hypermethylation and the scoring work, which makes the comparison much easier to explain.

Do you want me to do that next?
