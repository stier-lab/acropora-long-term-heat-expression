# RNA-seq

144 margin libraries (Days 1, 3 and 10; both temperatures and wound states; patches A, C and D) are being processed at UC Davis by Rachael Bay and Shreya. No expression counts have arrived yet; when they do, they go in `data/raw/sequencing/`.

| Need | File |
|---|---|
| Analysis plan, written before seeing the counts | `expression_integration_analysis_plan.md` |
| Candidate genes from the literature | `candidate_genes_reference.csv` |
| Candidate-gene web page (built from the CSV) | [`LTH_candidate_genes.html`](https://stier-lab.github.io/acropora-long-term-heat-expression/team_summary/LTH_candidate_genes.html); rebuild with `python3 docs/rnaseq/build_candidate_genes_page.py` |
| Per-library covariates: treatment, tank, day, symbiont density, host cluster, symbiont read shares | `output/tables/43_final_snp_rnaseq_covariates.csv` |
| Notes on the genetic and symbiont data | `docs/DATA_NOTES.md` |

The expression samples are destructively sampled fragments, not the 48 fragments followed to Day 15, and they come from before any fragment budded. Expression can test group differences (heat, wounding, patch), but it can't be matched to an individual fragment's regrowth.
