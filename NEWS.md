# popgeninfer 0.0.0.9000 (2023-05-08)

- Added function `af_test`, which implements a clean and well-tested version of the allele frequency test first implemented by Tiffany Tu for the NS GWAS project.  This tests if the data from two datasets have the same allele frequency (null hypothesis) or not (alternative hypothesis).  Supports combining data from different ancestries into a single p-value.

# popgeninfer 0.0.1.9000 (2023-05-09)

- Function `af_test` changed default behavior for vector and matrix inputs:
  - Previously, vector inputs were treated as data to combine, and matrices were flattened to vectors
  - Now, vectors are treated as column matrices, and for matrices in general rows are treated as separate tests and columns as data to combine
  - The new setup aligns better with the actual uses of this function in various lab projects.
- Added package vignette, confirming with simple simulations that `af_test` p-values are uniform under the null and concentrated near zero under the alternative hypothesis.

# popgeninfer 0.0.2.9000 (2024-05-15)

- Added function `af_test_or`, which complements `af_test` to calculate odds ratios, together with confidence intervals, and p-values.  However, this test has a different alternative model (more constrained because there is a single OR shared across ancestries) so the p-value are different when there is more than one column/ancestry.

# popgeninfer 0.0.3.9000 (2024-05-16)

- Function `af_test_or` corrected a minor bug that swapped the roles of the allele counts (which are response) and study (which is covariate).  Also internally optimized how this Binomial data is passed to function `glm`, and made minor clarifications to documentation
- Updated vignette to include `af_test_or`, and extended simulations to include new cases relevant to that function.  The vignette validates the accuracy of the new function and contrasts it to `af_test`.

# popgeninfer 0.0.4.9000 (2025-02-05)

- Added functions `bias_geno`, `revcomp`, and `flip_revcomps`, to simulate platform-specific genotyping biases and ambiguous strand flips.  All are based on originals from Tiffany Tu's project.

# popgeninfer 0.0.5.9000 (2025-02-11)

- Added functions `af_filter` and `filter_category` that implement the basic filter based on `af_test` considering flips, resulting in the three categories "remove", "flip", and "keep" depending on p-value threshold.
- Vignette increased sample size and introduced other minor corrections to avoid fatal plotting problems (involving log scales in particular).

# popgeninfer 0.0.6.9000 (2025-02-13)

- Added functions `binary_eval` and `filter_eval` to calculate precision and recall for binary and 3-category (flip, remove, keep) cases, respectively.

# popgeninfer 0.0.7.9000 (2025-03-21)

- Added function `gwas_eval`, and added it to the `filter_eval` wrapper too.

# popgeninfer 0.0.8.9000 (2025-03-31)

- Function `filter_eval` now returns a new column "removed" that simply counts total removals.

# platformbias 0.0.9.9000 (2026-08-18)

- Renamed package from `popgeninfer` to `platformbias`.
