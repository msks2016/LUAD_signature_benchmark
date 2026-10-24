# Reporting-completeness coding (76 papers)

Final coding used in the manuscript: `reporting_completeness_76papers_2026-10-12_v2.csv`
(`reporting_completeness_76papers.csv` is an identical copy). Numbers are computed by
scripts/41_reporting_completeness.R.

History of the coding files:
1. `reporting_completeness_76papers_autocoded.csv` - first pass; some items "unclear" or "none".
2. `reporting_completeness_76papers_2026-10-12.csv` - primary coder after checking cutoffs and external validation; this is the pre-discussion primary coding for the 10 dual-coded papers.
3. `reporting_completeness_76papers_2026-10-12_v2.csv` - scaling recoded with the definition agreed after dual coding ("implied" = RNA-seq mentioned without a unit; "not_stated" = no scale information).

Second coding:
- `dual_coding_moinul_original.csv` - second coder, 10 randomly drawn papers, all five items, before discussion; `dual_coding_filled.csv` and `dual_coding_merged.csv` - after discussion.
- `recheck_form_moinul_filled.csv` with `recheck_key.csv` - blinded second coding of the text-table item in the other 66 papers.

Bibliographic records:
- `records_76_pubmed.csv` - one row per paper (`paper_id` = screening record number): first author, year, title, journal, volume, pages, PMID and DOI, all taken from PubMed (checked 2026-10-24), and the benchmark signature where the paper is one of the 15. Script 41 merges it into Supplementary Table S19.
- In 9 papers the detailed extraction record differed from PubMed (paraphrased title in 4; wrong DOI in 4; wrong journal in 2; sixth author given as first author in 1); the `verification_note` column says what differed and how the paper was matched. Only the bibliographic fields were wrong; each record's title and content identify the paper, so the coding is not affected.

Definition of a text-table conflict: the paper's text (including the risk-score formula) and its tables disagree on a coefficient sign, a coefficient value (beyond rounding of the same value), or the number of genes or samples.
Any change to the final coding is recorded in its `notes` column with the reason.
