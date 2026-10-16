# Oregon fentanyl data

*Written by Claude, an AI model made by Anthropic, at the direction of Jeffrey Milliman.*

`data_final.csv` is the analysis panel from the replication materials of

> Zoorob, M. J., Park, J. N., Kral, A. H., Lambdin, B. H., and del Pozo, B. (2024). Drug decriminalization, fentanyl, and fatal overdoses in Oregon. *JAMA Network Open*, 7(9), e2431612. https://doi.org/10.1001/jamanetworkopen.2024.31612

Replication data: Harvard Dataverse, https://doi.org/10.7910/DVN/JT4CDS. The file is copied unchanged; it was built by the replication package's `Data_merging_file.R` from CDC WONDER mortality data and NFLIS drug seizure data.

One row per state (50 states + DC) and half-year, 2008 H1 to 2022 H2 (30 half-years). The columns used here:

| Column | Meaning |
|---|---|
| `State`, `ST` | state name and abbreviation |
| `Year`, `half` | calendar year and half (1 = January-June, 2 = July-December) |
| `trend`, `Year_half` | period index (1 = 2008 H1) and label |
| `od_death_rate_half` | drug overdose deaths per 100,000 residents in the half-year (CDC WONDER) |
| `nflis_fentanyl_percent_total` | percentage of drug items seized and analysed in the state that contained fentanyl or related substances (NFLIS) |
| `treatment_oregon` | 1 for Oregon from 2021 H1 (Measure 110 took effect on 1 February 2021) |

Following Zoorob et al., the analyses drop Washington, whose supreme court struck down the state's drug possession law in 2021, effectively decriminalizing possession for about four months. Two suppressed North Dakota death counts were replaced with random values below 10 by the original merging script.

Used by `docs/oregon_fentanyl/oregon_fentanyl_case_study.qmd` and `simulations/oregon_fentanyl_validation.R`, which read it from this path relative to the project root.

## License and sources

- **Data:** the replication data are released under **CC0 1.0** (public-domain dedication) on Harvard Dataverse: Zoorob, M., Park, J. N., Kral, A., Lambdin, B., and del Pozo, B. (2024). Replication Data for: Drug Decriminalization, Fentanyl, and Fatal Overdoses in Oregon. Version 1.0. https://doi.org/10.7910/DVN/JT4CDS. The file is redistributed here unchanged.
- **Article:** Zoorob et al. (2024), *JAMA Network Open*, is an open access article distributed under the terms of the CC-BY License (https://doi.org/10.1001/jamanetworkopen.2024.31612).
- **Underlying sources:** overdose deaths and population from the unrestricted (public) multiple cause of death files of CDC WONDER (https://wonder.cdc.gov); fentanyl share of seized drugs from the public data of the National Forensic Laboratory Information System (NFLIS, https://www.nflis.deadiversion.usdoj.gov). The file contains no suppressed CDC WONDER counts: the only two death counts below 10 (North Dakota) are random values the original merging script substituted for suppressed cells.

License status checked on 15 October 2026 (Dataverse metadata and the article's license statement); this is not legal advice.
