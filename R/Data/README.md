# Input data

`vaping_mccauley2023.csv` contains aggregate counts and percentages from
McCauley, Baiocchi, Cruse and Halpern-Felsher (2023), *Preventive Medicine Reports*,
33:102184 (doi:10.1016/j.pmedr.2023.102184). It contains no individual records.
The `source` column records provenance; `status` identifies published and
reconstructed counts.

`sample` distinguishes linked students (`linked`) from everyone observed
at each stage (`full`). `item` identifies the question; `n1`, `x1`, `n2`
and `x2` give the margins. `n00`, `n01`, `n10` and `n11` give paired
transition counts when available. Missing cells remain unknown.

For daily use and addiction and the definition of addiction, the full-sample
counts are inferred from rounded percentages and p-values. The selected
reconstruction minimises item nonresponse; `R/Data_preparation/reconstruct_counts.R`
enumerates all compatible alternatives for Table S2.
