# YSI update, October 27, 2026

Source: [YSI daily log](https://docs.google.com/spreadsheets/d/11MrOICcoIp-lOJ9HiHcB9xWzf-_MgeMLQbt79Bmbb84/edit), supplied in Molly's October 27 email. Formatted-value CSV snapshots of Data and READ.ME are retained under `data/raw/ysi/` as `daily_log_2026-10-27.csv` and `daily_log_units_2026-10-27.csv`.

SHA256: Data `0e5f7968d98c52fdbcf4ec43e749b985f2c674648f9d7c291f034efd18df5557`; units `127fc7159b1f008e300a72c0cc81a26735cdaccbf8f040f465484eb457c92028`.

The sheet specifies temperature in Fahrenheit (converted to Celsius), dissolved oxygen concentration in mg/L, and salinity as unitless. Dissolved oxygen saturation is separately recorded as percent. The pH scale and meter calibration procedure are not documented. Do not substitute mg/mL or invent a salinity unit or pH scale.

The inclusive June 4–19, 2025 window contains 112 complete measurements: eight tanks on 14 dates. June 4 and June 7 are absent. Measurements occur at varying times, not consistently at 18:00. No missing dates are imputed. The invalid pH string `8..02` for tank 6 on May 29 is preserved in raw data, logged in the parse audit, and left missing in processed data; it falls outside the study summary.

Script 09 replaces the old export as the analysis input without appending duplicate observations. Finite overlapping values agree with the legacy sheet. It writes study and treatment means, SD across tank-date readings, individual tank means, and SD/SE among those tank means. The latter SE uses eight tanks overall or four per temperature; repeated readings are not independent experimental replicates. Overall SE combines both temperature treatments and should not be interpreted as treatment-effect uncertainty. These are descriptive summaries, not tests of environmental equivalence.

Dissolved oxygen concentration is lower in heated tanks. Keep treatment-specific summaries visible rather than claiming that all conditions except temperature were identical. Continuous temperature coverage comes from the separate Apex log.
