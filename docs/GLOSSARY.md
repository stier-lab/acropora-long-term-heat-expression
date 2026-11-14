# Glossary

*Created: 2026-11-14 · Last updated: 2026-11-14*

Plain-language meanings for the terms used in this repo's code, tables and figures. Column names in the data are in `code font`.

## The experiment

| Term | Meaning |
|---|---|
| Fragment | One small branch of coral in its own spot in a tank. The unit everything is measured on. |
| Ambient / heated | The two temperature treatments: 28 °C and 31 °C. Heated tanks were ramped up 1 °C per day. |
| Tank | Temperature was set per tank (four tanks per temperature), so tanks, not fragments, are the true replicates for the temperature effect. |
| Wounded / clipped | About 1 cm cut from the branch tip on Day 0. Unwounded fragments were not cut. |
| Day 0, Day −1 | Day 0 is the day of clipping. Day −1 is just before clipping, not an unheated baseline. |
| Source patch | The field patch (A, C or D) a fragment's parent branch came from. Older column names say `thicket` or `genet`. These are field labels, not confirmed genotypes. |
| Host cluster | A group of fragments that DNA markers (SNPs) show are genetically similar. These come from Rachael Bay's lab and do not line up exactly with source patches. |
| Main experiment / microscope experiment | 192 fragments in tanks (48 followed daily), versus a separate set of 16 wounded fragments photographed under a microscope. |

## Coral anatomy and regeneration stages

| Term | Meaning |
|---|---|
| Polyp | One coral animal; a branch is covered in them. |
| Axial polyp | The polyp at the very tip of a branch. A new one forming on the cut surface is regeneration stage 1. |
| Corallite | The skeletal cup a polyp sits in. Axial = the cup at the tip; radial = the cups along the sides. |
| Budding / radial corallites form | New side cups appearing near the regrown tip (stage 5). This is the step heat delayed most. |
| Coenosarc | The thin sheet of living tissue between polyps. |
| Pigment over wound | Coloured tissue growing back over the cut surface. |
| AEF | *Acropora*-eating flatworm. Its signs were scored daily. |

## Measurements

| Term | Meaning |
|---|---|
| Fv/Fm | How efficiently the algae inside the coral turn light into energy, from 0 to 1. Healthy corals here read about 0.65–0.68; lower means heat stress. |
| PAM | The fluorometer (light meter) that measures Fv/Fm. |
| Symbiont / zooxanthellae | The single-celled algae living inside coral tissue that give it colour and food. Counted as million cells per cm² of coral surface. |
| Symbiont read group (A–D) | Which group of algae the RNA-seq reads matched. Only group A is confirmed (*Symbiodinium*); the others still need checking, and the shares are rough. |
| Color score / D-scale | A 1–6 score read off a colour card. Lower means paler (bleaching). |
| Buoyant weight | Weighing the coral underwater. With the water's density this gives skeleton dry mass without removing the coral. |
| Aragonite plug | The small mineral base each fragment was glued to. Its mass is subtracted before calculating growth. |
| Wax dipping | Dipping a fragment in wax and weighing the wax to estimate its surface area. |
| Apex | The tank controller that logged each tank's temperature every 10 minutes. |
| YSI | A handheld meter used for spot checks of temperature, oxygen, pH and salinity. |
| RNA-seq | Sequencing a sample's RNA to see which genes are switched on. Processed at UC Davis. |

## Statistics

| Term | Meaning |
|---|---|
| Mixed model | A regression that allows each tank (or fragment) its own baseline and trend, so repeated measures from the same tank aren't counted as independent. |
| Temperature × day | Whether the change over time differs between temperatures. This is the main test for Fv/Fm, color and symbionts. |
| Permutation test | Reshuffles which tanks were heated, many times, to see how often a difference as large as the real one appears by chance. With 4 vs 4 tanks there are only 70 possible shuffles. |
| Likelihood-ratio test | Asks whether adding a term (for example, source patch × temperature) makes a model fit clearly better. |
| Time-to-event / survival analysis | Methods for "how long until a fragment reaches a stage", which can handle fragments that never reached it by Day 15. |
| Censored | A fragment that hadn't reached a stage by the last check: we only know it would take longer than 15 days. *Interval-censored* means we know it happened between two checks but not exactly when. |
| Kaplan–Meier curve | A step plot of the share of fragments that had reached a stage by each day. |
| Hazard ratio | From a Cox model: how much more (or less) likely a heated fragment is to reach a stage on any given day than an ambient one. Below 1 means slower. |
| Time ratio | From the log-logistic model: how much longer heated fragments take to reach a stage. Above 1 means slower. |
| DHARMa checks | Simulation-based checks that a model's leftover scatter looks the way it should. Results can shift slightly between package versions, which is why versions are pinned. |
| Fisher's exact test | Compares two counts out of small totals, such as 12/12 vs 4/12 fragments. |
| Holm adjustment | A correction to p-values when several related tests are run at once. |
| renv | The R tool that records the exact package versions (`renv.lock`) so the analysis reruns the same way. |
