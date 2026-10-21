# 30-seed robustness analysis: interpretation

## Design

Thirty new A/B synthetic slice-pair realizations were prespecified using A seeds 1000–1029 and B seeds 2000–2029. FGW (γ=0) and LROT (γ=0.1) were evaluated on the same A/B pair in each replicate. The historical 42/123 run and the three earlier exploratory pairs were excluded from this primary set. Solver and data-generation settings were held fixed. The archived feature order and historical synthetic group membership were frozen to keep the comparison aligned with the reproducible 0.631 result.

## Results

Mean accuracy was 0.265 (sample SD 0.195) for FGW and 0.405 (sample SD 0.232) for LROT. The mean paired difference was +0.140 (sample SD 0.225; 95% paired t interval +0.056 to +0.224; percentile-bootstrap sensitivity interval +0.064 to +0.222). LROT was higher in 24/30 replicates and lower in 6/30. The median paired difference was +0.056, with a range from −0.191 to +0.533.

## What this supports

Under this exact frozen synthetic setup, the average paired accuracy favors LROT. The positive average does not mean the effect is uniform: one fifth of the prespecified realizations favor FGW, and the paired differences vary widely. The historical 0.631 versus 0.363 result is reproducible for its fixed realization, but its approximately 74% relative increase should not be presented as the typical or general effect. The 30-run mean difference is the relevant summary for this seed-robustness analysis.

The group membership is frozen from the historical simulator for forensic reproducibility; it has not been shown to represent biologically validated pathways. Consequently, this analysis tests technical sensitivity to simulation seeds under that historical design. It does not establish pathway-level validity or support a broad biological superiority claim.

## Manuscript reporting

The supplementary table and Figure S23 report the thirty prespecified paired runs. The manuscript should distinguish their mean paired effect and heterogeneity from the reproducible historical 42/123 realization. This technical simulation check does not validate biological pathway groupings.

