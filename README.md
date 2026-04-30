# Thermal Load Sensitivity (TLS) of Photosystem II via Chlorophyll Fluorescence

[![DOI](https://img.shields.io/badge/DOI-10.64898%2F2026.04.09.717599-blue)](https://doi.org/10.64898/2026.04.09.717599) [![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

## Overview

This repository contains code, data, and protocols for quantifying thermal tolerance of photosystem II (PSII) using the **Thermal Load Sensitivity (TLS)** approach, which integrates both temperature intensity and exposure duration to derive time-based critical temperature thresholds.

### Associated Publication

**Title:** Towards a standard approach to investigating the Thermal Load Sensitivity of photosystem II via chlorophyll fluorescence

**Authors:** Pieter A. Arnold, Rosalie J. Harris, Sabina M. Aitken, Max M. Hoek, Alicia M. Cook, Andy Leigh, Adrienne B. Nicotra

**Journal:** bioRxiv preprint

**Abstract:** Evaluating the drivers of variation in plant thermal tolerance limits requires a clearer understanding of how methodological matters can lead to different tolerance estimates. Chlorophyll fluorometry – to measure the temperature-dependent change in *F*<sub>V</sub>/*F*<sub>M</sub> – is a well-established approach to derive tolerance thresholds of photosystem II (PSII) in plants, but one-off, time-specific thermal exposures do not consider the fundamental dose-dependent effect of heat. The resurgent thermal death time (TDT) approach integrates both the temperature intensity and the exposure duration to derive time-based critical temperature thresholds and sensitivity parameters. We build upon this foundation to develop a protocol for evaluating thermal load sensitivity (TLS; non-lethal heat stress) of PSII in plants. Through five experiments across four diverse species, we tested the moderating effects of light, leaf sectioning, time since collection, and the temporal dynamics of *F*<sub>V</sub>/*F*<sub>M</sub> recovery. There were dramatic changes in tolerance threshold estimates based on thermal load (i.e. dose-dependent) effects on *F*<sub>V</sub>/*F*<sub>M</sub>, and strong effects of light intensity during heat and the presence of light post-heat. We offer recommendations pertaining to method implementation and discuss future empirical avenues. Appraising cumulative heat stress will enhance the utility of thermal tolerance estimates – the TLS approach outlined here moves us toward a new standard.

---

## Core concepts and glossary

**PPFD** (µmol m<sup>-2</sup> s<sup>-1</sup>) = photosynthetic photon flux density, i.e. light intensity.

***F*<sub>0</sub>** (arbitrary units) = basal or minimal fluorescence reemitted by a dark-adapted leaf sample (with reaction centres open) measured under a weak blue measuring light.

***F*<sub>M</sub>** (arbitrary units) = maximum fluorescence reemitted by a leaf sample under a saturating pulse of light.

***F*<sub>V</sub>/*F*<sub>M</sub>** (ratio) = maximum quantum yield of PSII, calculated as (*F*<sub>M</sub> - *F*<sub>0</sub>) / *F*<sub>M</sub>.

**Proportional *F*<sub>V</sub>/*F*<sub>M</sub>** = calculated as post-heat *F*<sub>V</sub>/*F*<sub>M</sub> / initial *F*<sub>V</sub>/*F*<sub>M</sub>, with values normalised between 0 and 1.

**Tolerance thresholds: *T*<sub>10</sub>, *T*<sub>50</sub>, *T*<sub>0.3</sub>, *T*<sub>90</sub>** (°C) = predicted temperature at which proportional *F*<sub>V</sub>/*F*<sub>M</sub> declines by the subscript percentage amount (e.g. *T*<sub>10</sub> is the temperature at which *F*<sub>V</sub>/*F*<sub>M</sub> declines by 10%). *T*<sub>0.3</sub> is calculated on the absolute *F*<sub>V</sub>/*F*<sub>M</sub> values, such that it is the temperature at *F*<sub>V</sub>/*F*<sub>M</sub> = 0.3. Typically derived from a dose-response curve (logistic, or quasi-binomial, generalised linear model).

**TLS parameter *CT*<sub>max_1m</sub>** (°C) = theoretical heat tolerance limit from 1-min exposure. This is the *y*-intercept derived from fitting a linear regression model to the tolerance thresholds calculated at each exposure duration when in minutes (log<sub>10</sub> scale).

**TLS parameter *CT*<sub>max_1h</sub>** (°C) = biologically relevant heat tolerance limit from 1 h exposure. Derived from fitting a linear regression model to the tolerance thresholds calculated at each exposure duration (log<sub>10</sub> scale).

**TLS parameter *z*** (dimensionless) = thermal sensitivity parameter that describes the change in tolerance threshold for 10-fold increase in exposure duration. This is the slope derived from fitting a linear regression model to the tolerance thresholds calculated at each exposure duration (log<sub>10</sub> scale).

### Threshold Calculations

For each sample and time point, threshold temperatures (e.g., *T*<sub>50</sub>) are calculated by fitting dose-response curves:
*F*<sub>V</sub>/*F*<sub>M</sub> = *d* / (1 + exp(*b* × (log(Temp) - log(ED<sub>50</sub>))))
where *d* is the upper asymptote, *b* is the slope, and ED<sub>50</sub> is the temperature causing 50% response. The temperature causing any threshold response level is then extracted by inverting this equation numerically using the ED() function from the drc package in R.

These threshold temperatures are then regressed against log<sub>10</sub>(time) to extract *CT*<sub>max_1m</sub> and z.

### Thermal Death Time (TDT) Model

The TDT approach fits a linear relationship between temperature and log-transformed time:
log<sub>10</sub>(time) = (Temperature - *CT*<sub>max_1m</sub>) / z

### Thermal Load Sensitivity (TLS) Model

Akin to the TDT approach at its core, it fits a linear relationship between temperature and log-transformed time. Focus on sublethal effects, with capacity to disentangle damage and repair dynamics. 

---

## Output Files

### Threshold Data
- One row per sample × time × threshold
- Raw threshold temperatures used for TDT fitting

### Diagnostic Plots
- Dose-response curves
- TDT regression plots
- Residual diagnostics
- Quality control summaries

### TLS Metrics
- One row per sample × threshold × response type
- Columns: `CTmax_1m`, `CTmax_1h`, `z`, `r_squared`, `valid_fit`

---

## Example

- Rmarkdown and html render demonstrates the workflow through applying the functions in a pipeline to example data.
---

## Related Resources
- **Chlorophyll fluorescence methods:** [Maxwell & Johnson 2000](https://doi.org/10.1093/jxb/51.345.659)
- **TDT framework:** [Rezende et al. 2014](https://doi.org/10.1111/1365-2435.12268)
- **TDT framework:** [Jørgensen et al. 2021](https://doi.org/10.1038/s41598-021-92004-6)
- **TLS framework:** [Arnold et al. 2025](https://doi.org/10.1111/gcb.70315)

