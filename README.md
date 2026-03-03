# MetaSprint Dose-Response

A single-file browser application for **dose-response meta-analysis**, guiding researchers through a structured 40-day sprint from evidence discovery to manuscript writing.

## Quick Start

Open `metasprint-dose-response.html` in any modern browser (Chrome, Firefox, Edge). No server, build step, or installation required.

## Features

### Dose-Response Models
- **8 models**: Linear, Quadratic, Emax, Log-Linear, Exponential, Hill, Fractional Polynomial (FP1/FP2), Restricted Cubic Splines (RCS-3/RCS-4)
- **Estimation**: ML, REML, and fixed-effects via profile likelihood with golden section optimization
- **One-stage mixed effects**: Random intercept + slope by study
- **Model averaging**: AIC-weighted predictions across all fitted models (Burnham & Anderson 2002)
- **Bayesian DR**: Laplace approximation posterior with credible interval bands

### Statistical Engine
- Greenland-Longnecker covariance reconstruction for correlated dose contrasts
- Block-diagonal GLS with full within-study covariance
- Prediction intervals (Higgins-Thompson-Spiegelhalter)
- Dose-finding metrics: MED, ED50, ED90, therapeutic window, NOAEL
- Publication bias: dose-arm funnel plot + Egger-style regression (t-distribution)
- Bootstrap breakpoint CI (500 cluster resamples)
- Leave-one-out sensitivity analysis

### Sprint Workflow
- **Discover**: RCT landscape explorer (ClinicalTrials.gov, PubMed, OpenAlex)
- **Protocol**: PICO definition with PROSPERO integration
- **Search**: Multi-source search with deduplication
- **Screen**: Keyboard-driven title/abstract screening (I/E/M shortcuts)
- **Extract**: Dose-level data entry with 2x2 auto-computation
- **Analyze**: Model fitting, comparison, and visualization
- **Write**: Auto-generated Methods and Results sections
- **Insights**: 10 evidence analysis tools (trust radar, harm-benefit balance, multiverse, NNT calculator, equity assessment)

### Evidence Landscape
- **Ayat Universe**: Interactive force-directed network of cardiology RCT clusters
- **Al-Burhan**: Living meta-analysis engine with automatic pooling
- **Gap-to-Protocol**: Identifies evidence gaps and auto-populates PICO

## Running Tests

### Requirements
- Python 3.10+
- Chrome browser
- `pip install selenium pytest`

### Selenium Integration Tests (81 tests)
```bash
cd metasprint-dose-response
python -m pytest test_dose_response_models.py -v
```

### Inline Self-Tests (446 tests)
Open the browser console and run:
```javascript
_runAllDRTests()  // Returns true if all 14 suites pass
```

### R Cross-Validation
Requires R with `dosresmeta` and `jsonlite`:
```bash
Rscript validate_vs_dosresmeta.R
```
This generates `validation_reference.json` with reference values from dosresmeta2 v2.2.0 for the alcohol_cvd and coffee_mort datasets.

## Validated Against R

The JS implementations are cross-validated against R `dosresmeta2` v2.2.0 with tolerance 1e-4:
- Linear, Quadratic, RCS-3 coefficients and standard errors
- ML and REML tau-squared estimation
- AIC/BIC model comparison
- Predictions at specified dose points
- Leave-one-out analysis (6 studies)
- Greenland-Longnecker covariance matrices

## Architecture

Single HTML file (~25,800 lines) with zero external dependencies:
- Vanilla JS (no frameworks)
- CSS custom properties for theming (light/dark mode)
- localStorage for persistence (with in-memory fallback)
- IndexedDB for trial universe cache
- Fetch API for ClinicalTrials.gov and PubMed

## License

Research and educational use. Not for clinical decision-making.
