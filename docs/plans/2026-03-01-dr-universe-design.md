# Dose-Response Universe Enhancement Design

**Date**: 2026-03-01
**App**: `metasprint-dose-response.html` (~20,905 lines)
**Scope**: Add dose-response awareness to 4 pipeline stages (Universe, Visualization, Search, Screening)

## Context

The Discover/Universe phase has 9 hardcoded cardio subcategories, a 3-layer trial classifier, and 8 visualization tabs. None of this is dose-response aware. Users cannot identify dose-ranging trials or generate DR-specific search strategies.

## Enhancement 1: Dose-Ranging Trial Detection

### What
Parse CT.gov arm descriptions to identify trials with multiple dose levels.

### Implementation
- After trial classification, run `detectDoseRanging(trial)` on each trial
- Regex-extract dose values from arm labels: `(\d+(?:\.\d+)?)\s*(mg|mcg|ug|g|mg\/kg|mg\/m2|iu|units?)`
- Count distinct dose values (excluding placebo arms)
- Flag as `doseRanging: true` if 3+ distinct dose levels found
- Store: `trial.doseRanging = bool`, `trial.doseLevels = [10, 20, 40]`, `trial.doseUnit = 'mg'`

### UI
- "DR" badge on trial cards in Universe grid
- Filter toggle: "Show only dose-ranging trials"
- Drug class cards show "X dose-ranging trials" count

## Enhancement 2: Dose-Response Landscape Visualization Tab

### What
A new 9th visualization tab in Discover: "DR Landscape"

### Content
- **Dose-Ranging Heatmap**: Drug class (rows) vs dose range buckets (columns), color = trial count
- **Dose Distribution Chart**: For a selected drug class, show dose levels across all dose-ranging trials (dot plot or histogram)
- **Top DR Opportunities**: List of drug classes with most dose-ranging trials and least existing dose-response meta-analyses (gap detection)
- **Summary stats**: Total DR trials found, % of universe that is dose-ranging, top 5 drug classes

### Interaction
- Click a drug class in the heatmap → drill down to dose distribution
- Click a specific trial → shows arm details with extracted dose levels
- "Start DR Review" button → populates PICO with dose-response framing

## Enhancement 3: DR-Specific Search Strategy

### What
When starting a review from the Universe (especially from a dose-ranging trial), auto-generate dose-response search terms.

### Implementation
- In `startReviewFromUniverse()`, detect if the source is a dose-ranging drug class
- Auto-add to PubMed search query:
  - MeSH: `"Dose-Response Relationship, Drug"[MeSH]`
  - Keywords: `dose-response OR dose-ranging OR dose-finding OR dose-titration OR dose-dependent`
- Auto-populate PICO:
  - I: "[Drug] at multiple doses ([dose levels])"
  - C: "placebo or lowest dose ([min dose])"
  - O: preserve existing outcome template
- Add DR search template as a search preset (dropdown option)

## Enhancement 4: DR Screening Relevance Boost

### What
Boost screening scores for abstracts containing dose-response signals.

### Implementation
- DR signal terms (regex): `dose.response`, `dose.rang`, `dose.find`, `dose.depend`, `dose.level`, `dose.titrat`, `mg/day`, `mg/kg`, `multiple.dos`, `graded.dos`
- In the PICO relevance scorer, add a DR bonus:
  - Match 1 DR term: +1 score
  - Match 2+ DR terms: +2 score
  - Match "dose-response" exactly: +3 score
- Add "DR" badge on references that score DR bonus
- DR score visible in screening table as a column or tooltip

## Safety / Constraints
- Dose extraction regex must handle: "10 mg", "10mg", "0.5 mg/kg", "100 mcg", "1000 IU"
- Must not false-positive on non-dose numbers (enrollment, duration, age)
- Arm label parsing: skip arms containing "placebo", "sham", "control", "standard of care"
- All changes are additive — do not break existing Universe functionality
