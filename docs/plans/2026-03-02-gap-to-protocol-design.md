# Gap-to-Protocol Pipeline — Design

**Goal:** Connect the DR Landscape gap detection to protocol auto-population and study search, so users can go from "this area needs a meta-analysis" to "here are the dose-ranging trials" in one click.

**Approach:** Wire existing components (gap scores, DR Landscape table, PICO builder, protocol generator, multi-source search) into a seamless pipeline. No new tabs or modals.

## Flow

1. DR Landscape "Top 10 Opportunities" table → each row gets a **"Plan Meta"** button
2. Click → auto-populates PICO fields from gap data (population from subcategory, intervention from drug class, comparator default, outcomes from cluster data)
3. Switches to Protocol tab with fields pre-filled
4. Protocol tab shows a **"Search Studies"** button
5. Search fires existing `searchCTGov()` + `searchPubMed()` with PICO terms
6. Results auto-filter to `doseRanging === true` trials only
7. Banner shows: "Showing X dose-ranging trials for [Drug Class] in [Population]"

## PICO Auto-Population Mapping

| Field | Source |
|---|---|
| Population | Subcategory label + keywords from `CARDIO_SUBCATEGORIES` |
| Intervention | Drug class name from landscape grouping |
| Comparator | Default: "placebo OR standard of care" |
| Outcome | From cluster `outcome_categories` if available, else blank |

## Components Modified

- `renderDRLandscape()` — add "Plan Meta" button to opportunity table rows
- New `planMetaFromGap(drugClass, subcategoryId, trials)` — extracts PICO, fills fields, switches tab
- Protocol tab — add "Search Studies" button that triggers search with current PICO
- Search results — add DR filter toggle + contextual banner

## No New Infrastructure

- Uses existing `savePICO()` for persistence
- Uses existing `searchCTGov()`, `searchPubMed()`, `buildPubMedQuery()`
- Uses existing `generateProtocol()` for export
- Uses existing gap scoring from `computeAllGapScores()`
