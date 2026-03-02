# Gap-to-Protocol Pipeline Implementation Plan

> **For Claude:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to implement this plan task-by-task.

**Goal:** Wire DR Landscape gap detection to PICO auto-population, protocol tab, and study search.

**Architecture:** Add "Plan Meta" button to DR Landscape opportunity table. New `planMetaFromGap()` function populates PICO fields from gap data, switches to Protocol tab. Add "Search DR Studies" button to Protocol tab that fires existing search with DR filter.

**Tech Stack:** Vanilla JS, existing HTML structure, existing search APIs

---

### Task 1: Add "Plan Meta" button to DR Landscape opportunity table

**Files:**
- Modify: `metasprint-dose-response.html:20112` (the button cell in renderDRLandscape)

**Step 1:** In the `renderDRLandscape()` function, find the table row that builds the "Start DR Review" button (line ~20112). Add a second button "Plan Meta" next to it that calls `planMetaFromGap(escapedDrug, subcatId, dTrials)`.

### Task 2: Implement `planMetaFromGap()` function

**Files:**
- Modify: `metasprint-dose-response.html` (add function near savePICO ~line 20270)

**Step 1:** Create `planMetaFromGap(drugClass, subcatId, trials)` that:
- Looks up population from `CARDIO_SUBCATEGORIES` using subcatId
- Sets `#picoP` / `#protP` = subcategory label + keywords
- Sets `#picoI` / `#protI` = drugClass
- Sets `#picoC` / `#protC` = "placebo OR standard of care"
- Sets `#picoO` / `#protO` = outcomes from trial cluster data (if available)
- Calls `savePICO()` to persist
- Switches to the Protocol tab by clicking `#tab-protocol`
- Shows a banner in the protocol section: "Planning DR meta-analysis for [drugClass] in [population]"

### Task 3: Add "Search DR Studies" button to Protocol tab

**Files:**
- Modify: `metasprint-dose-response.html:~1218` (phase-protocol section)

**Step 1:** After the PICO fields in the protocol section, add a "Search DR Studies" button that:
- Calls a new `searchDRStudies()` function
- This function calls `searchAll()` then applies a post-filter for `doseRanging === true`
- Shows a contextual banner with result count

### Task 4: Add DR filter + banner to search results

**Files:**
- Modify: `metasprint-dose-response.html` (add function near searchAll)

**Step 1:** Create `searchDRStudies()` that:
- Validates PICO P or I is non-empty
- Calls the existing search pipeline
- After results load, filters to `doseRanging === true` trials
- Displays banner: "Showing X dose-ranging trials for [Intervention] in [Population]"

### Task 5: Tests and verification

- Add inline self-test `_testGapToProtocol()` verifying PICO auto-population
- Add Selenium test in TestExceedDosresmeta class
- Verify div balance, no JS errors
