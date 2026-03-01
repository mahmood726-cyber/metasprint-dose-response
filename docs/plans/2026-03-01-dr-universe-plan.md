# DR-Aware Universe Enhancement Implementation Plan

> **For Claude:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to implement this plan task-by-task.

**Goal:** Add dose-response awareness to the Discover/Universe phase — detect dose-ranging trials from CT.gov, add a DR Landscape visualization tab, generate DR-specific search strategies, and boost DR-relevant papers in screening.

**Architecture:** Layer 4 enhancements onto the existing cardiovascular Universe infrastructure. The `detectDoseRanging(trial)` function parses CT.gov arm labels to extract dose levels. The DR Landscape tab visualizes dose-ranging trial distribution. Search and screening get DR-specific term injection.

**Tech Stack:** Vanilla JS (single-file HTML), SVG for DR Landscape visualization, CT.gov API arm data.

**Key file:** `C:\Users\user\Downloads\metasprint-dose-response\metasprint-dose-response.html` (~20,905 lines)

**Key line references:**
- Universe view tab buttons: ~line 1127-1134
- `classifyTrial()`: ~line 2455
- `computeAllGapScores()`: ~line 4101
- Screening scorer (picoScore/rctSignal/cardioSignal): ~line 4932, 5294
- `runAutoScreen()`: ~line 5494
- `renderUniverseGrid()`: ~line 13065
- `startReviewFromUniverse()`: ~line 13545
- `switchUniverseView()`: ~line 13589

**CRITICAL — Single-file HTML rules:**
- Never write literal `</script>` inside `<script>` blocks
- After structural HTML edits, verify div balance
- All element IDs globally unique
- All function names globally unique
- Use `?? fallback` not `|| fallback` for numeric values

---

### Task 1: Dose-Ranging Trial Detection Function

**Files:**
- Modify: `metasprint-dose-response.html` — insert after `classifyTrial()` (~line 2516)

**Step 1: Add `detectDoseRanging(trial)` function**

Insert after `classifyTrial()` ends:

```javascript
// ── Dose-Ranging Trial Detection ──
const DOSE_UNIT_REGEX = /(\d+(?:\.\d+)?)\s*(mg|mcg|ug|µg|g|mg\/kg|mg\/m2|iu|units?|ml|mmol|nmol|pmol)(?:\b|\/)/gi;
const PLACEBO_ARM_REGEX = /\bplacebo\b|\bsham\b|\bcontrol\b|\bstandard\s+of\s+care\b|\bSOC\b|\bno\s+treatment\b/i;

function detectDoseRanging(trial) {
  // Parse CT.gov arm descriptions to find distinct dose levels
  const arms = trial.arms || [];
  const interventions = trial.interventions || [];
  const allArmText = arms.map(a => {
    const label = typeof a === 'string' ? a : (a.label || a.name || '');
    return label;
  });

  // Also check intervention names
  const allIvText = interventions.map(iv => {
    const name = typeof iv === 'string' ? iv : (iv.name || '');
    return name;
  });

  const allText = [...allArmText, ...allIvText];
  const doseLevels = new Set();
  let doseUnit = '';
  let armCount = 0;

  for (const text of allText) {
    // Skip placebo/control arms
    if (PLACEBO_ARM_REGEX.test(text)) continue;
    armCount++;

    // Extract dose values
    let match;
    const regex = new RegExp(DOSE_UNIT_REGEX.source, DOSE_UNIT_REGEX.flags);
    while ((match = regex.exec(text)) !== null) {
      const doseVal = parseFloat(match[1]);
      const unit = match[2].toLowerCase();
      if (isFinite(doseVal) && doseVal > 0) {
        doseLevels.add(doseVal);
        if (!doseUnit) doseUnit = unit;
      }
    }
  }

  const levels = Array.from(doseLevels).sort((a, b) => a - b);
  const isDoseRanging = levels.length >= 3;

  trial.doseRanging = isDoseRanging;
  trial.doseLevels = levels;
  trial.doseUnit = doseUnit;
  trial.doseArmCount = armCount;

  return isDoseRanging;
}
```

**Step 2: Call detection after classification**

Find where `classifyTrial(trial)` is called (search for `classifyTrial(` — it should be in the trial loading/processing pipeline). After each `classifyTrial(trial)` call, add:

```javascript
detectDoseRanging(trial);
```

**Step 3: Add self-test**

```javascript
function _testDoseRanging() {
  const results = [];

  // Test 1: Multi-dose trial
  const t1 = {
    arms: [
      { label: 'Drug X 10 mg' },
      { label: 'Drug X 20 mg' },
      { label: 'Drug X 40 mg' },
      { label: 'Placebo' }
    ],
    interventions: []
  };
  detectDoseRanging(t1);
  results.push({ test: 'multi-dose detected', pass: t1.doseRanging === true });
  results.push({ test: '3 levels', pass: t1.doseLevels.length === 3 });
  results.push({ test: 'unit=mg', pass: t1.doseUnit === 'mg' });
  results.push({ test: 'levels sorted', pass: t1.doseLevels[0] === 10 && t1.doseLevels[2] === 40 });

  // Test 2: Not dose-ranging (only 2 arms)
  const t2 = {
    arms: [{ label: 'Drug Y 50 mg' }, { label: 'Placebo' }],
    interventions: []
  };
  detectDoseRanging(t2);
  results.push({ test: 'single-dose not DR', pass: t2.doseRanging === false });

  // Test 3: Dose in intervention names
  const t3 = {
    arms: [],
    interventions: [
      { name: 'Empagliflozin 10 mg' },
      { name: 'Empagliflozin 25 mg' },
      { name: 'Empagliflozin 50 mg' }
    ]
  };
  detectDoseRanging(t3);
  results.push({ test: 'intervention names detected', pass: t3.doseRanging === true });

  // Test 4: mg/kg units
  const t4 = {
    arms: [
      { label: '0.5 mg/kg' }, { label: '1.0 mg/kg' }, { label: '2.0 mg/kg' }
    ],
    interventions: []
  };
  detectDoseRanging(t4);
  results.push({ test: 'mg/kg detected', pass: t4.doseRanging === true });
  results.push({ test: 'mg/kg unit', pass: t4.doseUnit === 'mg/kg' });

  // Test 5: No dose info at all
  const t5 = { arms: [{ label: 'Active Treatment' }, { label: 'Control' }], interventions: [] };
  detectDoseRanging(t5);
  results.push({ test: 'no dose = not DR', pass: t5.doseRanging === false });
  results.push({ test: 'empty levels', pass: t5.doseLevels.length === 0 });

  console.table(results);
  return results.every(r => r.pass);
}
```

**Step 4: Commit**

```bash
git commit -m "feat: add dose-ranging trial detection from CT.gov arm labels"
```

---

### Task 2: Universe Grid DR Badges + Filter

**Files:**
- Modify: `metasprint-dose-response.html` — update `renderUniverseGrid()` (~line 13065)

**Step 1: Add DR count to drug class cards**

In `renderUniverseGrid()`, after the cluster building loop, compute the DR count per cluster. Find where each card's HTML is built (search for the card template inside `renderUniverseGrid`). Add a "DR" badge showing the count of dose-ranging trials:

```javascript
// Inside the card template, after existing badges:
const drCount = clusterTrials.filter(t => t.doseRanging).length;
const drBadge = drCount > 0
  ? '<span style="display:inline-block;background:#7c3aed;color:#fff;padding:1px 6px;border-radius:8px;font-size:0.7rem;margin-left:4px" title="' + drCount + ' dose-ranging trials">' + drCount + ' DR</span>'
  : '';
```

Insert `drBadge` into the card HTML where badges are rendered.

**Step 2: Add DR filter toggle**

Find the filter/sort controls area in the Universe grid (search for sort buttons or filter area near the grid). Add a toggle:

```html
<label style="font-size:0.8rem;cursor:pointer;margin-left:12px" title="Show only drug classes with dose-ranging trials">
  <input type="checkbox" id="drFilterToggle" onchange="filterUniverseGridDR(this.checked)"> DR only
</label>
```

**Step 3: Add filter function**

```javascript
function filterUniverseGridDR(onlyDR) {
  if (typeof currentGridTrials === 'undefined' || !currentGridTrials) return;
  if (onlyDR) {
    const drTrials = currentGridTrials.filter(t => t.doseRanging);
    renderUniverseGrid(drTrials, currentGridSubcat);
  } else {
    renderUniverseGrid(currentGridTrials, currentGridSubcat);
  }
}
```

**Step 4: Commit**

```bash
git commit -m "feat: add DR badges and filter to Universe grid cards"
```

---

### Task 3: DR Landscape Visualization Tab

**Files:**
- Modify: `metasprint-dose-response.html`:
  - HTML: add tab button after "Phase Pipeline" (~line 1134)
  - HTML: add container div for the view
  - JS: add `renderDRLandscape()` function
  - JS: update `switchUniverseView()` (~line 13589)

**Step 1: Add tab button**

After the "Phase Pipeline" button at ~line 1134, add:

```html
<button class="view-tab" data-view="drlandscape" onclick="switchUniverseView('drlandscape')" role="tab" aria-selected="false" title="DR Landscape: Dose-ranging trial distribution by drug class — find dose-response opportunities">DR Landscape</button>
```

**Step 2: Add container div**

Find the existing view containers (search for `id="universeViewAyat"` or similar). Add alongside them:

```html
<div id="universeViewDRLandscape" class="universe-view" style="display:none"></div>
```

**Step 3: Update `switchUniverseView()` to handle 'drlandscape'**

In `switchUniverseView()`, add a case for `viewName === 'drlandscape'` that shows the container and calls `renderDRLandscape()`.

**Step 4: Add `renderDRLandscape()` function**

```javascript
function renderDRLandscape(trials) {
  const el = document.getElementById('universeViewDRLandscape');
  if (!el) return;

  const drTrials = (trials || currentGridTrials || []).filter(t => t.doseRanging);
  const totalTrials = (trials || currentGridTrials || []).length;

  if (drTrials.length === 0) {
    el.innerHTML = '<p style="color:var(--text-muted);padding:24px">No dose-ranging trials found in the current universe. Try loading a different subspecialty.</p>';
    return;
  }

  // Group DR trials by drug class (intervention)
  const byDrug = {};
  for (const t of drTrials) {
    const ivNames = (t.interventions || [])
      .map(iv => typeof iv === 'string' ? iv : (iv.name || ''))
      .filter(n => !PLACEBO_ARM_REGEX.test(n));
    const drugLabel = ivNames[0] || 'Unknown';
    if (!byDrug[drugLabel]) byDrug[drugLabel] = [];
    byDrug[drugLabel].push(t);
  }

  // Sort by trial count descending
  const drugEntries = Object.entries(byDrug).sort((a, b) => b[1].length - a[1].length);
  const topN = drugEntries.slice(0, 20); // Top 20 drug classes

  // Summary stats
  const pctDR = totalTrials > 0 ? (drTrials.length / totalTrials * 100).toFixed(1) : '0';
  let html = '<div style="padding:16px">';
  html += '<h3 style="font-size:1rem;margin-bottom:12px">Dose-Response Landscape</h3>';
  html += '<div style="display:flex;gap:24px;margin-bottom:16px;flex-wrap:wrap">';
  html += '<div style="padding:12px 16px;background:var(--surface);border:1px solid var(--border);border-radius:var(--radius)"><div style="font-size:1.5rem;font-weight:700;color:#7c3aed">' + drTrials.length + '</div><div style="font-size:0.8rem;color:var(--text-muted)">Dose-ranging trials</div></div>';
  html += '<div style="padding:12px 16px;background:var(--surface);border:1px solid var(--border);border-radius:var(--radius)"><div style="font-size:1.5rem;font-weight:700;color:var(--primary)">' + pctDR + '%</div><div style="font-size:0.8rem;color:var(--text-muted)">of universe</div></div>';
  html += '<div style="padding:12px 16px;background:var(--surface);border:1px solid var(--border);border-radius:var(--radius)"><div style="font-size:1.5rem;font-weight:700;color:var(--success)">' + topN.length + '</div><div style="font-size:0.8rem;color:var(--text-muted)">Drug classes with DR data</div></div>';
  html += '</div>';

  // Bar chart: drug classes by DR trial count
  const maxCount = topN.length > 0 ? topN[0][1].length : 1;
  html += '<h4 style="font-size:0.9rem;margin:16px 0 8px">Top Drug Classes with Dose-Ranging Trials</h4>';
  html += '<div style="max-width:700px">';
  for (const [drug, trials] of topN) {
    const pct = (trials.length / maxCount * 100).toFixed(0);
    const doseSummary = trials.flatMap(t => t.doseLevels || []);
    const uniqueDoses = [...new Set(doseSummary)].sort((a, b) => a - b);
    const doseRange = uniqueDoses.length > 0
      ? uniqueDoses[0] + '-' + uniqueDoses[uniqueDoses.length - 1] + ' ' + (trials[0].doseUnit || '')
      : '';
    html += '<div style="margin-bottom:6px;display:flex;align-items:center;gap:8px">';
    html += '<div style="width:180px;font-size:0.8rem;text-align:right;overflow:hidden;text-overflow:ellipsis;white-space:nowrap" title="' + drug.replace(/"/g, '&quot;') + '">' + drug.slice(0, 30) + '</div>';
    html += '<div style="flex:1;height:20px;background:var(--bg);border-radius:4px;overflow:hidden">';
    html += '<div style="height:100%;width:' + pct + '%;background:linear-gradient(90deg,#7c3aed,#a78bfa);border-radius:4px;display:flex;align-items:center;padding-left:6px"><span style="font-size:0.7rem;color:#fff;font-weight:600">' + trials.length + '</span></div>';
    html += '</div>';
    html += '<div style="width:120px;font-size:0.7rem;color:var(--text-muted)">' + doseRange + '</div>';
    html += '</div>';
  }
  html += '</div>';

  // Top DR Opportunities table
  html += '<h4 style="font-size:0.9rem;margin:20px 0 8px">Top DR Review Opportunities</h4>';
  html += '<table style="font-size:0.8rem;border-collapse:collapse;width:100%;max-width:800px">';
  html += '<thead><tr style="border-bottom:2px solid var(--border)">';
  html += '<th style="padding:6px 10px;text-align:left">Drug Class</th>';
  html += '<th style="padding:6px 10px;text-align:right">DR Trials</th>';
  html += '<th style="padding:6px 10px;text-align:right">Total Enrollment</th>';
  html += '<th style="padding:6px 10px;text-align:left">Dose Range</th>';
  html += '<th style="padding:6px 10px;text-align:left">Phases</th>';
  html += '<th style="padding:6px 10px"></th>';
  html += '</tr></thead><tbody>';

  for (const [drug, trials] of topN.slice(0, 10)) {
    const enrollment = trials.reduce((s, t) => s + (parseInt(t.enrollment, 10) || 0), 0);
    const phases = [...new Set(trials.map(t => t.phase || '').filter(Boolean))].join(', ');
    const uniqueDoses = [...new Set(trials.flatMap(t => t.doseLevels || []))].sort((a, b) => a - b);
    const doseRange = uniqueDoses.length > 0 ? uniqueDoses[0] + '-' + uniqueDoses[uniqueDoses.length - 1] + ' ' + (trials[0].doseUnit || '') : '';
    const escapedDrug = drug.replace(/'/g, "\\'").replace(/"/g, '&quot;');

    html += '<tr style="border-bottom:1px solid var(--border)">';
    html += '<td style="padding:6px 10px">' + drug.slice(0, 40) + '</td>';
    html += '<td style="padding:6px 10px;text-align:right;font-weight:600;color:#7c3aed">' + trials.length + '</td>';
    html += '<td style="padding:6px 10px;text-align:right">' + (enrollment > 0 ? enrollment.toLocaleString() : '-') + '</td>';
    html += '<td style="padding:6px 10px">' + doseRange + '</td>';
    html += '<td style="padding:6px 10px">' + phases + '</td>';
    html += '<td style="padding:6px 10px"><button class="btn-sm btn-info" onclick="startReviewFromUniverse(\'' + escapedDrug + '\',currentGridSubcat?.id)">Start DR Review</button></td>';
    html += '</tr>';
  }
  html += '</tbody></table>';
  html += '</div>';

  el.innerHTML = html;
}
```

**Step 5: Commit**

```bash
git commit -m "feat: add DR Landscape visualization tab to Discover phase"
```

---

### Task 4: DR-Specific Search Strategy

**Files:**
- Modify: `metasprint-dose-response.html` — update `startReviewFromUniverse()` (~line 13545)

**Step 1: Add DR search terms to PICO population**

In `startReviewFromUniverse()`, after the PICO fields are populated, add dose-response search terms. Find where the PubMed search query is being built or where PICO fields are set. Add:

```javascript
// DR-specific search enhancement
const drTrials = (currentGridTrials || []).filter(t =>
  t.doseRanging && (t.interventions || []).some(iv =>
    (typeof iv === 'string' ? iv : (iv.name || '')).toLowerCase().includes(pickedIntervention.toLowerCase().slice(0, 10))
  )
);

if (drTrials.length > 0) {
  // Enhance search with DR MeSH and keywords
  const drSearchTerms = '"Dose-Response Relationship, Drug"[MeSH] OR dose-response OR dose-ranging OR dose-finding OR dose-titration OR dose-dependent';

  // If there's a search query field, append DR terms
  const searchInput = document.getElementById('searchQueryInput') || document.getElementById('pubmedQuery');
  if (searchInput) {
    const existing = searchInput.value || '';
    if (!existing.includes('dose-response')) {
      searchInput.value = existing ? '(' + existing + ') AND (' + drSearchTerms + ')' : drSearchTerms;
    }
  }

  // Update PICO intervention to include dose framing
  const picoI = document.getElementById('picoIntervention') || document.getElementById('pico-i');
  if (picoI && drTrials.length > 0) {
    const doseLevels = [...new Set(drTrials.flatMap(t => t.doseLevels || []))].sort((a, b) => a - b);
    const unit = drTrials[0].doseUnit || 'mg';
    if (doseLevels.length > 0) {
      picoI.value = pickedIntervention + ' at multiple doses (' + doseLevels.join(', ') + ' ' + unit + ')';
    }
  }

  // Update PICO comparator for dose-response
  const picoC = document.getElementById('picoComparator') || document.getElementById('pico-c');
  if (picoC) {
    picoC.value = 'placebo or lowest dose';
  }
}
```

**Step 2: Commit**

```bash
git commit -m "feat: add DR-specific search strategy to startReviewFromUniverse"
```

---

### Task 5: DR Screening Relevance Boost

**Files:**
- Modify: `metasprint-dose-response.html` — update screening scorer

**Step 1: Add DR signal detection**

Find the screening feature extraction (near line 5294 where `rctSignal`, `cardioSignal` are computed). Add a `drSignal` alongside:

```javascript
// DR signal detection
const DR_SIGNAL_TERMS = [
  /dose[\s-]?response/i,
  /dose[\s-]?rang/i,
  /dose[\s-]?find/i,
  /dose[\s-]?depend/i,
  /dose[\s-]?level/i,
  /dose[\s-]?titrat/i,
  /dose[\s-]?escala/i,
  /multiple\s+dos/i,
  /graded\s+dos/i,
  /\bmg\s*\/\s*day\b/i,
  /\bmg\s*\/\s*kg\b/i
];

function computeDRSignal(text) {
  if (!text) return 0;
  let matches = 0;
  for (const re of DR_SIGNAL_TERMS) {
    if (re.test(text)) matches++;
  }
  // Normalize: 0-1 scale
  if (matches >= 3) return 1.0;
  if (matches >= 2) return 0.7;
  if (matches >= 1) return 0.4;
  return 0;
}
```

**Step 2: Integrate DR signal into scoring**

Find where the final screening verdict is computed (the decision gate logic near line 5294-5320). Add `drSignal` to the features and use it as a boost:

```javascript
// In the feature extraction:
const drSignal = computeDRSignal(refText);

// In the scoring, add DR boost:
// If drSignal > 0, increase inclusion probability
if (drSignal >= 0.7) {
  // Strong DR signal: boost toward include
  features.drBoost = 0.15;
} else if (drSignal >= 0.4) {
  features.drBoost = 0.08;
} else {
  features.drBoost = 0;
}
```

**Step 3: Show DR badge in screening results**

In the screening results display (near line 4500-4512 where badges are shown), add a DR badge:

```javascript
// After the existing cardioSignal badge:
(autoScreenScores[r.id].drSignal > 0 ? ' <span style="color:#7c3aed">DR: ' + (autoScreenScores[r.id].drSignal * 100).toFixed(0) + '%</span>' : '') +
```

**Step 4: Add self-test**

```javascript
function _testDRSignal() {
  const results = [];
  results.push({ test: 'dose-response match', pass: computeDRSignal('A dose-response meta-analysis of X') >= 0.4 });
  results.push({ test: 'dose-ranging match', pass: computeDRSignal('This dose-ranging study evaluated 3 doses') >= 0.4 });
  results.push({ test: 'multi-term boost', pass: computeDRSignal('dose-response dose-finding dose-dependent study') >= 0.7 });
  results.push({ test: 'no match', pass: computeDRSignal('A randomized trial of drug X vs placebo') === 0 });
  results.push({ test: 'mg/kg match', pass: computeDRSignal('Patients received 0.5 mg/kg daily') >= 0.4 });
  console.table(results);
  return results.every(r => r.pass);
}
```

**Step 5: Commit**

```bash
git commit -m "feat: add DR relevance boost to screening scorer with DR signal detection"
```

---

### Task 6: Update Master Test Runner + Verification

**Files:**
- Modify: `metasprint-dose-response.html` — update `_runAllDRTests()`

**Step 1: Add new tests to master runner**

Update `_runAllDRTests()` to include the new test functions:

```javascript
try { results.doseRanging = _testDoseRanging(); } catch(e) { results.doseRanging = false; console.error('DoseRanging tests error:', e); }
try { results.drSignal = _testDRSignal(); } catch(e) { results.drSignal = false; console.error('DRSignal tests error:', e); }
```

**Step 2: Verify div balance**

```bash
grep -cE '<div[ >]' metasprint-dose-response.html
grep -c '</div>' metasprint-dose-response.html
```

**Step 3: Verify no literal `</script>`**

```bash
grep -n '</script>' metasprint-dose-response.html
```

**Step 4: Commit**

```bash
git commit -m "feat: update master test runner with DR detection and signal tests"
```

---

## Execution Order

| Task | Description | Depends On |
|------|-------------|------------|
| 1 | Dose-ranging trial detection | — |
| 2 | Universe grid DR badges + filter | Task 1 |
| 3 | DR Landscape visualization tab | Task 1 |
| 4 | DR-specific search strategy | Task 1 |
| 5 | DR screening relevance boost | — |
| 6 | Update tests + verification | Tasks 1-5 |

**Parallelizable:** Tasks 2, 3, 4 can run in parallel (all depend only on Task 1). Task 5 is independent.
