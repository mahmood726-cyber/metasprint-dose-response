# Dose-Response Improvements Implementation Plan

> **For Claude:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to implement this plan task-by-task.

**Goal:** Add AIC model averaging weights, leave-one-out influence table, dose-response heterogeneity Q test (Q_DR), and CSV import for dose-response data.

**Architecture:** All changes are additive edits to `metasprint-dose-response.html` (21,199 lines). Each feature modifies 1-2 existing functions and adds 1 new function. Tests added to `test_dose_response_models.py`.

**Tech Stack:** Vanilla JS (single-file HTML app), Selenium + pytest for integration tests.

---

### Task 1: AIC Model Averaging Weights

**Files:**
- Modify: `metasprint-dose-response.html:10050-10071` (`compareDoseResponseModels`)
- Modify: `metasprint-dose-response.html:10079-10116` (`renderModelComparisonTable`)
- Modify: `metasprint-dose-response.html:10230-10248` (best model stats in `renderDoseResponseAnalysis`)
- Test: `test_dose_response_models.py`

**Step 1: Add AIC weight computation to `compareDoseResponseModels`**

In `compareDoseResponseModels()` at line ~10069, after `models.sort((a, b) => a.aic - b.aic)`, add weight computation:

```javascript
    // AIC model averaging weights (Burnham & Anderson 2002)
    var minAIC = models[0].aic;
    var sumExp = 0;
    for (var i = 0; i < models.length; i++) {
      models[i].deltaAIC = models[i].aic - minAIC;
      models[i].aicLikelihood = Math.exp(-models[i].deltaAIC / 2);
      sumExp += models[i].aicLikelihood;
    }
    for (var i = 0; i < models.length; i++) {
      models[i].aicWeight = sumExp > 0 ? models[i].aicLikelihood / sumExp : 0;
    }
```

Insert between the existing `models.sort(...)` (line 10069) and `return { best: models[0], all: models }` (line 10070).

**Step 2: Add weight column to `renderModelComparisonTable`**

In `renderModelComparisonTable()`, add a new `<th>` for "AIC Weight" in the header row (line ~10088, before the Parameters column):

```javascript
      '<th style="padding:6px 10px;text-align:right" title="Akaike weight: probability this is the best model">w<sub>i</sub></th>' +
```

And in the row loop (line ~10112, before the params `<td>`):

```javascript
        '<td style="padding:5px 10px;text-align:right">' + (m.aicWeight != null ? (m.aicWeight * 100).toFixed(1) + '%' : '\u2014') + '</td>' +
```

**Step 3: Add ΔAIC column too**

In the same header, add between AIC and BIC:

```javascript
      '<th style="padding:6px 10px;text-align:right" title="Difference from best model AIC">\u0394AIC</th>' +
```

And in the row:

```javascript
        '<td style="padding:5px 10px;text-align:right">' + (m.deltaAIC != null ? m.deltaAIC.toFixed(1) : '\u2014') + '</td>' +
```

**Step 4: Write Selenium test for AIC weights**

In `test_dose_response_models.py`, add to `TestModelFitting`:

```python
    def test_aic_weights_sum_to_one(self, driver):
        result = driver.execute_script("""
            var pts = [];
            for (var i = 0; i <= 10; i++) pts.push({dose: i, effect: Math.log(i+1)*2, se: 0.2});
            var c = compareDoseResponseModels(pts);
            if (!c) return null;
            var sumW = 0;
            for (var j = 0; j < c.all.length; j++) sumW += c.all[j].aicWeight;
            return {sumW: sumW, bestDelta: c.best.deltaAIC, bestWeight: c.best.aicWeight};
        """)
        assert result is not None
        assert abs(result['sumW'] - 1.0) < 0.001, f"AIC weights sum to {result['sumW']}, expected 1.0"
        assert result['bestDelta'] == 0, "Best model should have deltaAIC = 0"
        assert result['bestWeight'] > 0, "Best model weight must be > 0"
```

**Step 5: Run test to verify it passes**

Run: `python -m pytest test_dose_response_models.py::TestModelFitting::test_aic_weights_sum_to_one -v`
Expected: PASS

**Step 6: Commit**

```bash
git add metasprint-dose-response.html test_dose_response_models.py
git commit -m "feat: add AIC model averaging weights (Burnham & Anderson 2002)"
```

---

### Task 2: Leave-One-Out Influence Table

**Files:**
- Modify: `metasprint-dose-response.html:10121-10249` (`renderDoseResponseAnalysis`)
- Add: new function `renderLeaveOneOutDR()` after `renderModelComparisonTable` (~line 10116)
- Test: `test_dose_response_models.py`

**Step 1: Create `renderLeaveOneOutDR()` function**

Insert after `renderModelComparisonTable` (after line 10116):

```javascript
  /**
   * Render leave-one-out influence analysis for dose-response.
   * Re-fits the selected model k times, each time omitting one study.
   * @param {Array} points - dose-response data points [{dose, effect, se, studyId}]
   * @param {Object} fullModel - the model fitted on all data
   * @param {string} modelType - model name to refit (e.g. 'linear', 'rcs3')
   * @returns {string} HTML table
   */
  function renderLeaveOneOutDR(points, fullModel, modelType) {
    // Group points by study
    var studyIds = [];
    var seen = {};
    for (var i = 0; i < points.length; i++) {
      var sid = points[i].studyId ?? ('Row' + i);
      if (!seen[sid]) { seen[sid] = true; studyIds.push(sid); }
    }
    if (studyIds.length < 3) return ''; // Need at least 3 studies

    var fitFn;
    switch (modelType) {
      case 'Linear': fitFn = fitLinearDR; break;
      case 'Quadratic': fitFn = fitQuadraticDR; break;
      case 'Emax': fitFn = fitEmaxDR; break;
      default:
        if (modelType.startsWith('RCS-3')) fitFn = function(p) { return fitRCSDR(p, 3); };
        else if (modelType.startsWith('RCS-4')) fitFn = function(p) { return fitRCSDR(p, 4); };
        else if (modelType.startsWith('FP1')) fitFn = fitFP1DR;
        else if (modelType.startsWith('FP2')) fitFn = fitFP2DR;
        else fitFn = fitLinearDR;
    }

    var results = [];
    for (var s = 0; s < studyIds.length; s++) {
      var subset = points.filter(function(p) { return (p.studyId ?? ('Row' + points.indexOf(p))) !== studyIds[s]; });
      if (subset.length < 3) continue;
      var refitted = fitFn(subset);
      if (refitted) {
        results.push({
          omitted: studyIds[s],
          aic: refitted.aic,
          bic: refitted.bic,
          R2: refitted.R2,
          deltaAIC: refitted.aic - fullModel.aic,
          b1: refitted.b1 ?? null
        });
      }
    }

    if (results.length === 0) return '';

    var html = '<h3 style="font-size:0.95rem;margin:16px 0 8px">Leave-One-Out Influence Analysis</h3>';
    html += '<p style="font-size:0.78rem;color:var(--text-muted);margin-bottom:6px">Each row shows model fit with that study removed. Large \u0394AIC indicates influential study.</p>';
    html += '<table style="font-size:0.82rem;border-collapse:collapse;width:100%;max-width:700px">';
    html += '<thead><tr style="border-bottom:2px solid var(--border)">';
    html += '<th style="padding:6px 10px;text-align:left">Omitted Study</th>';
    html += '<th style="padding:6px 10px;text-align:right">AIC</th>';
    html += '<th style="padding:6px 10px;text-align:right">\u0394AIC</th>';
    html += '<th style="padding:6px 10px;text-align:right">R\u00B2</th>';
    if (results[0].b1 != null) {
      html += '<th style="padding:6px 10px;text-align:right">Slope (\u03B2\u2081)</th>';
    }
    html += '</tr></thead><tbody>';

    // Add "All studies" baseline row
    html += '<tr style="border-bottom:1px solid var(--border);background:rgba(37,99,235,0.05);font-weight:600">';
    html += '<td style="padding:5px 10px">(All studies)</td>';
    html += '<td style="padding:5px 10px;text-align:right">' + fullModel.aic.toFixed(1) + '</td>';
    html += '<td style="padding:5px 10px;text-align:right">\u2014</td>';
    html += '<td style="padding:5px 10px;text-align:right">' + (fullModel.R2 * 100).toFixed(1) + '%</td>';
    if (results[0].b1 != null) {
      html += '<td style="padding:5px 10px;text-align:right">' + (fullModel.b1 != null ? fullModel.b1.toFixed(4) : '\u2014') + '</td>';
    }
    html += '</tr>';

    for (var r = 0; r < results.length; r++) {
      var row = results[r];
      var influence = Math.abs(row.deltaAIC) > 2 ? 'color:#ef4444;font-weight:600' : '';
      html += '<tr style="border-bottom:1px solid var(--border)">';
      html += '<td style="padding:5px 10px">' + escapeHtml(row.omitted) + '</td>';
      html += '<td style="padding:5px 10px;text-align:right">' + row.aic.toFixed(1) + '</td>';
      html += '<td style="padding:5px 10px;text-align:right;' + influence + '">' + (row.deltaAIC > 0 ? '+' : '') + row.deltaAIC.toFixed(1) + '</td>';
      html += '<td style="padding:5px 10px;text-align:right">' + (row.R2 * 100).toFixed(1) + '%</td>';
      if (row.b1 != null) {
        html += '<td style="padding:5px 10px;text-align:right">' + row.b1.toFixed(4) + '</td>';
      }
      html += '</tr>';
    }
    html += '</tbody></table>';
    return html;
  }
```

**Step 2: Wire into `renderDoseResponseAnalysis`**

In the standard (non-one-stage) path at line ~10248, change:

```javascript
    el.innerHTML = '<h3 ...' + svgHtml + statsHtml + tableHtml;
```

to:

```javascript
    // Leave-one-out influence
    const looHtml = renderLeaveOneOutDR(data.points, best, best.model);
    el.innerHTML = '<h3 ...' + svgHtml + statsHtml + tableHtml + looHtml;
```

Also wire it into the one-stage path at line ~10193-10194:

```javascript
    const looHtml = renderLeaveOneOutDR(data.points, oneStageResult, oneStageResult.model);
    el.innerHTML += tableHtml + looHtml;
```

**Step 3: Write Selenium test**

```python
    def test_leave_one_out_dr(self, driver):
        """renderLeaveOneOutDR should produce a table with one row per study."""
        result = driver.execute_script("""
            var pts = [];
            var studies = ['StudyA','StudyB','StudyC','StudyD'];
            for (var s = 0; s < studies.length; s++) {
                for (var d = 0; d <= 4; d++) {
                    pts.push({dose: d, effect: 1 + 0.3*d + (s-1.5)*0.2, se: 0.15, studyId: studies[s]});
                }
            }
            var model = fitLinearDR(pts);
            if (!model) return null;
            var html = renderLeaveOneOutDR(pts, model, 'Linear');
            return {
                hasTable: html.includes('<table'),
                hasAllStudies: html.includes('All studies'),
                hasStudyA: html.includes('StudyA'),
                length: html.length
            };
        """)
        assert result is not None
        assert result['hasTable'] is True
        assert result['hasAllStudies'] is True
        assert result['hasStudyA'] is True
```

**Step 4: Run test**

Run: `python -m pytest test_dose_response_models.py::TestModelFitting::test_leave_one_out_dr -v`
Expected: PASS

**Step 5: Commit**

```bash
git add metasprint-dose-response.html test_dose_response_models.py
git commit -m "feat: add leave-one-out influence table for dose-response models"
```

---

### Task 3: Dose-Response Heterogeneity Q Test (Q_DR)

**Files:**
- Modify: `metasprint-dose-response.html` — add `computeDRHeterogeneity()` after model fitting functions (~line 10050)
- Modify: `metasprint-dose-response.html:10230-10248` (stats section in `renderDoseResponseAnalysis`)
- Test: `test_dose_response_models.py`

**Step 1: Create `computeDRHeterogeneity()` function**

Insert before `compareDoseResponseModels` (before line 10050):

```javascript
  /**
   * Compute dose-response heterogeneity Q test.
   * Tests H0: all studies follow the same dose-response relationship.
   * Uses Cochran-type Q on study-level residuals from the pooled DR model.
   *
   * @param {Array} points - dose-response data points [{dose, effect, se, studyId}]
   * @param {Object} model - fitted model with predict(dose) method or coefficients
   * @returns {Object} {Q_DR, df, pValue, I2_DR}
   */
  function computeDRHeterogeneity(points, model) {
    if (!points || points.length < 3 || !model) return null;

    // Get residuals from model
    var residuals = model.residuals;
    if (!residuals || residuals.length !== points.length) return null;

    // Weighted sum of squared residuals
    var Q = 0;
    var sumW = 0;
    for (var i = 0; i < points.length; i++) {
      var w = 1 / (points[i].se * points[i].se);
      if (!isFinite(w) || w <= 0) continue;
      Q += w * residuals[i] * residuals[i];
      sumW += w;
    }

    // Degrees of freedom = n_points - n_parameters
    var nParams = 1; // intercept counted separately if exists
    if (model.model === 'Linear') nParams = 1; // b1
    else if (model.model === 'Quadratic') nParams = 2; // b1, b2
    else if (model.model === 'Emax') nParams = 2; // Emax, ED50
    else if (model.model && model.model.startsWith('RCS-3')) nParams = 2;
    else if (model.model && model.model.startsWith('RCS-4')) nParams = 3;
    else if (model.model && model.model.startsWith('FP1')) nParams = 1;
    else if (model.model && model.model.startsWith('FP2')) nParams = 2;

    var df = points.length - nParams;
    if (df <= 0) return null;

    // p-value from chi-squared distribution
    var pValue = 1 - chi2CDF(Q, df);

    // I²_DR: proportion of variability beyond sampling error
    var I2_DR = df > 0 && Q > df ? ((Q - df) / Q) * 100 : 0;

    return { Q_DR: Q, df: df, pValue: pValue, I2_DR: I2_DR };
  }
```

Note: `chi2CDF` should already exist (used by non-linearity Wald test). Verify with a grep. If not, use the Wilson-Hilferty approximation:

```javascript
  function chi2CDF(x, df) {
    if (x <= 0 || df <= 0) return 0;
    // Wilson-Hilferty approximation
    var z = Math.pow(x / df, 1/3) - (1 - 2 / (9 * df));
    z /= Math.sqrt(2 / (9 * df));
    return normalCDF(z);
  }
```

**Step 2: Wire into `renderDoseResponseAnalysis` stats section**

At line ~10247, before the closing `</div>` of `statsHtml`, add:

```javascript
    // Dose-response heterogeneity
    var drHet = computeDRHeterogeneity(data.points, best);
    if (drHet) {
      statsHtml += '<br><strong>DR Heterogeneity:</strong> Q<sub>DR</sub> = ' + drHet.Q_DR.toFixed(2) +
        ' (df=' + drHet.df + ', p=' + (drHet.pValue < 0.001 ? '<0.001' : drHet.pValue.toFixed(3)) + ')' +
        ' | I\u00B2<sub>DR</sub> = ' + drHet.I2_DR.toFixed(1) + '%';
    }
```

**Step 3: Write Selenium test**

```python
    def test_dr_heterogeneity(self, driver):
        """computeDRHeterogeneity should return Q, df, pValue, I2_DR."""
        result = driver.execute_script("""
            var pts = [];
            for (var i = 0; i <= 8; i++) pts.push({dose: i, effect: 1 + 0.5*i + (Math.random()-0.5)*0.3, se: 0.15});
            var model = fitLinearDR(pts);
            if (!model) return null;
            var het = computeDRHeterogeneity(pts, model);
            return het;
        """)
        assert result is not None
        assert 'Q_DR' in result
        assert result['Q_DR'] >= 0
        assert result['df'] > 0
        assert 0 <= result['pValue'] <= 1
        assert 0 <= result['I2_DR'] <= 100
```

**Step 4: Run test**

Run: `python -m pytest test_dose_response_models.py -k test_dr_heterogeneity -v`
Expected: PASS

**Step 5: Commit**

```bash
git add metasprint-dose-response.html test_dose_response_models.py
git commit -m "feat: add dose-response heterogeneity Q test (Q_DR, I²_DR)"
```

---

### Task 4: CSV Import for Dose-Response Data

**Files:**
- Modify: `metasprint-dose-response.html:1368` (add Import CSV button next to Export)
- Add: new function `importDoseResponseCSV()` after `exportStudiesCSV()` (~line 7220)
- Add: hidden `<input type="file">` element for CSV upload
- Test: `test_dose_response_models.py`

**Step 1: Add Import CSV button and file input to Extract phase**

At line ~1368, after the Export buttons, add:

```html
        <button class="btn-outline" onclick="document.getElementById('csvImportInput').click()">Import CSV</button>
        <input type="file" id="csvImportInput" accept=".csv,.tsv,.txt" style="display:none" onchange="importDoseResponseCSV(this)">
```

**Step 2: Create `importDoseResponseCSV()` function**

Insert after `exportStudiesCSV()` (after line ~7220):

```javascript
  /**
   * Import dose-response data from CSV file into extractedStudies.
   * Expects header row matching export format:
   *   Study ID, Trial ID, NCT ID, PMID, DOI, Outcome, Timepoint, Population,
   *   Dose, Dose Unit, Reference Dose, N, Effect, Lower CI, Upper CI, SE, Type, Subgroup, Notes
   * Also accepts minimal format: Study, Dose, N, Effect, LowerCI, UpperCI
   * @param {HTMLInputElement} input - file input element
   */
  function importDoseResponseCSV(input) {
    if (!input.files || !input.files[0]) return;
    var file = input.files[0];
    var reader = new FileReader();
    reader.onload = function(e) {
      var text = e.target.result;
      var lines = text.split(/\r?\n/).filter(function(l) { return l.trim().length > 0; });
      if (lines.length < 2) {
        showToast('CSV must have a header row and at least one data row', 'warning');
        input.value = '';
        return;
      }

      // Parse header — detect delimiter (comma or tab)
      var delim = lines[0].includes('\t') ? '\t' : ',';
      var header = parseCSVLine(lines[0], delim).map(function(h) { return h.trim().toLowerCase().replace(/[\s_-]+/g, ''); });

      // Map header columns to field names (flexible matching)
      var colMap = {};
      var fieldAliases = {
        studyid: 'authorYear', study: 'authorYear', author: 'authorYear', studyname: 'authorYear',
        trialid: 'trialId', trial: 'trialId',
        nctid: 'nctId', nct: 'nctId',
        pmid: 'pmid',
        doi: 'doi',
        outcome: 'outcomeId',
        timepoint: 'timepoint',
        population: 'analysisPopulation',
        dose: 'dose',
        doseunit: 'doseUnit', unit: 'doseUnit',
        referencedose: 'referenceDose', refdose: 'referenceDose',
        n: 'nTotal', samplesize: 'nTotal', ntotal: 'nTotal',
        effect: 'effectEstimate', estimate: 'effectEstimate', effectestimate: 'effectEstimate',
        lowerci: 'lowerCI', lower: 'lowerCI', cilo: 'lowerCI', lo: 'lowerCI', cilower: 'lowerCI',
        upperci: 'upperCI', upper: 'upperCI', cihi: 'upperCI', hi: 'upperCI', ciupper: 'upperCI',
        se: 'se', stderr: 'se', standarderror: 'se',
        type: 'effectType', effecttype: 'effectType',
        subgroup: 'subgroup',
        notes: 'notes',
        slope: 'slope',
        slopese: 'slopeSE',
        doserange: 'doseRange',
        kdoses: 'kDoses',
        mean: 'mean',
        sd: 'sd'
      };

      for (var c = 0; c < header.length; c++) {
        if (fieldAliases[header[c]]) {
          colMap[fieldAliases[header[c]]] = c;
        }
      }

      // Must have at minimum: Study + Dose + Effect
      if (colMap.authorYear == null || colMap.dose == null || (colMap.effectEstimate == null && colMap.slope == null && colMap.mean == null)) {
        showToast('CSV must have at least: Study, Dose, and Effect (or Slope, or Mean) columns', 'warning');
        input.value = '';
        return;
      }

      // Detect input mode
      if (colMap.slope != null && colMap.slopeSE != null) {
        setInputMode('2x2');
      } else if (colMap.mean != null && colMap.sd != null) {
        setInputMode('continuous');
      } else {
        setInputMode('effect');
      }

      // Parse rows
      var imported = [];
      for (var r = 1; r < lines.length; r++) {
        var vals = parseCSVLine(lines[r], delim);
        if (vals.length < 3) continue;

        var row = {
          authorYear: getCSVVal(vals, colMap.authorYear, ''),
          trialId: getCSVVal(vals, colMap.trialId, ''),
          nctId: getCSVVal(vals, colMap.nctId, ''),
          pmid: getCSVVal(vals, colMap.pmid, ''),
          doi: getCSVVal(vals, colMap.doi, ''),
          outcomeId: getCSVVal(vals, colMap.outcomeId, ''),
          timepoint: getCSVVal(vals, colMap.timepoint, ''),
          analysisPopulation: getCSVVal(vals, colMap.analysisPopulation, ''),
          dose: parseFloatSafe(getCSVVal(vals, colMap.dose, '')),
          doseUnit: getCSVVal(vals, colMap.doseUnit, ''),
          referenceDose: parseFloatSafe(getCSVVal(vals, colMap.referenceDose, '')),
          nTotal: parseIntSafe(getCSVVal(vals, colMap.nTotal, '')),
          effectEstimate: parseFloatSafe(getCSVVal(vals, colMap.effectEstimate, '')),
          lowerCI: parseFloatSafe(getCSVVal(vals, colMap.lowerCI, '')),
          upperCI: parseFloatSafe(getCSVVal(vals, colMap.upperCI, '')),
          se: parseFloatSafe(getCSVVal(vals, colMap.se, '')),
          effectType: getCSVVal(vals, colMap.effectType, 'MD'),
          subgroup: getCSVVal(vals, colMap.subgroup, ''),
          notes: getCSVVal(vals, colMap.notes, ''),
          slope: parseFloatSafe(getCSVVal(vals, colMap.slope, '')),
          slopeSE: parseFloatSafe(getCSVVal(vals, colMap.slopeSE, '')),
          doseRange: getCSVVal(vals, colMap.doseRange, ''),
          kDoses: parseIntSafe(getCSVVal(vals, colMap.kDoses, '')),
          mean: parseFloatSafe(getCSVVal(vals, colMap.mean, '')),
          sd: parseFloatSafe(getCSVVal(vals, colMap.sd, '')),
          verificationStatus: 'imported'
        };

        imported.push(row);
      }

      if (imported.length === 0) {
        showToast('No valid data rows found in CSV', 'warning');
        input.value = '';
        return;
      }

      // Append to existing studies (don't overwrite)
      extractedStudies = extractedStudies.concat(imported);
      renderExtractTable();
      saveExtractedStudies();
      showToast('Imported ' + imported.length + ' rows from CSV', 'success');
      input.value = '';
    };
    reader.readAsText(file);
  }

  function parseCSVLine(line, delim) {
    // Handle quoted fields with embedded delimiters
    var result = [];
    var current = '';
    var inQuotes = false;
    for (var i = 0; i < line.length; i++) {
      var ch = line[i];
      if (ch === '"') {
        if (inQuotes && i + 1 < line.length && line[i + 1] === '"') {
          current += '"';
          i++;
        } else {
          inQuotes = !inQuotes;
        }
      } else if (ch === delim && !inQuotes) {
        result.push(current);
        current = '';
      } else {
        current += ch;
      }
    }
    result.push(current);
    return result;
  }

  function getCSVVal(vals, idx, fallback) {
    if (idx == null || idx >= vals.length) return fallback;
    var v = vals[idx].trim().replace(/^"+|"+$/g, '');
    return v || fallback;
  }

  function parseFloatSafe(s) {
    if (s == null || s === '') return null;
    var v = parseFloat(s);
    return isFinite(v) ? v : null;
  }

  function parseIntSafe(s) {
    if (s == null || s === '') return null;
    var v = parseInt(s, 10);
    return isFinite(v) ? v : null;
  }
```

**Step 3: Write Selenium test**

```python
    def test_csv_import_parse(self, driver):
        """importDoseResponseCSV parser helpers should work correctly."""
        result = driver.execute_script("""
            // Test parseCSVLine
            var line1 = parseCSVLine('StudyA,10,mg,100,1.5,0.8,2.2,0.3,MD', ',');
            // Test with quoted field
            var line2 = parseCSVLine('"Study, B",20,mg,150,2.0,1.0,3.0,0.4,MD', ',');
            // Test parseFloatSafe
            var f1 = parseFloatSafe('3.14');
            var f2 = parseFloatSafe('');
            var f3 = parseFloatSafe('abc');
            var f4 = parseFloatSafe('0');
            return {
                line1Len: line1.length, line1Study: line1[0],
                line2Study: line2[0],
                f1: f1, f2: f2, f3: f3, f4: f4
            };
        """)
        assert result is not None
        assert result['line1Len'] == 9
        assert result['line1Study'] == 'StudyA'
        assert result['line2Study'] == 'Study, B'  # Quoted field preserved
        assert result['f1'] == 3.14
        assert result['f2'] is None
        assert result['f3'] is None
        assert result['f4'] == 0  # Must not drop valid zero
```

**Step 4: Run test**

Run: `python -m pytest test_dose_response_models.py -k test_csv_import_parse -v`
Expected: PASS

**Step 5: Commit**

```bash
git add metasprint-dose-response.html test_dose_response_models.py
git commit -m "feat: add CSV import for dose-response data with flexible column mapping"
```

---

### Task 5: Inline Self-Tests + Final Verification

**Files:**
- Modify: `metasprint-dose-response.html` — add `_testAICWeights()`, `_testDRHeterogeneity()`, `_testCSVParser()` to `_runAllDRTests()`
- Test: `test_dose_response_models.py` — run all 30+ tests

**Step 1: Add inline self-tests**

After the existing `_testDRSignal()` function (~line 10006), add:

```javascript
  function _testAICWeights() {
    var pts = [];
    for (var i = 0; i <= 8; i++) pts.push({dose: i, effect: 1 + 0.5*i, se: 0.15});
    var c = compareDoseResponseModels(pts);
    if (!c || !c.all || c.all.length === 0) return false;
    // All models must have aicWeight
    var sumW = 0;
    for (var j = 0; j < c.all.length; j++) {
      if (c.all[j].aicWeight == null) return false;
      sumW += c.all[j].aicWeight;
    }
    // Weights must sum to ~1
    if (Math.abs(sumW - 1.0) > 0.01) return false;
    // Best model must have deltaAIC = 0
    if (c.best.deltaAIC !== 0) return false;
    return true;
  }

  function _testDRHeterogeneity() {
    var pts = [];
    for (var i = 0; i <= 6; i++) pts.push({dose: i, effect: 1 + 0.5*i, se: 0.15});
    var model = fitLinearDR(pts);
    if (!model) return false;
    var het = computeDRHeterogeneity(pts, model);
    if (!het) return false;
    if (het.Q_DR < 0 || het.df <= 0) return false;
    if (het.pValue < 0 || het.pValue > 1) return false;
    if (het.I2_DR < 0 || het.I2_DR > 100) return false;
    return true;
  }

  function _testCSVParser() {
    var line = parseCSVLine('a,"b,c",d', ',');
    if (line.length !== 3 || line[1] !== 'b,c') return false;
    if (parseFloatSafe('0') !== 0) return false; // Must not drop zero
    if (parseFloatSafe('') !== null) return false;
    if (parseFloatSafe('abc') !== null) return false;
    if (parseIntSafe('42') !== 42) return false;
    return true;
  }
```

**Step 2: Register in `_runAllDRTests()`**

In the existing `_runAllDRTests()` function, add:

```javascript
    try { results.aicWeights = _testAICWeights(); } catch(e) { results.aicWeights = false; }
    try { results.drHeterogeneity = _testDRHeterogeneity(); } catch(e) { results.drHeterogeneity = false; }
    try { results.csvParser = _testCSVParser(); } catch(e) { results.csvParser = false; }
```

**Step 3: Run all Selenium tests**

Run: `python -m pytest test_dose_response_models.py -v`
Expected: All 34+ tests PASS

**Step 4: Verify div balance**

```bash
grep -cP '<div[\s>]' metasprint-dose-response.html
grep -c '</div>' metasprint-dose-response.html
```

Expected: similar to baseline ratio (off by ~22 due to JS-constructed divs).

**Step 5: Commit**

```bash
git add metasprint-dose-response.html test_dose_response_models.py
git commit -m "test: add inline self-tests for AIC weights, Q_DR, CSV parser (34+ tests pass)"
```
