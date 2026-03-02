# Advanced DR Features Design

**Goal:** Add 4 features: study-overlay interactive curve, threshold/breakpoint detection, network dose-response, and R/Python code export.

## Feature 1: Study-Overlay + Interactive Curve
- Extend renderDoseResponseCurveSVG to always draw per-study polylines
- Hover tooltip (foreignObject) on data points: study ID, dose, effect, SE, weight
- Click study in legend to toggle visibility
- Per-study curves as thin colored lines behind the pooled curve

## Feature 2: Threshold/Breakpoint Detection
- detectDRBreakpoints(model, maxDose): numerical 1st/2nd derivatives
- Plateau: |dy/dx| < 5% of max slope for 3+ consecutive points
- Inflection: d2y/dx2 sign change
- Render as vertical dashed lines + labels on SVG

## Feature 3: Network Dose-Response
- renderNetworkDRCurve(drugGroups, modelType, confLevel)
- Fits same model independently per drug, overlays on shared axes
- Color-coded CI bands + legend
- Comparison table: drug, fit quality, key params

## Feature 4: R/Python Code Export
- generateRCode(points, model, method): complete dosresmeta2 script
- generatePythonCode(points, model, method): numpy/scipy WLS script
- Download buttons next to existing export panel
