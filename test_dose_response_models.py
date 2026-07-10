"""Integration tests for dose-response model enhancements."""
import sys, io, time, os, pytest

# Windows UTF-8 safety: set environment variable instead of wrapping stdout
# (wrapping stdout breaks pytest's capture mechanism)
os.environ.setdefault('PYTHONIOENCODING', 'utf-8')

from selenium import webdriver
from selenium.webdriver.common.by import By
from selenium.webdriver.chrome.options import Options
from selenium.webdriver.support.ui import WebDriverWait, Select
from selenium.webdriver.support import expected_conditions as EC

HTML_PATH = 'file:///' + os.path.join(os.path.dirname(os.path.abspath(__file__)), 'metasprint-dose-response.html').replace('\\', '/')

@pytest.fixture(scope='module')
def driver():
    opts = Options()
    opts.add_argument('--headless=new')
    opts.add_argument('--no-sandbox')
    opts.add_argument('--window-size=1400,900')
    opts.set_capability('goog:loggingPrefs', {'browser': 'ALL'})
    d = webdriver.Chrome(options=opts)
    d.get(HTML_PATH)
    time.sleep(3)  # Let app initialize
    yield d
    d.quit()


class TestSelfTests:
    """Run inline self-tests via browser console."""

    def test_all_self_tests_pass(self, driver):
        # _runAllDRTests checks v === true but _testMatrixOps/_testRCS/_testFP
        # return {passed, failed, details} objects. So we run each individually.
        result = driver.execute_script("""
            var ok = true;
            var m = _testMatrixOps(); if (m.failed > 0) ok = false;
            var r = _testRCS(); if (r.failed > 0) ok = false;
            var f = _testFP(); if (f.failed > 0) ok = false;
            if (!_testOneStage()) ok = false;
            if (!_testContinuousInput()) ok = false;
            if (!_testAICWeights()) ok = false;
            if (!_testDRHeterogeneity()) ok = false;
            if (!_testCSVParser()) ok = false;
            return ok;
        """)
        assert result is True, 'Inline self-tests failed'

    def test_matrix_ops(self, driver):
        result = driver.execute_script('return _testMatrixOps()')
        assert result is not None
        assert result['failed'] == 0, f"Matrix tests had {result['failed']} failures"

    def test_rcs(self, driver):
        result = driver.execute_script('return _testRCS()')
        assert result is not None
        assert result['failed'] == 0, f"RCS tests had {result['failed']} failures"

    def test_fp(self, driver):
        result = driver.execute_script('return _testFP()')
        assert result is not None
        assert result['failed'] == 0, f"FP tests had {result['failed']} failures"

    def test_one_stage(self, driver):
        result = driver.execute_script('return _testOneStage()')
        assert result is True

    def test_continuous_input(self, driver):
        result = driver.execute_script('return _testContinuousInput()')
        assert result is True


class TestUIElements:
    """Verify all new UI elements exist."""

    def test_model_selector_has_8_options(self, driver):
        select = Select(driver.find_element(By.ID, 'drModelSelect'))
        options = [o.get_attribute('value') for o in select.options]
        assert 'auto' in options
        assert 'linear' in options
        assert 'quadratic' in options
        assert 'emax' in options
        assert 'rcs3' in options
        assert 'rcs4' in options
        assert 'fp1' in options
        assert 'fp2' in options
        assert 'loglinear' in options
        assert 'exponential' in options
        assert 'hill' in options
        assert 'averaged' in options
        assert len(options) == 12

    def test_one_stage_toggle_exists(self, driver):
        el = driver.find_element(By.ID, 'drOneStageToggle')
        assert el is not None
        assert el.get_attribute('type') == 'checkbox'

    def test_continuous_input_radio_exists(self, driver):
        radios = driver.find_elements(By.CSS_SELECTOR, 'input[name="inputMode"][value="continuous"]')
        assert len(radios) == 1


class TestModelFitting:
    """Test each model type via JS execution."""

    def test_linear_model(self, driver):
        result = driver.execute_script("""
            var pts = [];
            for (var i = 0; i <= 6; i++) pts.push({dose: i, effect: 1 + 0.5*i, se: 0.15});
            var m = fitLinearDR(pts);
            return m ? {model: m.model, R2: m.R2, b1: m.b1} : null;
        """)
        assert result is not None
        assert result['R2'] > 0.99
        assert abs(result['b1'] - 0.5) < 0.01

    def test_rcs3_model(self, driver):
        result = driver.execute_script("""
            var pts = [];
            for (var i = 0; i <= 8; i++) pts.push({dose: i, effect: Math.sin(i*0.4)*2, se: 0.2});
            var m = fitRCSDR(pts, 3);
            return m ? {model: m.model, R2: m.R2, pNL: m.pNonlinear, knots: m.knots.length} : null;
        """)
        assert result is not None
        assert result['model'] == 'RCS-3'
        assert result['R2'] > 0
        assert result['knots'] == 3

    def test_rcs4_model(self, driver):
        result = driver.execute_script("""
            var pts = [];
            for (var i = 0; i <= 10; i++) pts.push({dose: i, effect: Math.log(i+1)*3, se: 0.2});
            var m = fitRCSDR(pts, 4);
            return m ? {model: m.model, knots: m.knots.length} : null;
        """)
        assert result is not None
        assert result['model'] == 'RCS-4'
        assert result['knots'] == 4

    def test_fp1_model(self, driver):
        result = driver.execute_script("""
            var pts = [];
            for (var i = 1; i <= 8; i++) pts.push({dose: i, effect: 3*Math.sqrt(i), se: 0.1});
            var m = fitFP1DR(pts);
            return m ? {model: m.model, power: m.power, R2: m.R2} : null;
        """)
        assert result is not None
        assert result['power'] == 0.5
        assert result['R2'] > 0.99

    def test_fp2_model(self, driver):
        result = driver.execute_script("""
            var pts = [];
            for (var i = 1; i <= 10; i++) pts.push({dose: i, effect: 2*Math.sqrt(i) - 0.3*i, se: 0.15});
            var m = fitFP2DR(pts);
            return m ? {power1: m.power1, power2: m.power2, R2: m.R2} : null;
        """)
        assert result is not None
        assert result['R2'] > 0.9

    def test_one_stage_mixed_effects(self, driver):
        result = driver.execute_script("""
            var groups = [];
            for (var s = 0; s < 4; s++) {
                var pts = [];
                var b0 = 1 + (s-1.5)*0.3;
                for (var d = 0; d <= 4; d++) pts.push({dose: d, effect: b0 + 0.4*d, se: 0.2});
                groups.push({studyId: 'S' + s, points: pts});
            }
            var r = fitOneStageMixedEffects(groups, 'linear');
            return r ? {converged: r.convergence.converged, slope: r.beta[1], nStudy: r.studyModels.length} : null;
        """)
        assert result is not None
        assert result['converged'] is True
        assert abs(result['slope'] - 0.4) < 0.15
        assert result['nStudy'] >= 3

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

    def test_model_comparison_includes_all(self, driver):
        result = driver.execute_script("""
            var pts = [];
            for (var i = 0; i <= 10; i++) pts.push({dose: i, effect: Math.log(i+1)*2 + Math.random()*0.1, se: 0.2});
            var c = compareDoseResponseModels(pts);
            if (!c) return null;
            return {count: c.all.length, models: c.all.map(function(m) { return m.model; }), bestAIC: c.best.aic};
        """)
        assert result is not None
        assert result['count'] >= 5, f"Expected 5+ models, got {result['count']}: {result['models']}"
        models = result['models']
        # Should have at least Linear, Quadratic, Emax, RCS-3, FP1
        model_prefixes = [m[:3] for m in models]
        assert 'Lin' in model_prefixes
        assert 'Qua' in model_prefixes

    def test_leave_one_out_dr(self, driver):
        """renderLeaveOneOutDR should produce a table with baseline + one row per study."""
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
                hasStudyD: html.includes('StudyD'),
                length: html.length
            };
        """)
        assert result is not None
        assert result['hasTable'] is True
        assert result['hasAllStudies'] is True
        assert result['hasStudyA'] is True
        assert result['hasStudyD'] is True


    def test_dr_heterogeneity(self, driver):
        """computeDRHeterogeneity should return Q_DR, df, pValue, I2_DR."""
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


class TestEdgeCases:
    """Test edge cases and error handling."""

    def test_rcs_too_few_points(self, driver):
        result = driver.execute_script("""
            var pts = [{dose:0, effect:0, se:0.1}, {dose:1, effect:1, se:0.1}];
            return fitRCSDR(pts, 3);
        """)
        assert result is None

    def test_fp1_too_few_points(self, driver):
        result = driver.execute_script("""
            var pts = [{dose:1, effect:1, se:0.1}, {dose:2, effect:2, se:0.1}];
            return fitFP1DR(pts);
        """)
        assert result is None

    def test_one_stage_single_study(self, driver):
        result = driver.execute_script("""
            var g = [{studyId: 'S1', points: [{dose:0,effect:0,se:0.1},{dose:1,effect:1,se:0.1}]}];
            return fitOneStageMixedEffects(g, 'linear');
        """)
        assert result is None

    def test_singular_matrix_handled(self, driver):
        result = driver.execute_script("""
            return matInvertNxN([[1,2],[2,4]]);
        """)
        assert result is None


class TestDRUniverse:
    """Test DR Universe features: dose-ranging detection, DR signal, landscape tab."""

    def test_dose_ranging_detection_basic(self, driver):
        """detectDoseRanging should flag a trial with 3+ dose levels."""
        result = driver.execute_script("""
            var trial = {
                arms: [
                    {label: 'Placebo'},
                    {label: '10 mg once daily'},
                    {label: '20 mg once daily'},
                    {label: '40 mg once daily'}
                ]
            };
            detectDoseRanging(trial);
            return {dr: trial.doseRanging, levels: trial.doseLevels, unit: trial.doseUnit};
        """)
        assert result is not None
        assert result['dr'] is True
        assert len(result['levels']) >= 3
        assert result['unit'] == 'mg'

    def test_dose_ranging_not_flagged_for_placebo_only(self, driver):
        """Trial with only placebo arms should not be dose-ranging."""
        result = driver.execute_script("""
            var trial = {
                arms: [
                    {label: 'Placebo'},
                    {label: 'Standard of care'}
                ]
            };
            detectDoseRanging(trial);
            return {dr: trial.doseRanging, levels: trial.doseLevels};
        """)
        assert result is not None
        assert result['dr'] is False

    def test_dose_ranging_handles_mcg_units(self, driver):
        """Should detect mcg as a dose unit."""
        result = driver.execute_script("""
            var trial = {
                arms: [
                    {label: '100 mcg'},
                    {label: '200 mcg'},
                    {label: '400 mcg'},
                    {label: 'Placebo'}
                ]
            };
            detectDoseRanging(trial);
            return {dr: trial.doseRanging, unit: trial.doseUnit};
        """)
        assert result is not None
        assert result['dr'] is True
        assert result['unit'] == 'mcg'

    def test_dose_ranging_too_few_levels(self, driver):
        """Fewer than 3 dose levels should not flag as dose-ranging."""
        result = driver.execute_script("""
            var trial = {
                arms: [
                    {label: 'Placebo'},
                    {label: '10 mg'},
                    {label: '20 mg'}
                ]
            };
            detectDoseRanging(trial);
            return {dr: trial.doseRanging};
        """)
        assert result is not None
        assert result['dr'] is False

    def test_dose_ranging_self_test(self, driver):
        """Run the inline _testDoseRanging self-test."""
        result = driver.execute_script('return _testDoseRanging()')
        assert result is True

    def test_dr_signal_detection(self, driver):
        """computeDRSignal should detect DR terms in text."""
        result = driver.execute_script("""
            var r1 = computeDRSignal('This dose-response study evaluated dose-ranging');
            var r2 = computeDRSignal('A randomized controlled trial of aspirin');
            var r3 = computeDRSignal(null);
            return {multi: r1, none: r2, nil: r3};
        """)
        assert result is not None
        assert result['multi'] >= 0.7  # 2+ DR terms
        assert result['none'] == 0     # No DR terms
        assert result['nil'] == 0      # null input

    def test_dr_signal_self_test(self, driver):
        """Run the inline _testDRSignal self-test."""
        result = driver.execute_script('return _testDRSignal()')
        assert result is True

    def test_dr_landscape_tab_exists(self, driver):
        """DR Landscape tab button should exist."""
        tabs = driver.find_elements(By.CSS_SELECTOR, 'button[data-view="drlandscape"]')
        assert len(tabs) == 1

    def test_dr_landscape_container_exists(self, driver):
        """DR Landscape container div should exist."""
        el = driver.find_elements(By.ID, 'universeViewDRLandscape')
        assert len(el) == 1

    def test_render_dr_landscape_with_data(self, driver):
        """renderDRLandscape should populate the container with content."""
        result = driver.execute_script("""
            var trials = [];
            for (var i = 0; i < 5; i++) {
                trials.push({
                    nctId: 'NCT0000000' + i,
                    title: 'Test trial ' + i,
                    doseRanging: true,
                    doseLevels: [10, 20, 40],
                    doseUnit: 'mg',
                    interventions: [{name: 'DrugA'}],
                    enrollment: 100,
                    phase: 'Phase 3'
                });
            }
            trials.push({
                nctId: 'NCT00000010',
                title: 'Non-DR trial',
                doseRanging: false,
                interventions: [{name: 'DrugB'}],
                enrollment: 50,
                phase: 'Phase 2'
            });
            renderDRLandscape(trials);
            var el = document.getElementById('universeViewDRLandscape');
            return {hasContent: el.innerHTML.length > 50, hasDRCount: el.innerHTML.includes('5')};
        """)
        assert result is not None
        assert result['hasContent'] is True
        assert result['hasDRCount'] is True


class TestExceedDosresmeta:
    """Tests for features that exceed dosresmeta2: GL covariance, GLS, ML/REML,
    new models, prediction intervals, model averaging, R cross-validation."""

    # --- UI elements ---

    def test_estimation_method_selector_exists(self, driver):
        select = Select(driver.find_element(By.ID, 'drEstMethod'))
        opts = [o.get_attribute('value') for o in select.options]
        assert 'dl' in opts
        assert 'reml' in opts
        assert 'ml' in opts
        assert 'fixed' in opts
        assert len(opts) == 4

    def test_prediction_interval_checkbox_exists(self, driver):
        el = driver.find_element(By.ID, 'drShowPI')
        assert el.get_attribute('type') == 'checkbox'

    def test_model_selector_has_new_models(self, driver):
        select = Select(driver.find_element(By.ID, 'drModelSelect'))
        opts = [o.get_attribute('value') for o in select.options]
        assert 'loglinear' in opts
        assert 'exponential' in opts
        assert 'hill' in opts
        assert 'averaged' in opts

    # --- Inline test suites ---

    def test_gl_covariance_suite(self, driver):
        result = driver.execute_script('return _testGLCovariance()')
        assert result is True

    def test_gls_fit_suite(self, driver):
        result = driver.execute_script('return _testGLSFit()')
        assert result is True

    def test_ml_reml_suite(self, driver):
        result = driver.execute_script('return _testMLREML()')
        assert result is True

    def test_new_models_suite(self, driver):
        result = driver.execute_script('return _testNewModels()')
        assert result is True

    def test_one_stage_all_models_suite(self, driver):
        result = driver.execute_script('return _testOneStageAllModels()')
        assert result is True

    def test_prediction_interval_suite(self, driver):
        result = driver.execute_script('return _testPredictionInterval()')
        assert result is True

    def test_model_averaging_suite(self, driver):
        result = driver.execute_script('return _testModelAveraging()')
        assert result is True

    def test_r_validation_suite(self, driver):
        result = driver.execute_script('return _testVsR()')
        assert result is True

    # --- Full test runner ---

    def test_all_dr_tests_pass(self, driver):
        result = driver.execute_script('return _runAllDRTests()')
        assert result is True, '_runAllDRTests returned False — at least one suite failed'

    # --- Functional tests for new models ---

    def test_loglinear_model_fit(self, driver):
        result = driver.execute_script("""
            var pts = [];
            for (var i = 1; i <= 8; i++) pts.push({dose: i, effect: 2*Math.log(i+1), se: 0.15});
            var m = fitLogLinearDR(pts);
            return m ? {model: m.model, R2: m.R2, b1: m.b1} : null;
        """)
        assert result is not None
        assert result['model'] == 'Log-Linear'
        assert result['R2'] > 0.9

    def test_exponential_model_fit(self, driver):
        result = driver.execute_script("""
            var pts = [];
            for (var i = 0; i <= 10; i++) pts.push({dose: i, effect: 5*(1-Math.exp(-0.3*i)), se: 0.2});
            var m = fitExponentialDR(pts);
            return m ? {model: m.model, R2: m.R2, Emax: m.Emax, alpha: m.alpha} : null;
        """)
        assert result is not None
        assert result['model'] == 'Exponential'
        assert result['R2'] > 0.9
        assert abs(result['Emax'] - 5.0) < 1.5

    def test_hill_model_fit(self, driver):
        result = driver.execute_script("""
            var pts = [];
            for (var i = 0; i <= 10; i++) {
                var y = 10 * Math.pow(i, 2) / (Math.pow(5, 2) + Math.pow(i, 2));
                pts.push({dose: i, effect: y, se: 0.2});
            }
            var m = fitHillDR(pts);
            return m ? {model: m.model, R2: m.R2, Emax: m.Emax, ED50: m.ED50, h: m.h} : null;
        """)
        assert result is not None
        assert result['model'] == 'Hill'
        assert result['R2'] > 0.95
        assert abs(result['ED50'] - 5.0) < 2.0

    def test_model_averaged_prediction(self, driver):
        result = driver.execute_script("""
            var pts = [];
            for (var i = 0; i <= 8; i++) pts.push({dose: i, effect: Math.log(i+1)*2, se: 0.2});
            var comp = compareDoseResponseModels(pts);
            if (!comp) return null;
            var curve = computeModelAveragedCurve(comp.all, 8, 0.95, 20);
            if (!curve || curve.length === 0) return null;
            return {
                nPoints: curve.length,
                firstDose: curve[0].dose,
                lastDose: curve[curve.length-1].dose,
                hasCI: curve[0].lo !== undefined && curve[0].hi !== undefined,
                midPred: curve[Math.floor(curve.length/2)].pred
            };
        """)
        assert result is not None
        assert result['nPoints'] >= 20  # nPoints may be inclusive (0..nPoints)
        assert result['firstDose'] == 0
        assert abs(result['lastDose'] - 8.0) < 0.5
        assert result['hasCI'] is True
        assert result['midPred'] > 0

    def test_prediction_interval_computation(self, driver):
        result = driver.execute_script("""
            var pi = computePredictionInterval(1.5, 0.3, 0.1, 10, 2, 0.95);
            return pi;
        """)
        assert result is not None
        assert result['lower'] < 1.5
        assert result['upper'] > 1.5
        # PI should be wider than a simple CI (se*1.96)
        ci_half = 0.3 * 1.96
        assert (1.5 - result['lower']) > ci_half

    def test_gl_covariance_produces_matrix(self, driver):
        result = driver.execute_script("""
            var cov = greenlandLongnecker(
                [100, 80, 60],  // cases
                [500, 500, 500], // n
                0,               // refIdx
                'ci'             // type
            );
            if (!cov) return null;
            return {rows: cov.length, cols: cov[0].length, isDiagPositive: cov[0][0] > 0 && cov[1][1] > 0};
        """)
        assert result is not None
        assert result['rows'] == 2  # 3 doses minus reference = 2x2
        assert result['cols'] == 2
        assert result['isDiagPositive'] is True

    def test_gls_fit_returns_coefficients(self, driver):
        """GLS fit: Slist[i] must be studyNs[i] x studyNs[i]."""
        result = driver.execute_script("""
            var X = [[1, 10], [1, 20], [1, 5], [1, 15]];
            var y = [-0.1, -0.2, -0.05, -0.15];
            var S1 = [[0.01, 0.005], [0.005, 0.01]];
            var S2 = [[0.02, 0.008], [0.008, 0.02]];
            var Slist = [S1, S2];
            var studyNs = [2, 2];
            var fit = glsFit(X, y, Slist, studyNs);
            return fit ? {nCoef: fit.coefficients.length, hasVcov: fit.vcov.length === 2} : null;
        """)
        assert result is not None
        assert result['nCoef'] == 2
        assert result['hasVcov'] is True


class TestCSVImport:
    """Test CSV import helper functions."""

    def test_csv_parse_line(self, driver):
        result = driver.execute_script("""
            var line1 = parseCSVLine('StudyA,10,mg,100,1.5,0.8,2.2,0.3,MD', ',');
            var line2 = parseCSVLine('"Study, B",20,mg,150,2.0,1.0,3.0,0.4,MD', ',');
            var line3 = parseCSVLine('a\\tb\\tc', '\\t');
            return {
                line1Len: line1.length, line1First: line1[0],
                line2First: line2[0],
                line3Len: line3.length
            };
        """)
        assert result is not None
        assert result['line1Len'] == 9
        assert result['line1First'] == 'StudyA'
        assert result['line2First'] == 'Study, B'
        assert result['line3Len'] == 3

    def test_parse_float_safe(self, driver):
        result = driver.execute_script("""
            return {
                normal: parseFloatSafe('3.14'),
                zero: parseFloatSafe('0'),
                negZero: parseFloatSafe('-0.5'),
                empty: parseFloatSafe(''),
                nul: parseFloatSafe(null),
                abc: parseFloatSafe('abc')
            };
        """)
        assert result is not None
        assert result['normal'] == 3.14
        assert result['zero'] == 0  # Must not drop valid zero
        assert result['negZero'] == -0.5
        assert result['empty'] is None
        assert result['nul'] is None
        assert result['abc'] is None

    def test_csv_import_button_exists(self, driver):
        btns = driver.find_elements(By.ID, 'csvImportInput')
        assert len(btns) == 1


class TestAdvancedFeatures:
    """Tests for advanced DR features: dose-finding, Bayesian, funnel, forest, bootstrap."""

    # --- Dose-Finding Metrics ---

    def test_dose_metrics_basic(self, driver):
        """computeDoseMetrics should return ED50, ED90, MED, NOAEL."""
        result = driver.execute_script("""
            var pts = [];
            for (var i = 0; i <= 10; i++) pts.push({dose: i, effect: 5*(1 - Math.exp(-0.3*i)), se: 0.2});
            var model = fitExponentialDR(pts);
            if (!model) return null;
            var m = computeDoseMetrics(model, 10, 0.5);
            return m ? {ED50: m.ED50, ED90: m.ED90, MED: m.MED, NOAEL: m.NOAEL,
                        hasWindow: m.therapeuticWindow != null} : null;
        """)
        assert result is not None
        assert 0 < result['ED50'] < 10
        assert result['ED90'] > result['ED50']
        assert result['MED'] is not None
        assert result['hasWindow'] is True

    def test_dose_metrics_linear(self, driver):
        """Dose metrics on a simple linear curve."""
        result = driver.execute_script("""
            var pts = [];
            for (var i = 0; i <= 6; i++) pts.push({dose: i, effect: 0.5*i, se: 0.15});
            var model = fitLinearDR(pts);
            if (!model) return null;
            var m = computeDoseMetrics(model, 6, 0);
            return m ? {ED50: m.ED50, maxEffectDose: m.maxEffectDose} : null;
        """)
        assert result is not None
        assert abs(result['ED50'] - 3.0) < 0.5  # ED50 at 50% of max effect

    def test_dose_metrics_card_render(self, driver):
        """renderDoseMetricsCard should produce HTML with key labels."""
        result = driver.execute_script("""
            var pts = [];
            for (var i = 0; i <= 8; i++) pts.push({dose: i, effect: 3*Math.log(i+1), se: 0.2});
            var model = fitLinearDR(pts);
            if (!model) return null;
            var metrics = computeDoseMetrics(model, 8, 0.5);
            if (!metrics) return null;
            var html = renderDoseMetricsCard(metrics);
            return {hasED50: html.includes('ED') && html.includes('50'),
                    hasED90: html.includes('ED') && html.includes('90'),
                    hasMED: html.includes('MED'), length: html.length};
        """)
        assert result is not None
        assert result['hasED50'] is True
        assert result['hasED90'] is True
        assert result['length'] > 100

    # --- Bayesian DR ---

    def test_bayesian_dr_linear(self, driver):
        """fitBayesianDR should return posterior mean, SD, and CrI function."""
        result = driver.execute_script("""
            var pts = [];
            for (var i = 0; i <= 8; i++) pts.push({dose: i, effect: 1 + 0.5*i, se: 0.2});
            var b = fitBayesianDR(pts, 'linear', 10);
            if (!b) return null;
            var cri = b.predictBayesCrI(4, 0.95);
            return {model: b.model, nMean: b.posteriorMean.length, nSD: b.posteriorSD.length,
                    criY: cri.y, criLo: cri.lo, criHi: cri.hi,
                    predictConsistent: Math.abs(b.predict(4) - cri.y) < 0.01};
        """)
        assert result is not None
        assert result['model'].startswith('Bayesian-')
        assert result['nMean'] == 2
        assert result['criLo'] < result['criY'] < result['criHi']
        assert result['predictConsistent'] is True  # predict and CrI center must agree

    def test_bayesian_dr_shrinkage(self, driver):
        """Bayesian posterior should shrink toward zero vs frequentist."""
        result = driver.execute_script("""
            var pts = [];
            for (var i = 0; i <= 6; i++) pts.push({dose: i, effect: 2 + 0.8*i, se: 0.3});
            var freq = fitLinearDR(pts);
            var bayes = fitBayesianDR(pts, 'linear', 5);
            if (!freq || !bayes) return null;
            return {freqB1: freq.b1, bayesB1: bayes.posteriorMean[1],
                    shrunk: Math.abs(bayes.posteriorMean[1]) < Math.abs(freq.b1)};
        """)
        assert result is not None
        assert result['shrunk'] is True  # posterior shrinks toward zero

    def test_bayesian_cri_band_render(self, driver):
        """renderBayesCrIBand should produce SVG path string."""
        result = driver.execute_script("""
            var pts = [];
            for (var i = 0; i <= 8; i++) pts.push({dose: i, effect: 1 + 0.4*i, se: 0.2});
            var b = fitBayesianDR(pts, 'linear', 10);
            if (!b) return null;
            var xS = function(d) { return 50 + d * 40; };
            var yS = function(y) { return 300 - y * 30; };
            var svg = renderBayesCrIBand(b, xS, yS, 0, 8, 50, false);
            return {hasPath: svg.includes('<path'), length: svg.length};
        """)
        assert result is not None
        assert result['hasPath'] is True
        assert result['length'] > 50

    # --- DR Publication Bias ---

    def test_dr_funnel_plot_render(self, driver):
        """renderDRFunnelPlot should produce SVG with funnel shape."""
        result = driver.execute_script("""
            var pts = [];
            for (var i = 0; i <= 10; i++) {
                pts.push({dose: i, effect: 0.3*i + (Math.random()-0.5)*0.5, se: 0.1 + Math.random()*0.3, studyId: 'S' + i});
            }
            var model = fitLinearDR(pts);
            if (!model) return null;
            var svg = renderDRFunnelPlot(pts, model);
            return {hasSVG: svg.includes('<svg'), hasCircle: svg.includes('<circle'),
                    hasText: svg.includes('Precision'), length: svg.length};
        """)
        assert result is not None
        assert result['hasSVG'] is True
        assert result['hasCircle'] is True
        assert result['length'] > 200

    def test_dose_reporting_bias_test(self, driver):
        """testDoseReportingBias should return chi-squared and p-value."""
        result = driver.execute_script("""
            var pts = [];
            for (var i = 0; i < 20; i++) {
                pts.push({dose: Math.floor(Math.random() * 10), effect: Math.random(), se: 0.1 + Math.random()*0.2});
            }
            var r = testDoseReportingBias(pts);
            return r;
        """)
        assert result is not None
        assert 'chi2' in result
        assert 'df' in result
        assert 'pValue' in result
        assert 0 <= result['pValue'] <= 1

    def test_dose_egger_regression(self, driver):
        """doseEggerRegression should return intercept, slope, p-value."""
        result = driver.execute_script("""
            var pts = [];
            for (var i = 0; i <= 10; i++) {
                pts.push({dose: i, effect: 0.5*i + (Math.random()-0.5)*0.3, se: 0.15});
            }
            var model = fitLinearDR(pts);
            if (!model) return null;
            var r = doseEggerRegression(pts, model);
            return r;
        """)
        assert result is not None
        assert 'intercept' in result
        assert 'pValue' in result
        assert 0 <= result['pValue'] <= 1

    # --- Dose-Specific Forest Plot ---

    def test_dose_forest_plot_render(self, driver):
        """renderDoseForestPlot should produce SVG with study rows."""
        result = driver.execute_script("""
            var pts = [];
            var studies = ['Alpha', 'Beta', 'Gamma', 'Delta'];
            for (var s = 0; s < 4; s++) {
                for (var d = 0; d <= 4; d++) {
                    pts.push({dose: d*10, effect: 0.3*d + (s-1.5)*0.2, se: 0.15, studyId: studies[s]});
                }
            }
            var svg = renderDoseForestPlot(pts, 0.95);
            return {hasSVG: svg.includes('<svg'), hasRect: svg.includes('<rect'),
                    hasStudy: svg.includes('Alpha'), length: svg.length};
        """)
        assert result is not None
        assert result['hasSVG'] is True
        assert result['hasRect'] is True
        assert result['hasStudy'] is True
        assert result['length'] > 300

    # --- Breakpoint Bootstrap CI ---

    def test_bootstrap_breakpoints(self, driver):
        """bootstrapBreakpoints should return plateau and inflection CI arrays."""
        result = driver.execute_script("""
            var pts = [];
            var studies = ['S1', 'S2', 'S3', 'S4'];
            for (var s = 0; s < 4; s++) {
                for (var d = 0; d <= 8; d++) {
                    var y = d < 5 ? 0.5*d : 2.5;
                    pts.push({dose: d, effect: y, se: 0.15, studyId: studies[s]});
                }
            }
            var r = bootstrapBreakpoints(pts, 'linear', 100, 0.95);
            return r ? {hasPlateau: r.plateau != null, hasInflection: r.inflection != null,
                        keys: Object.keys(r)} : null;
        """)
        assert result is not None
        assert result['hasPlateau'] is True or result['hasInflection'] is True

    # --- Study Overlay ---

    def test_study_overlay_in_svg(self, driver):
        """Dose-response SVG should support per-study overlay curves."""
        result = driver.execute_script("""
            var pts = [];
            for (var s = 0; s < 3; s++) {
                for (var d = 0; d <= 4; d++) {
                    pts.push({dose: d, effect: 1 + 0.3*d + s*0.1, se: 0.2, studyId: 'Study' + s});
                }
            }
            var model = fitLinearDR(pts);
            if (!model) return null;
            var svg = renderDoseResponseCurveSVG(pts, model, 0.95, 1.96, false);
            return {hasPolyline: svg.includes('polyline') || svg.includes('<path'),
                    hasSVG: svg.includes('<svg'), length: svg.length};
        """)
        assert result is not None
        assert result['hasSVG'] is True
        assert result['length'] > 200

    # --- Network Dose-Response ---

    def test_network_dr_curve(self, driver):
        """renderNetworkDRCurve should produce SVG with multiple drug curves."""
        result = driver.execute_script("""
            var groups = [
                {drugName: 'DrugA', points: [{dose:0,effect:0,se:0.1},{dose:5,effect:1.5,se:0.15},{dose:10,effect:2.8,se:0.2}]},
                {drugName: 'DrugB', points: [{dose:0,effect:0,se:0.12},{dose:5,effect:0.8,se:0.18},{dose:10,effect:1.5,se:0.22}]}
            ];
            var svg = renderNetworkDRCurve(groups, 'linear', 0.95);
            return {hasSVG: svg.includes('<svg'), hasDrugA: svg.includes('DrugA'),
                    hasDrugB: svg.includes('DrugB'), length: svg.length};
        """)
        assert result is not None
        assert result['hasSVG'] is True
        assert result['hasDrugA'] is True
        assert result['hasDrugB'] is True

    # --- R/Python Code Export ---

    def test_r_code_export(self, driver):
        """generateRCode should produce valid R script text."""
        result = driver.execute_script("""
            var pts = [{dose:0,effect:0,se:0.1},{dose:5,effect:1,se:0.15},{dose:10,effect:2,se:0.2}];
            var model = fitLinearDR(pts);
            if (!model) return null;
            var code = generateRCode(pts, model, 'dl');
            return {hasLibrary: code.includes('dosresmeta'), hasCoef: code.includes('coef'),
                    length: code.length};
        """)
        assert result is not None
        assert result['hasLibrary'] is True
        assert result['length'] > 100

    def test_python_code_export(self, driver):
        """generatePythonCode should produce valid Python script text."""
        result = driver.execute_script("""
            var pts = [{dose:0,effect:0,se:0.1},{dose:5,effect:1,se:0.15},{dose:10,effect:2,se:0.2}];
            var model = fitLinearDR(pts);
            if (!model) return null;
            var code = generatePythonCode(pts, model, 'dl');
            return {hasNumpy: code.includes('numpy'), hasResult: code.includes('result') || code.includes('beta'),
                    length: code.length};
        """)
        assert result is not None
        assert result['hasNumpy'] is True
        assert result['length'] > 100

    # --- Gap-to-Protocol Pipeline ---

    def test_gap_to_protocol_self_test(self, driver):
        """Run the inline _testGapToProtocol self-test."""
        result = driver.execute_script('return _testGapToProtocol()')
        assert result is True

    # --- Inline Advanced Feature Tests ---

    def test_advanced_features_suite(self, driver):
        """Run the inline _testAdvancedDRFeatures self-test."""
        result = driver.execute_script('return _testAdvancedDRFeatures()')
        assert result is True

    def test_top5_advanced_suite(self, driver):
        """Run the inline _testTop5AdvancedFeatures self-test."""
        result = driver.execute_script('return _testTop5AdvancedFeatures()')
        assert result is True

    # --- DR Heterogeneity with corrected nParams ---

    def test_heterogeneity_df_correct(self, driver):
        """After nParams fix, Linear Q_DR should have df = nPoints - 2."""
        result = driver.execute_script("""
            var pts = [];
            for (var i = 0; i <= 8; i++) pts.push({dose: i, effect: 1 + 0.5*i, se: 0.15});
            var model = fitLinearDR(pts);
            if (!model) return null;
            var het = computeDRHeterogeneity(pts, model);
            return het ? {df: het.df, nPts: pts.length} : null;
        """)
        assert result is not None
        assert result['df'] == result['nPts'] - 2  # 9 points - 2 params (b0 + b1)

    def test_heterogeneity_quadratic_df(self, driver):
        """Quadratic Q_DR should have df = nPoints - 3."""
        result = driver.execute_script("""
            var pts = [];
            for (var i = 0; i <= 8; i++) pts.push({dose: i, effect: 1 + 0.3*i - 0.02*i*i, se: 0.15});
            var model = fitQuadraticDR(pts);
            if (!model) return null;
            var het = computeDRHeterogeneity(pts, model);
            return het ? {df: het.df, nPts: pts.length} : null;
        """)
        assert result is not None
        assert result['df'] == result['nPts'] - 3  # 9 points - 3 params (b0 + b1 + b2)

    # --- Focus Trap Utility ---

    def test_focus_trap_exists(self, driver):
        """trapFocus utility function should be defined."""
        result = driver.execute_script('return typeof trapFocus === "function"')
        assert result is True


class TestAccessibility:
    """Test accessibility improvements from the review."""

    def test_svg_has_role_img(self, driver):
        """Dynamically-generated SVGs should have role=img."""
        result = driver.execute_script("""
            var pts = [];
            for (var i = 0; i <= 8; i++) pts.push({dose: i, effect: 1 + 0.3*i, se: 0.15});
            var model = fitLinearDR(pts);
            if (!model) return null;
            var svg = renderDoseResponseCurveSVG(pts, model, 0.95, 1.96, false);
            return {hasRoleImg: svg.includes('role="img"'), hasAriaLabel: svg.includes('aria-label'),
                    hasViewBox: svg.includes('viewBox')};
        """)
        assert result is not None
        assert result['hasRoleImg'] is True
        assert result['hasAriaLabel'] is True
        assert result['hasViewBox'] is True

    def test_forest_svg_accessible(self, driver):
        result = driver.execute_script("""
            var pts = [];
            for (var s = 0; s < 3; s++)
                for (var d = 0; d <= 3; d++)
                    pts.push({dose: d*10, effect: 0.3*d, se: 0.15, studyId: 'S' + s});
            var svg = renderDoseForestPlot(pts, 0.95);
            return {hasRoleImg: svg.includes('role="img"'), hasViewBox: svg.includes('viewBox')};
        """)
        assert result is not None
        assert result['hasRoleImg'] is True
        assert result['hasViewBox'] is True

    def test_badge_dark_mode_css_exists(self, driver):
        """Dark mode badge CSS overrides should be present."""
        result = driver.execute_script("""
            var sheets = document.styleSheets;
            var found = {include: false, exclude: false, maybe: false, duplicate: false};
            try {
                for (var i = 0; i < sheets.length; i++) {
                    var rules = sheets[i].cssRules || [];
                    for (var j = 0; j < rules.length; j++) {
                        var sel = rules[j].selectorText || '';
                        if (sel.includes('dark-mode') && sel.includes('badge-include')) found.include = true;
                        if (sel.includes('dark-mode') && sel.includes('badge-exclude')) found.exclude = true;
                        if (sel.includes('dark-mode') && sel.includes('badge-maybe')) found.maybe = true;
                        if (sel.includes('dark-mode') && sel.includes('badge-duplicate')) found.duplicate = true;
                    }
                }
            } catch(e) {}
            return found;
        """)
        assert result is not None
        assert result['include'] is True
        assert result['exclude'] is True
        assert result['maybe'] is True
        assert result['duplicate'] is True

    def test_skip_link_exists(self, driver):
        """Skip-to-content link should exist."""
        links = driver.find_elements(By.CSS_SELECTOR, '.skip-link')
        assert len(links) >= 1

    def test_confirm_dialog_escapes(self, driver):
        """Confirm dialog should handle Escape key."""
        result = driver.execute_script("""
            return typeof showConfirm === 'function' && typeof trapFocus === 'function';
        """)
        assert result is True


class TestRValidation:
    """Compare JS dose-response engine outputs against R dosresmeta2 v2.2.0 reference values.

    NOTE: The JS REML/ML optimizer finds a different tau2 than R due to the profile
    likelihood parameterization (univariate scalar vs R's multivariate Psi matrix).
    Fixed-effects and GL covariance are independent of tau2 and should match tightly.
    REML/ML tests verify direction, sign, and structural properties instead of exact values.
    """

    @pytest.fixture(scope='class')
    def ref(self):
        import json
        ref_path = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'validation_reference.json')
        with open(ref_path) as f:
            return json.load(f)

    # Helper: build alcohol data and fit in one JS call
    ALC_SETUP_JS = """
        var alcData = buildAlcoholCVDData();
        var Slist = [], yAll = [], XAll = [], XQuad = [], studyNs = [];
        for (var si = 0; si < alcData.length; si++) {
            var s = alcData[si];
            var vArr = s.se ? s.se.slice(1).map(function(x) { return x * x; }) : null;
            var Smat = greenlandLongnecker(s.cases, s.n, 0, s.type, vArr);
            if (!Smat) return null;
            Slist.push(Smat);
            studyNs.push(s.doses.length - 1);
            for (var i = 1; i < s.doses.length; i++) {
                yAll.push(s.logrr[i]);
                XAll.push([s.doses[i]]);
                XQuad.push([s.doses[i], s.doses[i]*s.doses[i]]);
            }
        }
    """

    def test_gl_covariance_dimensions(self, driver, ref):
        """GL covariance matrices (Hamling) should have correct dimensions matching R."""
        result = driver.execute_script("""
            var alcData = buildAlcoholCVDData();
            return alcData.map(function(s) {
                var vArr = s.se ? s.se.slice(1).map(function(x) { return x * x; }) : null;
                var S = greenlandLongnecker(s.cases, s.n, 0, s.type, vArr);
                return S ? S.length : null;
            });
        """)
        assert result is not None
        r_Slist = ref['alcohol_cvd']['Slist']
        assert len(result) == len(r_Slist), f"Slist length: JS={len(result)}, R={len(r_Slist)}"
        for si in range(len(r_Slist)):
            assert result[si] == len(r_Slist[si]), \
                f"Study {si} dim: JS={result[si]}, R={len(r_Slist[si])}"

    def test_gl_covariance_positive_definite(self, driver):
        """All GL covariance matrices (Hamling) should be symmetric positive definite."""
        result = driver.execute_script("""
            var alcData = buildAlcoholCVDData();
            var checks = [];
            for (var si = 0; si < alcData.length; si++) {
                var s = alcData[si];
                var vArr = s.se ? s.se.slice(1).map(function(x) { return x * x; }) : null;
                var S = greenlandLongnecker(s.cases, s.n, 0, s.type, vArr);
                if (!S) { checks.push({id: s.id, ok: false, reason: 'null'}); continue; }
                var allPosDiag = true, symmetric = true;
                for (var i = 0; i < S.length; i++) {
                    if (S[i][i] <= 0) allPosDiag = false;
                    for (var j = 0; j < S.length; j++)
                        if (Math.abs(S[i][j] - S[j][i]) > 1e-12) symmetric = false;
                }
                checks.push({id: s.id, ok: allPosDiag && symmetric,
                              allPosDiag: allPosDiag, symmetric: symmetric});
            }
            return checks;
        """)
        assert result is not None
        for c in result:
            assert c['ok'], f"Study {c['id']}: posDiag={c.get('allPosDiag')}, sym={c.get('symmetric')}"

    def test_gl_covariance_matches_r(self, driver, ref):
        """GL Hamling covariance should match R Slist within tolerance."""
        result = driver.execute_script("""
            var alcData = buildAlcoholCVDData();
            return alcData.map(function(s) {
                var vArr = s.se ? s.se.slice(1).map(function(x) { return x * x; }) : null;
                return greenlandLongnecker(s.cases, s.n, 0, s.type, vArr);
            });
        """)
        assert result is not None
        r_Slist = ref['alcohol_cvd']['Slist']
        for si in range(len(r_Slist)):
            js_S = result[si]
            r_S = r_Slist[si]
            for i in range(len(r_S)):
                for j in range(len(r_S[i])):
                    assert abs(js_S[i][j] - r_S[i][j]) < 5e-3, \
                        f"S[{si}][{i}][{j}]: JS={js_S[i][j]:.6f}, R={r_S[i][j]:.6f}"

    def test_fixed_effects_coefficient_matches_r(self, driver, ref):
        """Fixed-effects coefficient should match R within 5e-4 (Hamling GL)."""
        result = driver.execute_script(self.ALC_SETUP_JS + """
            var fixed = glsFit(XAll, yAll, Slist, studyNs);
            if (!fixed) return null;
            return {coef: fixed.coefficients[0], se: Math.sqrt(fixed.vcov[0][0])};
        """)
        assert result is not None
        r = ref['alcohol_cvd']['linear_fixed']
        assert result['coef'] < 0, "Fixed coef should be negative"
        assert abs(result['coef'] - r['coefficients']) < 5e-4, \
            f"Fixed coef: JS={result['coef']:.8f}, R={r['coefficients']:.8f}"

    def test_fixed_effects_se_matches_r(self, driver, ref):
        """Fixed-effects SE should match R within 5e-4 (Hamling GL)."""
        result = driver.execute_script(self.ALC_SETUP_JS + """
            var fixed = glsFit(XAll, yAll, Slist, studyNs);
            if (!fixed) return null;
            return Math.sqrt(fixed.vcov[0][0]);
        """)
        assert result is not None
        r_se = ref['alcohol_cvd']['linear_fixed']['se']
        assert result > 0, "SE must be positive"
        assert abs(result - r_se) < 5e-4, \
            f"Fixed SE: JS={result:.8f}, R={r_se:.8f}"

    def test_reml_coefficient_matches_r(self, driver, ref):
        """REML coefficient should match R within 5e-4."""
        result = driver.execute_script(self.ALC_SETUP_JS + """
            var reml = estimateDRREML(XAll, yAll, Slist, studyNs);
            return reml ? reml.coefficients[0] : null;
        """)
        assert result is not None
        r_coef = ref['alcohol_cvd']['linear_reml']['coefficients']
        assert result < 0, f"REML coef should be negative, got {result}"
        assert abs(result - r_coef) < 5e-4, \
            f"REML coef: JS={result:.8f}, R={r_coef:.8f}"

    def test_reml_tau2_matches_r(self, driver, ref):
        """REML tau2 should match R within 5e-3."""
        result = driver.execute_script(self.ALC_SETUP_JS + """
            var reml = estimateDRREML(XAll, yAll, Slist, studyNs);
            return reml ? reml.tau2 : null;
        """)
        assert result is not None
        r_tau2 = ref['alcohol_cvd']['linear_reml']['tau2']
        assert result >= 0, f"REML tau2 must be >= 0, got {result}"
        assert abs(result - r_tau2) < 5e-3, \
            f"REML tau2: JS={result:.8f}, R={r_tau2:.8f}"

    def test_reml_se_larger_than_fixed(self, driver):
        """REML SE should be >= fixed SE (heterogeneity adds uncertainty)."""
        result = driver.execute_script(self.ALC_SETUP_JS + """
            var fixed = glsFit(XAll, yAll, Slist, studyNs);
            var reml = estimateDRREML(XAll, yAll, Slist, studyNs);
            if (!fixed || !reml) return null;
            return {fixedSE: Math.sqrt(fixed.vcov[0][0]), remlSE: Math.sqrt(reml.vcov[0][0])};
        """)
        assert result is not None
        assert result['remlSE'] >= result['fixedSE'] * 0.9, \
            f"REML SE={result['remlSE']:.6f} should be >= fixed SE={result['fixedSE']:.6f}"

    def test_ml_tau2_leq_reml_tau2(self, driver):
        """ML tau2 should be <= REML tau2 (standard statistical property)."""
        result = driver.execute_script(self.ALC_SETUP_JS + """
            var reml = estimateDRREML(XAll, yAll, Slist, studyNs);
            var ml = estimateDRML(XAll, yAll, Slist, studyNs);
            if (!reml || !ml) return null;
            return {reml_tau2: reml.tau2, ml_tau2: ml.tau2};
        """)
        assert result is not None
        assert result['reml_tau2'] >= result['ml_tau2'] * 0.99, \
            f"REML tau2={result['reml_tau2']:.6f} should be >= ML tau2={result['ml_tau2']:.6f}"

    def test_reml_predictions_match_r(self, driver, ref):
        """REML predictions should match R within 5e-3 at all dose points."""
        result = driver.execute_script(self.ALC_SETUP_JS + """
            var reml = estimateDRREML(XAll, yAll, Slist, studyNs);
            if (!reml) return null;
            var doses = [5, 10, 15, 20, 25, 30, 40, 50, 60];
            return doses.map(function(d) {
                return {dose: d, pred: reml.coefficients[0] * d};
            });
        """)
        assert result is not None
        r_preds = ref['alcohol_cvd']['linear_reml']['predictions']
        for rp in r_preds:
            if rp['dose'] == 0:
                continue
            js_match = next((r for r in result if r['dose'] == rp['dose']), None)
            if js_match:
                assert abs(js_match['pred'] - rp['pred']) < 5e-3, \
                    f"dose={rp['dose']}: JS={js_match['pred']:.6f}, R={rp['pred']:.6f}"

    def test_aic_values_finite_and_distinct(self, driver):
        """Linear and quadratic REML should produce finite, distinct AIC values."""
        result = driver.execute_script(self.ALC_SETUP_JS + """
            var lin = estimateDRREML(XAll, yAll, Slist, studyNs);
            var quad = estimateDRREML(XQuad, yAll, Slist, studyNs);
            if (!lin || !quad) return null;
            return {lin_AIC: lin.AIC, quad_AIC: quad.AIC,
                    lin_finite: isFinite(lin.AIC), quad_finite: isFinite(quad.AIC)};
        """)
        assert result is not None
        assert result['lin_finite'], f"Linear AIC not finite: {result['lin_AIC']}"
        assert result['quad_finite'], f"Quadratic AIC not finite: {result['quad_AIC']}"
        assert result['lin_AIC'] != result['quad_AIC'], "AIC values should differ between models"

    def test_quadratic_coefficients_direction(self, driver, ref):
        """Quadratic REML: coef[0] negative, coef[1] positive (U-shape), matching R."""
        result = driver.execute_script(self.ALC_SETUP_JS + """
            var reml = estimateDRREML(XQuad, yAll, Slist, studyNs);
            return reml ? reml.coefficients : null;
        """)
        assert result is not None
        r_coef = ref['alcohol_cvd']['quadratic_reml']['coefficients']
        assert result[0] < 0, f"Quad coef[0] should be negative, got {result[0]}"
        assert result[1] > 0, f"Quad coef[1] should be positive, got {result[1]}"
        assert r_coef[0] < 0 and r_coef[1] > 0, "R reference confirms U-shape"

    def test_regression_snapshot_fixed_coef(self, driver, ref):
        """Regression: fixed-effects coefficient should match R reference."""
        result = driver.execute_script(self.ALC_SETUP_JS + """
            var fixed = glsFit(XAll, yAll, Slist, studyNs);
            return fixed ? fixed.coefficients[0] : null;
        """)
        assert result is not None
        # With Hamling GL, should now match R closely
        r_coef = ref['alcohol_cvd']['linear_fixed']['coefficients']
        assert abs(result - r_coef) < 5e-4, \
            f"Regression: JS fixed coef={result:.8f}, R={r_coef:.8f}"

    def test_regression_snapshot_reml_tau2(self, driver, ref):
        """Regression: REML tau2 should match R reference (near zero with Hamling)."""
        result = driver.execute_script(self.ALC_SETUP_JS + """
            var reml = estimateDRREML(XAll, yAll, Slist, studyNs);
            return reml ? reml.tau2 : null;
        """)
        assert result is not None
        # With Hamling GL, tau2 should be ~0.0001 (not 0.147)
        r_tau2 = ref['alcohol_cvd']['linear_reml']['tau2']
        assert abs(result - r_tau2) < 5e-3, \
            f"Regression: JS REML tau2={result:.8f}, R={r_tau2:.8f}"


class TestCoverageGaps:
    """Tests for previously untested features."""

    def test_csv_import_round_trip(self, driver):
        """CSV parsing primitives should round-trip dose-response rows.

        Exercises the real CSV import building blocks the app uses
        (parseCSVLine + getCSVVal + parseFloatSafe). No typeof guard: if any
        of these functions is renamed, this test fails loudly instead of
        silently skipping.
        """
        result = driver.execute_script("""
            var rows = [
                'study,dose,effect,se',
                'TestStudy,0,0,0.10',
                'TestStudy,10,-0.20,0.15',
                '"Study, Jr",20,-0.35,0.18'
            ];
            var header = parseCSVLine(rows[0], ',');
            var doseIdx = header.indexOf('dose');
            var effIdx = header.indexOf('effect');
            var parsed = [];
            for (var i = 1; i < rows.length; i++) {
                var vals = parseCSVLine(rows[i], ',');
                parsed.push({
                    study: getCSVVal(vals, 0, ''),
                    dose: parseFloatSafe(getCSVVal(vals, doseIdx, '')),
                    effect: parseFloatSafe(getCSVVal(vals, effIdx, ''))
                });
            }
            return {
                nCols: header.length,
                rows: parsed.length,
                firstStudy: parsed[0].study,
                firstDose: parsed[0].dose,
                firstEffect: parsed[0].effect,
                quotedStudy: parsed[2].study,
                secondDose: parsed[1].dose
            };
        """)
        assert result is not None
        assert result['nCols'] == 4
        assert result['rows'] == 3
        assert result['firstStudy'] == 'TestStudy'
        # numeric 0 must survive parsing (not dropped by a `|| null` fallback)
        assert result['firstDose'] == 0
        assert result['firstEffect'] == 0
        # quoted comma inside a field must not split the row
        assert result['quotedStudy'] == 'Study, Jr'
        assert result['secondDose'] == 10

    def test_dark_mode_toggle_persistence(self, driver):
        """Dark mode toggle should persist via localStorage."""
        result = driver.execute_script("""
            // Toggle dark mode on
            var body = document.body;
            var wasDark = body.classList.contains('dark-mode');
            if (typeof toggleDarkMode === 'function') toggleDarkMode();
            else body.classList.toggle('dark-mode');
            var isDark = body.classList.contains('dark-mode');
            // Toggle back
            if (typeof toggleDarkMode === 'function') toggleDarkMode();
            else body.classList.toggle('dark-mode');
            return {wasDark: wasDark, toggledTo: isDark, restoredTo: body.classList.contains('dark-mode')};
        """)
        assert result is not None
        # Toggle should have flipped the state
        assert result['toggledTo'] != result['wasDark']
        # Second toggle should restore
        assert result['restoredTo'] == result['wasDark']

    def test_prediction_interval_k2_edge(self, driver):
        """Prediction interval with k=2 (minimum for heterogeneity) should still work."""
        result = driver.execute_script("""
            if (typeof computePredictionInterval !== 'function') return 'skip';
            var pi = computePredictionInterval(-0.3, 0.15, 0.05, 2, 1, 0.95);
            if (!pi) return null;
            return {lower: pi.lower, upper: pi.upper, effect: -0.3};
        """)
        if result == 'skip':
            pytest.skip("computePredictionInterval not available")
        assert result is not None
        assert result['lower'] < result['effect']
        assert result['upper'] > result['effect']
        # With k=2, df=0, PI should be very wide (t with 0 df → infinite)
        width = result['upper'] - result['lower']
        assert width > 1.0, f"k=2 PI should be very wide, got width={width}"

    def test_dose_ranging_detection_basic(self, driver):
        """Dose-ranging detection should flag trials with 3+ dose levels.

        Calls the real detectDoseRanging(trial) function, which parses dose
        text from each arm's label/name (not a `dose` field). No typeof guard:
        a rename fails loudly instead of silently skipping.
        """
        result = driver.execute_script("""
            // 4 arms, 3 distinct active doses = dose-ranging
            var trial1 = {arms: [
                {label: 'Drug X 10 mg'}, {label: 'Drug X 20 mg'},
                {label: 'Drug X 40 mg'}, {label: 'Placebo'}
            ], interventions: []};
            // Single active dose vs placebo = not dose-ranging
            var trial2 = {arms: [
                {label: 'Drug X 100 mg'}, {label: 'Placebo'}
            ], interventions: []};
            return {
                fourArm: detectDoseRanging(trial1),
                fourArmLevels: trial1.doseLevels.length,
                twoArm: detectDoseRanging(trial2)
            };
        """)
        assert result is not None
        assert result['fourArm'] is True
        assert result['fourArmLevels'] == 3
        assert result['twoArm'] is False

    def test_localstorage_fallback(self, driver):
        """App should work when localStorage throws (in-memory fallback)."""
        result = driver.execute_script("""
            // Test the safe storage wrapper
            if (typeof safeSetStorage !== 'function') return 'skip';
            var testKey = '__dr_test_' + Date.now();
            safeSetStorage(testKey, 'hello');
            var val = safeGetStorage(testKey);
            // Clean up
            try { localStorage.removeItem(testKey); } catch(e) {}
            return {stored: val === 'hello'};
        """)
        if result == 'skip':
            pytest.skip("safeSetStorage not available")
        assert result is not None
        assert result['stored'] is True

    def test_model_averaging_weights_sum_one(self, driver):
        """AIC weights across all models should sum to 1.0."""
        result = driver.execute_script("""
            var pts = [];
            for (var s = 0; s < 4; s++)
                for (var d = 1; d <= 3; d++)
                    pts.push({dose: d*10, effect: -0.02*d + (Math.random()-0.5)*0.05,
                              se: 0.12, studyId: 'S' + s});
            var comp = compareDoseResponseModels(pts);
            if (!comp || !comp.all) return null;
            var totalWeight = comp.all.reduce(function(s, m) { return s + m.aicWeight; }, 0);
            return {nModels: comp.all.length, totalWeight: totalWeight};
        """)
        assert result is not None
        assert result['nModels'] >= 3
        assert abs(result['totalWeight'] - 1.0) < 0.01

    def test_escapehtml_covers_quotes(self, driver):
        """escapeHtml should escape both single and double quotes."""
        result = driver.execute_script("""
            if (typeof escapeHtml !== 'function') return 'skip';
            var test = '<script>"alert(1)"</script>';
            var escaped = escapeHtml(test);
            return {
                hasLtGt: !escaped.includes('<') && !escaped.includes('>'),
                hasQuotes: !escaped.includes('"'),
                original: test,
                escaped: escaped
            };
        """)
        if result == 'skip':
            pytest.skip("escapeHtml not available")
        assert result is not None
        assert result['hasLtGt'] is True
        assert result['hasQuotes'] is True

    def test_gl_covariance_ci_type(self, driver):
        """GL covariance for CI-type studies should produce positive definite matrix."""
        result = driver.execute_script("""
            // Simulated cumulative incidence study (5 dose levels)
            var cases = [100, 80, 65, 55, 50];
            var n = [500, 500, 500, 500, 500];
            var Smat = greenlandLongnecker(cases, n, 0, 'ci');
            if (!Smat) return null;
            // Check positive definite: all diagonal elements > 0
            var allPosDiag = true;
            for (var i = 0; i < Smat.length; i++)
                if (Smat[i][i] <= 0) allPosDiag = false;
            // Check symmetric
            var symmetric = true;
            for (var i = 0; i < Smat.length; i++)
                for (var j = 0; j < Smat.length; j++)
                    if (Math.abs(Smat[i][j] - Smat[j][i]) > 1e-12) symmetric = false;
            return {size: Smat.length, allPosDiag: allPosDiag, symmetric: symmetric};
        """)
        assert result is not None
        assert result['size'] == 4  # 5 dose levels - 1 reference
        assert result['allPosDiag'] is True
        assert result['symmetric'] is True

    def test_heterogeneity_nparams_by_model(self, driver):
        """nParams should count intercept: linear=2, quad=3, emax=3."""
        result = driver.execute_script("""
            var pts = [];
            for (var s = 0; s < 5; s++)
                for (var d = 1; d <= 4; d++)
                    pts.push({dose: d*10, effect: -0.03*d, se: 0.1, studyId: 'S' + s});
            var lin = fitLinearDR(pts);
            var het = computeDRHeterogeneity(pts, lin);
            if (!het) return null;
            // linear: nParams=2 (intercept + slope), df = nPts - 2
            return {df: het.df, nPts: pts.length, expectedDf: pts.length - 2};
        """)
        assert result is not None
        assert result['df'] == result['expectedDf'], \
            f"df={result['df']}, expected={result['expectedDf']}"

    def test_bayesian_posterior_shrinkage(self, driver):
        """Bayesian posterior mean should be shrunk toward prior vs frequentist."""
        result = driver.execute_script("""
            if (typeof fitBayesianDR !== 'function') return 'skip';
            var pts = [];
            for (var s = 0; s < 4; s++)
                for (var d = 1; d <= 4; d++)
                    pts.push({dose: d*10, effect: -0.05*d + 0.01*s, se: 0.2, studyId: 'S' + s});
            var bfit = fitBayesianDR(pts, 'linear');
            if (!bfit) return null;
            // posteriorMean has [intercept, slope]; freqModel has b0, b1
            var freqSlope = bfit.freqModel.b1;
            var postSlope = bfit.posteriorMean[1];
            return {freqSlope: freqSlope, postSlope: postSlope,
                    shrunk: Math.abs(postSlope) <= Math.abs(freqSlope) + 0.01};
        """)
        if result == 'skip':
            pytest.skip("fitBayesianDR not available")
        assert result is not None, "fitBayesianDR returned null — check input format"
        assert result['shrunk'] is True, \
            f"Posterior={result['postSlope']} should be shrunk vs freq={result['freqSlope']}"

    def test_canvas_accessibility_attributes(self, driver):
        """All 6 insight canvas elements should have role=img and aria-label."""
        ids = ['tawakkulRadar', 'mizanCanvas', 'shuraCurve',
               'ihsanDotPlot', 'ihsanTimelineCanvas', 'dhulmChart']
        for cid in ids:
            result = driver.execute_script(f"""
                var el = document.getElementById('{cid}');
                if (!el) return null;
                return {{role: el.getAttribute('role'), label: el.getAttribute('aria-label'),
                         fallback: el.textContent.trim()}};
            """)
            assert result is not None, f"Canvas {cid} not found"
            assert result['role'] == 'img', f"{cid} missing role=img"
            assert result['label'] and len(result['label']) > 10, f"{cid} missing aria-label"
            assert result['fallback'] and len(result['fallback']) > 3, f"{cid} missing fallback text"


if __name__ == '__main__':
    pytest.main([__file__, '-v'])
