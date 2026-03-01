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

HTML_PATH = 'file:///C:/Users/user/Downloads/metasprint-dose-response/metasprint-dose-response.html'

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
        assert len(options) == 8

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


if __name__ == '__main__':
    pytest.main([__file__, '-v'])
