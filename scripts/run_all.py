"""End-to-end reproducibility driver.

Default mode is fast: rebuild summaries, panel CSVs, and figures from bundled data.
Use --recompute to rerun deterministic theory/dynamics/validation first.
Use --stochastic to rerun the R=64 stochastic pilot.
"""
from pathlib import Path
import subprocess, sys, argparse
ROOT = Path(__file__).resolve().parents[1]
parser = argparse.ArgumentParser()
parser.add_argument('--recompute', action='store_true')
parser.add_argument('--stochastic', action='store_true')
args = parser.parse_args()

def run(script):
    subprocess.run([sys.executable, str(ROOT/'scripts'/script)], cwd=ROOT, check=True)

if args.recompute:
    run('run_theory_static.py')
    run('run_dynamic.py')
    run('analyze_dynamic.py')
    run('run_validation.py')
if args.stochastic:
    run('run_stochastic_pilot.py')
run('build_summaries.py')
run('make_figures.py')
run('export_panel_csvs.py')
print('Package rebuild complete.')
