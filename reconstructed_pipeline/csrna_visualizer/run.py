"""Local-only launcher, independent of the analysis pipeline."""
import argparse
from pathlib import Path
import subprocess
import sys

p = argparse.ArgumentParser(description='Open csRNA pipeline data locally. Reads pipeline JSON; never runs the pipeline.')
p.add_argument('--report', help='Path to a pipeline JSON report; can also be chosen in the app')
p.add_argument('--port', type=int, default=8501)
args = p.parse_args()
if not 1024 <= args.port <= 65535:
    p.error('Choose a port between 1024 and 65535')
command = [sys.executable, '-m', 'streamlit', 'run', str(Path(__file__).resolve().with_name('app.py')),
           '--server.address=127.0.0.1', f'--server.port={args.port}', '--server.headless=true',
           '--browser.gatherUsageStats=false', '--theme.base=light', '--theme.primaryColor=#487C97']
if args.report:
    command += ['--', '--report', str(Path(args.report).expanduser().resolve())]
raise SystemExit(subprocess.call(command))
