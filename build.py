#!/usr/bin/env python3
"""Build the single-file web app:  python3 build.py [--lint]  ->  dist/index.html

Layout
  data/            generated reference data (elements, densities, stopping tables)
  src/engine/      physics + Monte Carlo (runs in a Web Worker built from a blob)
  src/ui/          interface; files are concatenated in name order into one script
  tools/           scripts that regenerate data/ from libdedx
  tests/           Playwright / node regression tests (see tests/README.md)
"""
import glob, os, subprocess, sys

ROOT = os.path.dirname(os.path.abspath(__file__))
def read(p): return open(os.path.join(ROOT, p), encoding='utf-8').read()
def cat(pattern): return ''.join(read(os.path.relpath(f, ROOT)) for f in sorted(glob.glob(os.path.join(ROOT, pattern))))

TITLE = 'Charged Particle Tracker'
HEAD = f'''<title>{TITLE}</title>
<meta name="claude:viewport-fit-sides" content="cover">
<link rel="preconnect" href="https://fonts.googleapis.com">
<link rel="preconnect" href="https://fonts.gstatic.com" crossorigin>
<link rel="stylesheet" href="https://fonts.googleapis.com/css2?family=Bricolage+Grotesque:opsz,wght@12..96,600;12..96,700&family=IBM+Plex+Sans:wght@400;500;600&family=IBM+Plex+Mono:wght@400;500&display=swap">
<style id="appcss">
'''

def build():
    engine = read('data/elements.js') + read('data/densities.js') + read('data/stopping-tables.js') + cat('src/engine/*.js')
    ui = cat('src/ui/*.js')
    html = (HEAD + read('src/ui/style.css') + '</style>\n' + read('src/ui/body.html') +
            '<script id="engine" type="text/plain">\n' + engine + '\n</script>\n<script id="appjs">\n' + ui + '\n</script>\n')
    os.makedirs(os.path.join(ROOT, 'dist'), exist_ok=True)
    open(os.path.join(ROOT, 'dist/index.html'), 'w', encoding='utf-8').write(html)
    open(os.path.join(ROOT, 'dist/engine.js'), 'w', encoding='utf-8').write(engine)   # for node tests
    open(os.path.join(ROOT, 'dist/ui.js'), 'w', encoding='utf-8').write(ui)
    return engine, ui

def check(engine, ui):
    for name in ('engine', 'ui'):
        r = subprocess.run(['node', '--check', os.path.join(ROOT, f'dist/{name}.js')], capture_output=True, text=True)
        if r.returncode: sys.exit(f'Syntax error in {name}:\n{r.stderr}')

if __name__ == '__main__':
    e, u = build(); check(e, u)
    print(f'dist/index.html  {os.path.getsize(os.path.join(ROOT, "dist/index.html"))/1e3:.0f} kB')
