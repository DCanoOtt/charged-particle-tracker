#!/usr/bin/env python3
"""Build the single-file web app:  python3 build.py [--lint] [--debug] [--review]
  dist/index.html     complete standalone document (doctype, charset, viewport): open it, host it, or ship it
  dist/artifact.html  the same page as a bare fragment (title/style/body/scripts) for the Claude Artifact tool,
                      which wraps it in its own document skeleton at publish time

Layout
  data/            generated reference data (elements, densities, stopping tables)
  src/engine/      physics + Monte Carlo (runs in a Web Worker built from a blob)
  src/ui/          interface; files are concatenated in name order into one script
  tools/           scripts that regenerate data/ from libdedx
  tests/           Playwright / node regression tests (see tests/README.md)

Every module starts with a banner line  /* ===== src/ui/49-gamma.js ===== */  (search for it, or for the function name).
--review additionally writes dist/review.txt: all the code (README, build, engine and UI modules, help text, styles, body, types, tests) in one text file for a reviewer or an AI assistant,
with the large generated data tables replaced by one-line placeholders (see REVIEW_STUBS).
--debug additionally writes dist/index.debug.html: the same page with an inline source map (original files and line numbers in the browser's debugger)
and  //# sourceURL  names for the UI script (cpt-ui.js) and the worker (cpt-engine.js, always present).
"""
import glob, json, os, subprocess, sys

ROOT = os.path.dirname(os.path.abspath(__file__))
def read(p): return open(os.path.join(ROOT, p), encoding='utf-8').read()
def banner(name): return f'/* ===== {name} ===== */\n'
def module(name): s = read(name); return banner(name) + (s if s.endswith('\n') else s + '\n')
def files(pattern): return [os.path.relpath(f, ROOT).replace(os.sep, '/') for f in sorted(glob.glob(os.path.join(ROOT, pattern)))]
def cat(pattern): return ''.join(module(f) for f in files(pattern))

import re
def flatten(name):
    """An ES module of the engine as a plain script: imports and export keywords removed (all top-level names are unique, checked in build())."""
    s = read(name)
    s = re.sub(r'^import [^;]*;[ \t]*\n', '', s, flags=re.M)
    s = re.sub(r'^export (?=(?:async )?(?:function|const|let|class)\b)', '', s, flags=re.M)
    if re.search(r'^(?:import|export)\b', s, flags=re.M): sys.exit(f'{name}: unsupported import/export syntax')
    return s
def top_names(s):
    names = []
    for m in re.finditer(r'^(?:function\s+(\w+)|(?:const|let)\s+([^=;\n]+?)\s*(?==))', s, re.M):
        names += [m.group(1)] if m.group(1) else [x.strip() for x in m.group(2).split(',') if re.fullmatch(r'\w+', x.strip())]
    return names

B64 = 'ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789+/'
def vlq(n):
    n = (-n << 1) | 1 if n < 0 else n << 1; out = ''
    while True:
        d = n & 31; n >>= 5
        if n: d |= 32
        out += B64[d]
        if not n: return out
def sourcemap(names, texts):
    """V3 source map of a concatenation in which every module is preceded by a banner line (line to line mapping)."""
    maps, prev_src, prev_line = [], 0, 0
    for i, s in enumerate(texts):
        n = s.count('\n') + (0 if s.endswith('\n') else 1)
        maps.append(';')                                             # banner line: no mapping
        for k in range(n):
            maps.append(vlq(0) + vlq(i - prev_src) + vlq(k - prev_line) + vlq(0) + ';'); prev_src, prev_line = i, k
    return {'version': 3, 'file': 'cpt-ui.js', 'sources': names, 'sourcesContent': texts, 'names': [], 'mappings': ''.join(maps)}

TITLE = 'Charged Particle Tracker'
FONTS = '''<link rel="preconnect" href="https://fonts.googleapis.com">
<link rel="preconnect" href="https://fonts.gstatic.com" crossorigin>
<link rel="stylesheet" href="https://fonts.googleapis.com/css2?family=Bricolage+Grotesque:opsz,wght@12..96,600;12..96,700&family=IBM+Plex+Sans:wght@400;500;600&family=IBM+Plex+Mono:wght@400;500&display=swap">
<style id="appcss">
'''
HEAD = f'''<title>{TITLE}</title>
<meta name="claude:viewport-fit-sides" content="cover">
''' + FONTS                                   # fragment head (the Artifact tool adds doctype/charset/viewport itself)
DOC_HEAD = f'''<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>{TITLE}</title>
''' + FONTS

def build():
    eng = files('src/engine/*.js'); seen = {}
    for n in eng:
        for nm in top_names(flatten(n)):
            if nm in seen: sys.exit(f'engine: "{nm}" is declared in {seen[nm]} and in {n} (the engine modules share one scope; names must be unique)')
            seen[nm] = n
    engine = read('data/elements.js') + read('data/densities.js') + read('data/stopping-tables.js') + ''.join(banner(n) + flatten(n) + ('' if flatten(n).endswith('\n') else '\n') for n in eng)
    help_html = ''.join(read(n) for n in files('src/help/*.html'))              # the text of the instructions: one file per chapter / section
    ui_names = ['data/xray.js', 'src/help/*.html'] + files('src/ui/*.js'); ui_texts = [read('data/xray.js'), 'const HELP_HTML = ' + json.dumps(help_html, ensure_ascii=False).replace('</', '<\\/') + ';\n'] + [read(n) for n in ui_names[2:]]
    ui = ''.join(banner(n) + (s if s.endswith('\n') else s + '\n') for n, s in zip(ui_names, ui_texts))
    css, body = read('src/ui/style.css'), read('src/ui/body.html')
    scripts = '<script id="engine" type="text/plain">\n' + engine + '\n</script>\n<script id="appjs">\n' + ui + '\n</script>\n'
    fragment = HEAD + css + '</style>\n' + body + scripts
    standalone = DOC_HEAD + css + '</style>\n</head>\n<body>\n' + body + scripts + '</body>\n</html>\n'
    os.makedirs(os.path.join(ROOT, 'dist'), exist_ok=True)
    open(os.path.join(ROOT, 'dist/index.html'), 'w', encoding='utf-8').write(standalone)
    open(os.path.join(ROOT, 'dist/artifact.html'), 'w', encoding='utf-8').write(fragment)
    open(os.path.join(ROOT, 'dist/engine.js'), 'w', encoding='utf-8').write(engine)   # for node tests
    open(os.path.join(ROOT, 'dist/ui.js'), 'w', encoding='utf-8').write(ui)
    if '--debug' in sys.argv:
        import base64
        sm = base64.b64encode(json.dumps(sourcemap(ui_names, [s if s.endswith('\n') else s + '\n' for s in ui_texts]), ensure_ascii=False).encode('utf-8')).decode()
        dbg = scripts.replace('<script id="appjs">\n' + ui, '<script id="appjs">\n' + ui + '\n//# sourceURL=cpt-ui.js\n//# sourceMappingURL=data:application/json;base64,' + sm + '\n')
        open(os.path.join(ROOT, 'dist/index.debug.html'), 'w', encoding='utf-8').write(DOC_HEAD + css + '</style>\n</head>\n<body>\n' + body + dbg + '</body>\n</html>\n')
    return engine, ui

REVIEW_STUBS = ['data/elements.js', 'data/densities.js', 'data/stopping-tables.js', 'data/xray.js', 'src/ui/48-half-lives.js', 'src/ui/48-decay-icrp.js']     # generated tables
def review():
    parts = []
    def add(name, text): parts.append(f'\n\n################ {name} ################\n' + text)
    for name in ['README.md', 'build.py', 'jsconfig.json', 'src/types.d.ts', 'src/ui/body.html', 'src/ui/style.css'] + files('data/*.js') + files('src/engine/*.js') + files('src/ui/*.js') + files('src/help/*.html') + files('tests/*.js') + ['tests/test_app.py'] + files('tools/*.py'):
        if name in REVIEW_STUBS:
            s = read(name); first = s.split('\n', 1)[0][:300]
            add(name, f'[generated data table, {len(s) / 1e3:.0f} kB, {s.count(chr(10)) + 1} lines, omitted here. First line of the file: {first}]\n')
        else: add(name, read(name))
    out = '# Charged Particle Tracker: complete source for review (generated by build.py --review)\n# Sections are separated by ################ <path> ################. Generated tables are replaced by a placeholder; the page needs them to run.\n' + ''.join(parts)
    open(os.path.join(ROOT, 'dist/review.txt'), 'w', encoding='utf-8').write(out)
    return len(out)

def check(engine, ui):
    for name in ('engine', 'ui'):
        r = subprocess.run(['node', '--check', os.path.join(ROOT, f'dist/{name}.js')], capture_output=True, text=True)
        if r.returncode: sys.exit(f'Syntax error in {name}:\n{r.stderr}')

if __name__ == '__main__':
    e, u = build(); check(e, u)
    if '--review' in sys.argv: print(f'dist/review.txt  {review() / 1e3:.0f} kB')
    for f in ('index.html', 'artifact.html') + (('index.debug.html',) if '--debug' in sys.argv else ()): print(f'dist/{f}  {os.path.getsize(os.path.join(ROOT, "dist", f))/1e3:.0f} kB')
