// Copies the single-file web app (built by ../build.py) into ./app/ so that Electron can load it.
const fs = require('fs'), path = require('path');
const src = path.join(__dirname, '..', 'dist', 'index.html'), dst = path.join(__dirname, 'app');
if (!fs.existsSync(src)) { console.error('dist/index.html not found — run "python3 build.py" in the repository root first.'); process.exit(1); }
fs.mkdirSync(dst, { recursive: true }); fs.copyFileSync(src, path.join(dst, 'index.html'));
console.log('app/index.html updated (' + (fs.statSync(path.join(dst, 'index.html')).size / 1024).toFixed(0) + ' kB)');
