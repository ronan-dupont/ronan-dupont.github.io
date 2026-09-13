/* Dependency-free bundling. Vendor sources remain unchanged and licensed. */
const fs = require('node:fs');
const path = require('node:path');
const root = path.resolve(__dirname, '..');
const sources = [
  'assets/js/vendor/jquery/jquery-1.12.4.min.js',
  'assets/js/plugins/jquery.fitvids.js',
  'assets/js/plugins/jquery.magnific-popup.js',
  'assets/js/_main.js'
];
const output = sources.map(file => fs.readFileSync(path.join(root, file), 'utf8')).join('\n;\n');
fs.writeFileSync(path.join(root, 'assets/js/main.min.js'), output);
console.log(`Built assets/js/main.min.js from ${sources.length} preserved sources (${Buffer.byteLength(output)} bytes).`);
