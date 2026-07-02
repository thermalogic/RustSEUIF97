// node demo_wasm_simple.js
import init, { pt } from '../../pkg/seuif97.js';
import { readFileSync } from 'fs';
import { fileURLToPath } from 'url';
import { dirname, join } from 'path';

const __dirname = dirname(fileURLToPath(import.meta.url));

const wasmPath = join(__dirname, '..', '..', 'pkg', 'seuif97_bg.wasm');
const wasmBuffer = readFileSync(wasmPath);

await init({ module_or_path: wasmBuffer });

const p = 16.0;  // MPa
const t = 535.1; // °C

const h = pt(p, t, 4);

const result = `p = ${p} MPa, t = ${t} °C\nh = ${h.toFixed(3)} kJ/kg`;
console.log(result);