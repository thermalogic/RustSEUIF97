// node demo_wasm_npm_simple.js
//import init, { pt } from './node_modules/seuif97/pkg/seuif97.js'; // node_modules m
import init, { pt } from '../../pkg/seuif97.js';  // local pkg

import { readFileSync } from 'fs';
import { fileURLToPath } from 'url';
import { dirname, join } from 'path';

const __dirname = dirname(fileURLToPath(import.meta.url));

//const wasmPath = join(__dirname, 'node_modules', 'seuif97', 'pkg', 'seuif97_bg.wasm'); //node_modules 
const wasmPath = join(__dirname, '..', '..', 'pkg', 'seuif97_bg.wasm'); //local pkg

const wasmBuffer = readFileSync(wasmPath);

await init({ module_or_path: wasmBuffer });

const p = 16.0;
const t = 535.1;
const h = pt(p, t, 4);

console.log(`p = ${p} MPa, t = ${t} °C`);
console.log(`h = ${h.toFixed(3)} kJ/kg`);