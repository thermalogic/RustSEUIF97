# The WASM and NPM Package  

![npm version](https://img.shields.io/npm/v/seuif97)![NPM Downloads](https://img.shields.io/npm/dm/seuif97)![NPM Downloads](https://img.shields.io/npm/dt/seuif97)

The WebAssembly (ES modules) implementation of the high-speed IAPWS-IF97 package SEUIF97, written in Rust, enabling fast and accurate thermodynamic property calculations for water and steam in the browser.

This package supports **12 distinct input state pairs** for calculating **36 thermodynamic, transport, and derived properties**, plus **thermodynamic process functions** for isentropic enthalpy drop and efficiency calculations.

## Building the WASM with wasm-bindgen 

```bash
cargo build --release --features wasm --target wasm32-unknown-unknown
```

```bash
wasm-bindgen target/wasm32-unknown-unknown/release/seuif97.wasm --out-dir pkg --target web
```

## Property Calculation Functions

The package provides two types of API for property calculation.

### Universal Property Functions

The following 12 input pairs are implemented:

```bash
  (p,t), (p,h), (p,s), (p,v)
  (h,s)
  (t,h), (t,s), (t,v)
  (h,x), (t,x), (v,x), (s,x)
```            
Each function accepts an input pair, an output property ID ([o_id](#properties)), and an `optional` region parameter for faster computation.
For example: the input pair (p,t): `pt(p,t,o_id)`, `pt(p,t,(o_id,region))`

**Note:**
 * The `region` parameter is Rust-only. C, Python and WASM bindings support the `o_id` form only.
 * Only `linearly` related thermodynamic properties are calculable in the `wet` steam region.

### Direct Property Functions

Function naming convention: `{input1}{input2}2{output}`.

| Input | Outputs | Input | Outputs | Input | Outputs |
|-------|---------|-------|---------|-------|---------|
| (p,t) | h,s,v,x | (t,h) | p,s,v,x | (p,x) | t,h,s,v |
| (p,h) | t,s,v,x | (t,s) | p,h,v,x | (t,x) | p,h,s,v |
| (p,s) | t,h,v,x | (t,v) | p,h,s,x | (h,x) | p,t,s,v |
| (p,v) | t,h,s,x | (h,s) | p,t,v,x | (s,x) | p,t,h,v |

Total: 48 functions(e.g. `pt2h(p,t)`, `ph2t(p,h)`, `hs2p(h,s)`)

### Thermodynamic Process Functions

The following thermodynamic process functions are implemented:

- `ishd(pi, ti, pe)`: isentropic enthalpy drop for steam expansion (kJ/kg)
- `ief(pi, ti, pe, te)`: isentropic efficiency for superheated steam expansion (%)

## WASM(ES Modules)

* [./demo_using_lib/demo_wasm_npm](./demo_using_lib/demo_wasm_npm/)

```bash
python -m http.server 8080
```

```
http://localhost:8080/
```

## NPM Package 

* NPM package: [seuif97](https://www.npmjs.com/seuif97)

* NPM Package example: [./demo_using_lib/demo_wasm_npm](./demo_using_lib/demo_wasm_npm/)

```bash
npm install
```

```bash
npm install seuif97
```

```javascript
// node demo_npm_simple.js
import init, { pt } from './node_modules/seuif97/pkg/seuif97.js';
import { readFileSync } from 'fs';
import { fileURLToPath } from 'url';
import { dirname, join } from 'path';

const __dirname = dirname(fileURLToPath(import.meta.url));

const wasmPath = join(__dirname, 'node_modules', 'seuif97', 'pkg', 'seuif97_bg.wasm');
const wasmBuffer = readFileSync(wasmPath);

await init({ module_or_path: wasmBuffer });

const p = 16.0;
const t = 535.1;
const h = pt(p, t, 4);

console.log(`p = ${p} MPa, t = ${t} °C`);
console.log(`h = ${h.toFixed(3)} kJ/kg`);
```

## T-s Diagram

![](img/ts_wasm.jpg)
