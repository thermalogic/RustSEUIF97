# The WASM and NPM Package  

The WebAssembly (ES modules) implementation of the high-speed IAPWS-IF97 package SEUIF97, written in Rust, enabling fast and accurate thermodynamic property calculations for water and steam in the browser.

## Building the WASM with wasm-bindgen 

```bash
cargo build --release --features wasm --target wasm32-unknown-unknown
```

* Local package: `pkg`

```bash
wasm-bindgen target/wasm32-unknown-unknown/release/seuif97.wasm --out-dir pkg --target web
```

## WASM(ES Modules)

* NPM package: [seuif97](https://www.npmjs.com/seuif97)

```bash
npm install
```
```bash
npm install seuif97
```

## Examples

* [./demo_using_lib/demo_wasm](./demo_using_lib/demo_wasm/)

### JavaScript Example

```javascript
// node demo_wasm_simple.js
// import init, { pt } from './node_modules/seuif97/pkg/seuif97.js'; // node_modules
import init, { pt } from '../../pkg/seuif97.js';  // local pkg
import { readFileSync } from 'fs';
import { fileURLToPath } from 'url';
import { dirname, join } from 'path';

const __dirname = dirname(fileURLToPath(import.meta.url));

// const wasmPath = join(__dirname, 'node_modules', 'seuif97', 'pkg', 'seuif97_bg.wasm'); // node_modules
const wasmPath = join(__dirname, '..', '..', 'pkg', 'seuif97_bg.wasm'); //local pkg
const wasmBuffer = readFileSync(wasmPath);

await init({ module_or_path: wasmBuffer });

const p = 16.0;
const t = 535.1;
const h = pt(p, t, 4);

console.log(`p = ${p} MPa, t = ${t} °C`);
console.log(`h = ${h.toFixed(3)} kJ/kg`);
```

### HTML Example

* [./demo_using_lib/demo_wasm/demo_wasm_simple.html](./demo_using_lib/demo_wasm/demo_wasm_simple.html)

```bash
python -m http.server 8080
```

```
http://localhost:8080/
```

### T-s Diagram

* [./demo_using_lib/demo_wasm/demo_wasm_simple.html](./demo_using_lib/demo_wasm/T-S.html)

![](img/ts_wasm.jpg)

