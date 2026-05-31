# The WASM and NPM Package  

![npm version](https://img.shields.io/npm/v/seuif97)![NPM Downloads](https://img.shields.io/npm/dm/seuif97)![NPM Downloads](https://img.shields.io/npm/dt/seuif97)

The WebAssembly implementation of the high-speed IAPWS-IF97 package SEUIF97 in Rust, enabling fast and accurate thermodynamic property calculations for water and steam directly in the browser or Node.js. 

## Building the WASM file

```bash
cargo build --release --features wasm --target wasm32-unknown-unknown
```

```bash
wasm-bindgen target/wasm32-unknown-unknown/release/seuif97.wasm --out-dir demo_wasm/pkg --target web
```

## Basic Usage (ES Modules)

```javascript
import { pt } from './pkg/seuif97.js';
await init();

const p = 16.0;  // MPa
const t = 535.1; // °C

const h = pt(p, t, 4);

console.log(`p = ${p} MPa, t = ${t} °C`);
console.log(`h = ${h.toFixed(3)} kJ/kg`);
```
## Using in Web Browsers

* Example: [./demo_wasm](./demo_wasm/)

```bash
python -m http.server 8080
```

```
http://localhost:8080/
```

## Using in Node.js With NPM Package 

* NPM package: [seuif97](https://www.npmjs.com/seuif97)

```bash
npm install seuif97
```

* NPM Package example: [./demo_npm](./demo_npm/)

```javascript
import init, { pt2h } from 'seuif97';

await init();

const p = 3.0;    // MPa
const t = 250.0;  // °C

const h = pt2h(p, t);

console.log(`p = ${p} MPa, t = ${t} °C`);
console.log(`h = ${h.toFixed(5)} kJ/kg`);
```

## T-s Diagram

![](img/ts_wasm.jpg)
