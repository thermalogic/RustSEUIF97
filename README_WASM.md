# The WASM and NPM Package  

This package is the WebAssembly implementation of the high-speed IAPWS-IF97 package SEUIF97 in Rust. enabling fast and accurate thermodynamic property calculations for water and steam directly in the browser or Node.js. 

## Local WASM 

### Building the WASM file

```bash
cargo build --release --features wasm --target wasm32-unknown-unknown
```

```bash
wasm-bindgen target/wasm32-unknown-unknown/release/seuif97.wasm --out-dir demo_html/pkg --target web
```

### Example

```javascript
import { pt2h, pt2s, ph2t, hs2p } from './pkg/seuif97.js';

const h = pt(16.0, 535.1,4);   // p=16MPa, t=535.1°C → h
const s = pt2s(16.0, 535.1);   // → s
```

* Local WASM example: [./demo_wasm](./demo_wasm/)

```bash
python -m http.server 8080
```

```
http://localhost:8080/
```

## NPM Package 

* NPM package: [seuif97](https://www.npmjs.com/seuif97)

```
npm install seuif97
```

* NPM Package example: [./demo_npm](./demo_npm/)

## T-s Diagram

![](img/ts_wasm.jpg)
