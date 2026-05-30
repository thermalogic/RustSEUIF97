# The WASM binding 

## Building the WASM file

```bash
cargo build --release --features wasm --target wasm32-unknown-unknown
```
```bash
wasm-bindgen target/wasm32-unknown-unknown/release/seuif97.wasm --out-dir demo_html/pkg --target web
```

## Example

```javascript
import { pt2h, pt2s, ph2t, hs2p } from './pkg/seuif97.js';

const h = pt(16.0, 535.1,4);   // p=16MPa, t=535.1°C → h
const s = pt2s(16.0, 535.1);   // → s
```

## Example Web

* [./demo_html/](./demo_html/)

```bash
python -m http.server 8080
```

```
http://localhost:8080/
```
