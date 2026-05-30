import init, {
    pt, ph, ps, hs,
    pt2h, pt2s, pt2v,
    ph2t, ph2s, ph2v,
    ps2t, ps2h, ps2v,
    hs2p, hs2t, hs2v,
    tv, th, ts,
    tv2h, tv2p, tv2s,
    th2p, th2s,
    ts2p, ts2h
} from 'seuif97';

const statusEl = document.getElementById('wasmStatus');

try {
    await init();
    
    statusEl.className = 'status success';
    statusEl.textContent = '✅ seuif97 npm package loaded successfully!';

    const infoEl = document.getElementById('packageInfo');
    infoEl.innerHTML = `
        <p><strong>Package:</strong> seuif97 (npm)</p>
        <p><strong>Functions available:</strong></p>
        <p style="font-size:12px; font-family:monospace;">
            pt, ph, ps, hs, tv, th, ts | 
            pt2h, pt2s, pt2v, pt2x | 
            ph2t, ph2s, ph2v, ph2x |
            ps2t, ps2h, ps2v, ps2x |
            hs2p, hs2t, hs2v, hs2x |
            tv2h, tv2p, tv2s |
            th2p, th2s |
            ts2p, ts2h
        </p>
    `;

    document.getElementById('calcBtn').disabled = false;
    document.getElementById('hsBtn').disabled = false;
} catch (e) {
    statusEl.className = 'status error';
    statusEl.textContent = '❌ Failed to load seuif97: ' + e.message;
    console.error('Package load error:', e);
}

function calculate() {
    const p = parseFloat(document.getElementById('pInput').value);
    const t = parseFloat(document.getElementById('tInput').value);

    const region = pt(p, t, 16);
    const density = pt(p, t, 2);
    const v = pt2v(p, t);
    const h = pt2h(p, t);
    const s = pt2s(p, t);
    const u = pt(p, t, 7);
    const cp = pt(p, t, 8);
    const cv = pt(p, t, 9);
    const w = pt(p, t, 10);
    const k = pt(p, t, 11);
    const viscosity = pt(p, t, 24);
    const thermalCond = pt(p, t, 26);
    const prandtl = pt(p, t, 28);

    document.getElementById('resultPT').style.display = 'block';
    document.getElementById('resultPT').innerHTML = `
        <strong>Results for p=${p} MPa, t=${t}°C:</strong><br>
        <table class="result-table">
            <tr><th>Property</th><th>Value</th><th>Unit</th></tr>
            <tr><td>Region</td><td class="value">${region}</td><td class="unit">-</td></tr>
            <tr><td>Density (ρ)</td><td class="value">${density.toFixed(3)}</td><td class="unit">kg/m³</td></tr>
            <tr><td>Specific Volume (v)</td><td class="value">${v.toFixed(6)}</td><td class="unit">m³/kg</td></tr>
            <tr><td>Specific Enthalpy (h)</td><td class="value">${h.toFixed(3)}</td><td class="unit">kJ/kg</td></tr>
            <tr><td>Specific Entropy (s)</td><td class="value">${s.toFixed(5)}</td><td class="unit">kJ/(kg·K)</td></tr>
            <tr><td>Internal Energy (u)</td><td class="value">${u.toFixed(3)}</td><td class="unit">kJ/kg</td></tr>
            <tr><td>Isobaric Heat Capacity (cp)</td><td class="value">${cp.toFixed(4)}</td><td class="unit">kJ/(kg·K)</td></tr>
            <tr><td>Isochoric Heat Capacity (cv)</td><td class="value">${cv.toFixed(4)}</td><td class="unit">kJ/(kg·K)</td></tr>
            <tr><td>Speed of Sound (w)</td><td class="value">${w.toFixed(2)}</td><td class="unit">m/s</td></tr>
            <tr><td>Isentropic Exponent (k)</td><td class="value">${k.toFixed(6)}</td><td class="unit">-</td></tr>
            <tr><td>Dynamic Viscosity (η)</td><td class="value">${viscosity.toExponential(3)}</td><td class="unit">Pa·s</td></tr>
            <tr><td>Thermal Conductivity (λ)</td><td class="value">${thermalCond.toFixed(4)}</td><td class="unit">W/(m·K)</td></tr>
            <tr><td>Prandtl Number (Pr)</td><td class="value">${prandtl.toFixed(6)}</td><td class="unit">-</td></tr>
        </table>
    `;
}

function calculateHS() {
    const h = parseFloat(document.getElementById('hInput').value);
    const s = parseFloat(document.getElementById('sInput').value);

    const p = hs2p(h, s);
    const t = hs2t(h, s);
    const v = hs2v(h, s);

    if (p < 0 || t < -273.15) {
        document.getElementById('resultHS').style.display = 'block';
        document.getElementById('resultHS').innerHTML = `
            <strong>Invalid (h, s) combination:</strong> h=${h} kJ/kg, s=${s} kJ/(kg·K)<br>
            The given enthalpy-entropy pair is outside the valid region.
        `;
        return;
    }

    const region = pt(p, t, 16);
    const density = pt(p, t, 2);
    const u = pt(p, t, 7);
    const cp = pt(p, t, 8);
    const cv = pt(p, t, 9);
    const w = pt(p, t, 10);
    const k = pt(p, t, 11);
    const viscosity = pt(p, t, 24);
    const thermalCond = pt(p, t, 26);
    const prandtl = pt(p, t, 28);

    document.getElementById('resultHS').style.display = 'block';
    document.getElementById('resultHS').innerHTML = `
        <strong>Results for h=${h} kJ/kg, s=${s} kJ/(kg·K):</strong><br>
        <table class="result-table">
            <tr><th>Property</th><th>Value</th><th>Unit</th></tr>
            <tr><td>Pressure (p)</td><td class="value">${p.toFixed(4)}</td><td class="unit">MPa</td></tr>
            <tr><td>Temperature (t)</td><td class="value">${t.toFixed(2)}</td><td class="unit">°C</td></tr>
            <tr><td>Region</td><td class="value">${region}</td><td class="unit">-</td></tr>
            <tr><td>Density (ρ)</td><td class="value">${density.toFixed(3)}</td><td class="unit">kg/m³</td></tr>
            <tr><td>Specific Volume (v)</td><td class="value">${v.toFixed(6)}</td><td class="unit">m³/kg</td></tr>
            <tr><td>Internal Energy (u)</td><td class="value">${u.toFixed(3)}</td><td class="unit">kJ/kg</td></tr>
            <tr><td>Isobaric Heat Capacity (cp)</td><td class="value">${cp.toFixed(4)}</td><td class="unit">kJ/(kg·K)</td></tr>
            <tr><td>Isochoric Heat Capacity (cv)</td><td class="value">${cv.toFixed(4)}</td><td class="unit">kJ/(kg·K)</td></tr>
            <tr><td>Speed of Sound (w)</td><td class="value">${w.toFixed(2)}</td><td class="unit">m/s</td></tr>
            <tr><td>Isentropic Exponent (k)</td><td class="value">${k.toFixed(6)}</td><td class="unit">-</td></tr>
            <tr><td>Dynamic Viscosity (η)</td><td class="value">${viscosity.toExponential(3)}</td><td class="unit">Pa·s</td></tr>
            <tr><td>Thermal Conductivity (λ)</td><td class="value">${thermalCond.toFixed(4)}</td><td class="unit">W/(m·K)</td></tr>
            <tr><td>Prandtl Number (Pr)</td><td class="value">${prandtl.toFixed(6)}</td><td class="unit">-</td></tr>
        </table>
    `;
}

function loadPreset(region) {
    const presets = {
        1: { p: 3.0, t: 250.0 },
        2: { p: 3.0, t: 500.0 },
        3: { p: 20.0, t: 400.0 }
    };
    const preset = presets[region];
    document.getElementById('pInput').value = preset.p;
    document.getElementById('tInput').value = preset.t;
    calculate();
}

window.calculate = calculate;
window.calculateHS = calculateHS;
window.loadPreset = loadPreset;
