/* @ts-self-types="./seuif97.d.ts" */

/**
 * @param {number} h
 * @param {number} s
 * @param {number} o_id
 * @returns {number}
 */
export function hs(h, s, o_id) {
    const ret = wasm.hs(h, s, o_id);
    return ret;
}

/**
 * @param {number} h
 * @param {number} s
 * @returns {number}
 */
export function hs2p(h, s) {
    const ret = wasm.hs2p(h, s);
    return ret;
}

/**
 * @param {number} h
 * @param {number} s
 * @returns {number}
 */
export function hs2t(h, s) {
    const ret = wasm.hs2t(h, s);
    return ret;
}

/**
 * @param {number} h
 * @param {number} s
 * @returns {number}
 */
export function hs2v(h, s) {
    const ret = wasm.hs2v(h, s);
    return ret;
}

/**
 * @param {number} h
 * @param {number} s
 * @returns {number}
 */
export function hs2x(h, s) {
    const ret = wasm.hs2x(h, s);
    return ret;
}

/**
 * @param {number} h
 * @param {number} x
 * @param {number} o_id
 * @returns {number}
 */
export function hx(h, x, o_id) {
    const ret = wasm.hx(h, x, o_id);
    return ret;
}

/**
 * @param {number} h
 * @param {number} x
 * @returns {number}
 */
export function hx2p(h, x) {
    const ret = wasm.hx2p(h, x);
    return ret;
}

/**
 * @param {number} h
 * @param {number} x
 * @returns {number}
 */
export function hx2s(h, x) {
    const ret = wasm.hx2s(h, x);
    return ret;
}

/**
 * @param {number} h
 * @param {number} x
 * @returns {number}
 */
export function hx2t(h, x) {
    const ret = wasm.hx2t(h, x);
    return ret;
}

/**
 * @param {number} h
 * @param {number} x
 * @returns {number}
 */
export function hx2v(h, x) {
    const ret = wasm.hx2v(h, x);
    return ret;
}

/**
 * @param {number} p
 * @param {number} h
 * @param {number} o_id
 * @returns {number}
 */
export function ph(p, h, o_id) {
    const ret = wasm.ph(p, h, o_id);
    return ret;
}

/**
 * @param {number} p
 * @param {number} h
 * @returns {number}
 */
export function ph2s(p, h) {
    const ret = wasm.ph2s(p, h);
    return ret;
}

/**
 * @param {number} p
 * @param {number} h
 * @returns {number}
 */
export function ph2t(p, h) {
    const ret = wasm.ph2t(p, h);
    return ret;
}

/**
 * @param {number} p
 * @param {number} h
 * @returns {number}
 */
export function ph2v(p, h) {
    const ret = wasm.ph2v(p, h);
    return ret;
}

/**
 * @param {number} p
 * @param {number} h
 * @returns {number}
 */
export function ph2x(p, h) {
    const ret = wasm.ph2x(p, h);
    return ret;
}

/**
 * @param {number} p
 * @param {number} s
 * @param {number} o_id
 * @returns {number}
 */
export function ps(p, s, o_id) {
    const ret = wasm.ps(p, s, o_id);
    return ret;
}

/**
 * @param {number} p
 * @param {number} s
 * @returns {number}
 */
export function ps2h(p, s) {
    const ret = wasm.ps2h(p, s);
    return ret;
}

/**
 * @param {number} p
 * @param {number} s
 * @returns {number}
 */
export function ps2t(p, s) {
    const ret = wasm.ps2t(p, s);
    return ret;
}

/**
 * @param {number} p
 * @param {number} s
 * @returns {number}
 */
export function ps2v(p, s) {
    const ret = wasm.ps2v(p, s);
    return ret;
}

/**
 * @param {number} p
 * @param {number} s
 * @returns {number}
 */
export function ps2x(p, s) {
    const ret = wasm.ps2x(p, s);
    return ret;
}

/**
 * @param {number} p
 * @param {number} t
 * @param {number} o_id
 * @returns {number}
 */
export function pt(p, t, o_id) {
    const ret = wasm.pt(p, t, o_id);
    return ret;
}

/**
 * @param {number} p
 * @param {number} t
 * @returns {number}
 */
export function pt2h(p, t) {
    const ret = wasm.pt2h(p, t);
    return ret;
}

/**
 * @param {number} p
 * @param {number} t
 * @returns {number}
 */
export function pt2s(p, t) {
    const ret = wasm.pt2s(p, t);
    return ret;
}

/**
 * @param {number} p
 * @param {number} t
 * @returns {number}
 */
export function pt2v(p, t) {
    const ret = wasm.pt2v(p, t);
    return ret;
}

/**
 * @param {number} p
 * @param {number} t
 * @returns {number}
 */
export function pt2x(p, t) {
    const ret = wasm.pt2x(p, t);
    return ret;
}

/**
 * @param {number} p
 * @param {number} v
 * @param {number} o_id
 * @returns {number}
 */
export function pv(p, v, o_id) {
    const ret = wasm.pv(p, v, o_id);
    return ret;
}

/**
 * @param {number} p
 * @param {number} v
 * @returns {number}
 */
export function pv2h(p, v) {
    const ret = wasm.pv2h(p, v);
    return ret;
}

/**
 * @param {number} p
 * @param {number} v
 * @returns {number}
 */
export function pv2s(p, v) {
    const ret = wasm.pv2s(p, v);
    return ret;
}

/**
 * @param {number} p
 * @param {number} v
 * @returns {number}
 */
export function pv2t(p, v) {
    const ret = wasm.pv2t(p, v);
    return ret;
}

/**
 * @param {number} p
 * @param {number} v
 * @returns {number}
 */
export function pv2x(p, v) {
    const ret = wasm.pv2x(p, v);
    return ret;
}

/**
 * @param {number} p
 * @param {number} x
 * @param {number} o_id
 * @returns {number}
 */
export function px(p, x, o_id) {
    const ret = wasm.px(p, x, o_id);
    return ret;
}

/**
 * @param {number} p
 * @param {number} x
 * @returns {number}
 */
export function px2h(p, x) {
    const ret = wasm.px2h(p, x);
    return ret;
}

/**
 * @param {number} p
 * @param {number} x
 * @returns {number}
 */
export function px2s(p, x) {
    const ret = wasm.px2s(p, x);
    return ret;
}

/**
 * @param {number} p
 * @param {number} x
 * @returns {number}
 */
export function px2t(p, x) {
    const ret = wasm.px2t(p, x);
    return ret;
}

/**
 * @param {number} p
 * @param {number} x
 * @returns {number}
 */
export function px2v(p, x) {
    const ret = wasm.px2v(p, x);
    return ret;
}

/**
 * @param {number} s
 * @param {number} x
 * @param {number} o_id
 * @returns {number}
 */
export function sx(s, x, o_id) {
    const ret = wasm.sx(s, x, o_id);
    return ret;
}

/**
 * @param {number} s
 * @param {number} x
 * @returns {number}
 */
export function sx2h(s, x) {
    const ret = wasm.sx2h(s, x);
    return ret;
}

/**
 * @param {number} s
 * @param {number} x
 * @returns {number}
 */
export function sx2p(s, x) {
    const ret = wasm.sx2p(s, x);
    return ret;
}

/**
 * @param {number} s
 * @param {number} x
 * @returns {number}
 */
export function sx2t(s, x) {
    const ret = wasm.sx2t(s, x);
    return ret;
}

/**
 * @param {number} s
 * @param {number} x
 * @returns {number}
 */
export function sx2v(s, x) {
    const ret = wasm.sx2v(s, x);
    return ret;
}

/**
 * @param {number} t
 * @param {number} h
 * @param {number} o_id
 * @returns {number}
 */
export function th(t, h, o_id) {
    const ret = wasm.th(t, h, o_id);
    return ret;
}

/**
 * @param {number} t
 * @param {number} h
 * @returns {number}
 */
export function th2p(t, h) {
    const ret = wasm.th2p(t, h);
    return ret;
}

/**
 * @param {number} t
 * @param {number} h
 * @returns {number}
 */
export function th2s(t, h) {
    const ret = wasm.th2s(t, h);
    return ret;
}

/**
 * @param {number} t
 * @param {number} h
 * @returns {number}
 */
export function th2v(t, h) {
    const ret = wasm.th2v(t, h);
    return ret;
}

/**
 * @param {number} t
 * @param {number} h
 * @returns {number}
 */
export function th2x(t, h) {
    const ret = wasm.th2x(t, h);
    return ret;
}

/**
 * @param {number} t
 * @param {number} s
 * @param {number} o_id
 * @returns {number}
 */
export function ts(t, s, o_id) {
    const ret = wasm.ts(t, s, o_id);
    return ret;
}

/**
 * @param {number} t
 * @param {number} s
 * @returns {number}
 */
export function ts2h(t, s) {
    const ret = wasm.ts2h(t, s);
    return ret;
}

/**
 * @param {number} t
 * @param {number} s
 * @returns {number}
 */
export function ts2p(t, s) {
    const ret = wasm.ts2p(t, s);
    return ret;
}

/**
 * @param {number} t
 * @param {number} s
 * @returns {number}
 */
export function ts2v(t, s) {
    const ret = wasm.ts2v(t, s);
    return ret;
}

/**
 * @param {number} t
 * @param {number} s
 * @returns {number}
 */
export function ts2x(t, s) {
    const ret = wasm.ts2x(t, s);
    return ret;
}

/**
 * @param {number} t
 * @param {number} v
 * @param {number} o_id
 * @returns {number}
 */
export function tv(t, v, o_id) {
    const ret = wasm.tv(t, v, o_id);
    return ret;
}

/**
 * @param {number} t
 * @param {number} v
 * @returns {number}
 */
export function tv2h(t, v) {
    const ret = wasm.tv2h(t, v);
    return ret;
}

/**
 * @param {number} t
 * @param {number} v
 * @returns {number}
 */
export function tv2p(t, v) {
    const ret = wasm.tv2p(t, v);
    return ret;
}

/**
 * @param {number} t
 * @param {number} v
 * @returns {number}
 */
export function tv2s(t, v) {
    const ret = wasm.tv2s(t, v);
    return ret;
}

/**
 * @param {number} t
 * @param {number} v
 * @returns {number}
 */
export function tv2x(t, v) {
    const ret = wasm.tv2x(t, v);
    return ret;
}

/**
 * @param {number} t
 * @param {number} x
 * @param {number} o_id
 * @returns {number}
 */
export function tx(t, x, o_id) {
    const ret = wasm.tx(t, x, o_id);
    return ret;
}

/**
 * @param {number} t
 * @param {number} x
 * @returns {number}
 */
export function tx2h(t, x) {
    const ret = wasm.tx2h(t, x);
    return ret;
}

/**
 * @param {number} t
 * @param {number} x
 * @returns {number}
 */
export function tx2p(t, x) {
    const ret = wasm.tx2p(t, x);
    return ret;
}

/**
 * @param {number} t
 * @param {number} x
 * @returns {number}
 */
export function tx2s(t, x) {
    const ret = wasm.tx2s(t, x);
    return ret;
}

/**
 * @param {number} t
 * @param {number} x
 * @returns {number}
 */
export function tx2v(t, x) {
    const ret = wasm.tx2v(t, x);
    return ret;
}
function __wbg_get_imports() {
    const import0 = {
        __proto__: null,
        __wbindgen_init_externref_table: function() {
            const table = wasm.__wbindgen_externrefs;
            const offset = table.grow(4);
            table.set(0, undefined);
            table.set(offset + 0, undefined);
            table.set(offset + 1, null);
            table.set(offset + 2, true);
            table.set(offset + 3, false);
        },
    };
    return {
        __proto__: null,
        "./seuif97_bg.js": import0,
    };
}

let wasmModule, wasmInstance, wasm;
function __wbg_finalize_init(instance, module) {
    wasmInstance = instance;
    wasm = instance.exports;
    wasmModule = module;
    wasm.__wbindgen_start();
    return wasm;
}

async function __wbg_load(module, imports) {
    if (typeof Response === 'function' && module instanceof Response) {
        if (typeof WebAssembly.instantiateStreaming === 'function') {
            try {
                return await WebAssembly.instantiateStreaming(module, imports);
            } catch (e) {
                const validResponse = module.ok && expectedResponseType(module.type);

                if (validResponse && module.headers.get('Content-Type') !== 'application/wasm') {
                    console.warn("`WebAssembly.instantiateStreaming` failed because your server does not serve Wasm with `application/wasm` MIME type. Falling back to `WebAssembly.instantiate` which is slower. Original error:\n", e);

                } else { throw e; }
            }
        }

        const bytes = await module.arrayBuffer();
        return await WebAssembly.instantiate(bytes, imports);
    } else {
        const instance = await WebAssembly.instantiate(module, imports);

        if (instance instanceof WebAssembly.Instance) {
            return { instance, module };
        } else {
            return instance;
        }
    }

    function expectedResponseType(type) {
        switch (type) {
            case 'basic': case 'cors': case 'default': return true;
        }
        return false;
    }
}

function initSync(module) {
    if (wasm !== undefined) return wasm;


    if (module !== undefined) {
        if (Object.getPrototypeOf(module) === Object.prototype) {
            ({module} = module)
        } else {
            console.warn('using deprecated parameters for `initSync()`; pass a single object instead')
        }
    }

    const imports = __wbg_get_imports();
    if (!(module instanceof WebAssembly.Module)) {
        module = new WebAssembly.Module(module);
    }
    const instance = new WebAssembly.Instance(module, imports);
    return __wbg_finalize_init(instance, module);
}

async function __wbg_init(module_or_path) {
    if (wasm !== undefined) return wasm;


    if (module_or_path !== undefined) {
        if (Object.getPrototypeOf(module_or_path) === Object.prototype) {
            ({module_or_path} = module_or_path)
        } else {
            console.warn('using deprecated parameters for the initialization function; pass a single object instead')
        }
    }

    if (module_or_path === undefined) {
        module_or_path = new URL('seuif97_bg.wasm', import.meta.url);
    }
    const imports = __wbg_get_imports();

    if (typeof module_or_path === 'string' || (typeof Request === 'function' && module_or_path instanceof Request) || (typeof URL === 'function' && module_or_path instanceof URL)) {
        module_or_path = fetch(module_or_path);
    }

    const { instance, module } = await __wbg_load(await module_or_path, imports);

    return __wbg_finalize_init(instance, module);
}

export { initSync, __wbg_init as default };
