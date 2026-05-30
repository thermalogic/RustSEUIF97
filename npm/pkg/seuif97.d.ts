/* tslint:disable */
/* eslint-disable */

export function hs(h: number, s: number, o_id: number): number;

export function hs2p(h: number, s: number): number;

export function hs2t(h: number, s: number): number;

export function hs2v(h: number, s: number): number;

export function hs2x(h: number, s: number): number;

export function hx(h: number, x: number, o_id: number): number;

export function hx2p(h: number, x: number): number;

export function hx2s(h: number, x: number): number;

export function hx2t(h: number, x: number): number;

export function hx2v(h: number, x: number): number;

export function ph(p: number, h: number, o_id: number): number;

export function ph2s(p: number, h: number): number;

export function ph2t(p: number, h: number): number;

export function ph2v(p: number, h: number): number;

export function ph2x(p: number, h: number): number;

export function ps(p: number, s: number, o_id: number): number;

export function ps2h(p: number, s: number): number;

export function ps2t(p: number, s: number): number;

export function ps2v(p: number, s: number): number;

export function ps2x(p: number, s: number): number;

export function pt(p: number, t: number, o_id: number): number;

export function pt2h(p: number, t: number): number;

export function pt2s(p: number, t: number): number;

export function pt2v(p: number, t: number): number;

export function pt2x(p: number, t: number): number;

export function pv(p: number, v: number, o_id: number): number;

export function pv2h(p: number, v: number): number;

export function pv2s(p: number, v: number): number;

export function pv2t(p: number, v: number): number;

export function pv2x(p: number, v: number): number;

export function px(p: number, x: number, o_id: number): number;

export function px2h(p: number, x: number): number;

export function px2s(p: number, x: number): number;

export function px2t(p: number, x: number): number;

export function px2v(p: number, x: number): number;

export function sx(s: number, x: number, o_id: number): number;

export function sx2h(s: number, x: number): number;

export function sx2p(s: number, x: number): number;

export function sx2t(s: number, x: number): number;

export function sx2v(s: number, x: number): number;

export function th(t: number, h: number, o_id: number): number;

export function th2p(t: number, h: number): number;

export function th2s(t: number, h: number): number;

export function th2v(t: number, h: number): number;

export function th2x(t: number, h: number): number;

export function ts(t: number, s: number, o_id: number): number;

export function ts2h(t: number, s: number): number;

export function ts2p(t: number, s: number): number;

export function ts2v(t: number, s: number): number;

export function ts2x(t: number, s: number): number;

export function tv(t: number, v: number, o_id: number): number;

export function tv2h(t: number, v: number): number;

export function tv2p(t: number, v: number): number;

export function tv2s(t: number, v: number): number;

export function tv2x(t: number, v: number): number;

export function tx(t: number, x: number, o_id: number): number;

export function tx2h(t: number, x: number): number;

export function tx2p(t: number, x: number): number;

export function tx2s(t: number, x: number): number;

export function tx2v(t: number, x: number): number;

export type InitInput = RequestInfo | URL | Response | BufferSource | WebAssembly.Module;

export interface InitOutput {
    readonly memory: WebAssembly.Memory;
    readonly sx: (a: number, b: number, c: number) => number;
    readonly pt: (a: number, b: number, c: number) => number;
    readonly px2t: (a: number, b: number) => number;
    readonly hs2p: (a: number, b: number) => number;
    readonly hs2t: (a: number, b: number) => number;
    readonly hs2v: (a: number, b: number) => number;
    readonly ph2s: (a: number, b: number) => number;
    readonly ph2t: (a: number, b: number) => number;
    readonly ph2v: (a: number, b: number) => number;
    readonly ps2h: (a: number, b: number) => number;
    readonly ps2t: (a: number, b: number) => number;
    readonly ps2v: (a: number, b: number) => number;
    readonly pt2h: (a: number, b: number) => number;
    readonly pt2s: (a: number, b: number) => number;
    readonly pt2v: (a: number, b: number) => number;
    readonly pv2h: (a: number, b: number) => number;
    readonly pv2s: (a: number, b: number) => number;
    readonly pv2t: (a: number, b: number) => number;
    readonly sx2h: (a: number, b: number) => number;
    readonly sx2p: (a: number, b: number) => number;
    readonly sx2v: (a: number, b: number) => number;
    readonly th2p: (a: number, b: number) => number;
    readonly th2s: (a: number, b: number) => number;
    readonly th2v: (a: number, b: number) => number;
    readonly ts2h: (a: number, b: number) => number;
    readonly ts2p: (a: number, b: number) => number;
    readonly ts2v: (a: number, b: number) => number;
    readonly tv2h: (a: number, b: number) => number;
    readonly tv2p: (a: number, b: number) => number;
    readonly tv2s: (a: number, b: number) => number;
    readonly px2h: (a: number, b: number) => number;
    readonly px2s: (a: number, b: number) => number;
    readonly px2v: (a: number, b: number) => number;
    readonly sx2t: (a: number, b: number) => number;
    readonly tx2h: (a: number, b: number) => number;
    readonly tx2s: (a: number, b: number) => number;
    readonly tx2v: (a: number, b: number) => number;
    readonly pt2x: (a: number, b: number) => number;
    readonly tx: (a: number, b: number, c: number) => number;
    readonly hs2x: (a: number, b: number) => number;
    readonly ph2x: (a: number, b: number) => number;
    readonly ps2x: (a: number, b: number) => number;
    readonly pv2x: (a: number, b: number) => number;
    readonly th2x: (a: number, b: number) => number;
    readonly ts2x: (a: number, b: number) => number;
    readonly tv2x: (a: number, b: number) => number;
    readonly hx2p: (a: number, b: number) => number;
    readonly hx: (a: number, b: number, c: number) => number;
    readonly px: (a: number, b: number, c: number) => number;
    readonly hx2t: (a: number, b: number) => number;
    readonly hs: (a: number, b: number, c: number) => number;
    readonly ph: (a: number, b: number, c: number) => number;
    readonly ps: (a: number, b: number, c: number) => number;
    readonly pv: (a: number, b: number, c: number) => number;
    readonly th: (a: number, b: number, c: number) => number;
    readonly ts: (a: number, b: number, c: number) => number;
    readonly tv: (a: number, b: number, c: number) => number;
    readonly hx2s: (a: number, b: number) => number;
    readonly hx2v: (a: number, b: number) => number;
    readonly tx2p: (a: number, b: number) => number;
    readonly __wbindgen_externrefs: WebAssembly.Table;
    readonly __wbindgen_start: () => void;
}

export type SyncInitInput = BufferSource | WebAssembly.Module;

/**
 * Instantiates the given `module`, which can either be bytes or
 * a precompiled `WebAssembly.Module`.
 *
 * @param {{ module: SyncInitInput }} module - Passing `SyncInitInput` directly is deprecated.
 *
 * @returns {InitOutput}
 */
export function initSync(module: { module: SyncInitInput } | SyncInitInput): InitOutput;

/**
 * If `module_or_path` is {RequestInfo} or {URL}, makes a request and
 * for everything else, calls `WebAssembly.instantiate` directly.
 *
 * @param {{ module_or_path: InitInput | Promise<InitInput> }} module_or_path - Passing `InitInput` directly is deprecated.
 *
 * @returns {Promise<InitOutput>}
 */
export default function __wbg_init (module_or_path?: { module_or_path: InitInput | Promise<InitInput> } | InitInput | Promise<InitInput>): Promise<InitOutput>;
