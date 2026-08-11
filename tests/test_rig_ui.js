/* test_rig_ui.js -- headless smoke test of the rig page's UI layer.
 * The model tests cannot catch UI-layer crashes (a thrown exception in a
 * render freezes the page while looking merely idle -- exactly the
 * shipped marks-vs-atoms indexing bug this guards against). A minimal
 * DOM shim loads rigapp.js, drags the box, banks poses, and walks both
 * phases, asserting no handler throws and the visible counters move. */
'use strict';

var failures = 0;
function ok(cond, name) {
  console.log((cond ? 'ok:   ' : 'FAIL: ') + name);
  if (!cond) failures++;
}

/* ---- minimal DOM shim -------------------------------------------------- */
function makeCtx() {
  return new Proxy({}, {
    get: function (t, k) {
      if (k === 'measureText') return function () { return { width: 0 }; };
      if (k === 'getImageData')
        return function () { return { data: new Uint8ClampedArray(4 * 160 * 120) }; };
      if (k === 'createImageData')
        return function (w, h) { return { data: new Uint8ClampedArray(4 * w * h) }; };
      if (typeof k === 'string' && k !== 'canvas') return function () {};
      return undefined;
    },
    set: function () { return true; }
  });
}
var elements = {};
function el(id) {
  if (!elements[id]) elements[id] = {
    id: id, width: 320, height: 240, value: '16', checked: false,
    textContent: '', disabled: false, listeners: {},
    classList: { add: function () {}, remove: function () {} },
    style: {},
    getContext: function () { return makeCtx(); },
    addEventListener: function (ev, fn) {
      (this.listeners[ev] = this.listeners[ev] || []).push(fn);
    },
    getBoundingClientRect: function () { return { left: 0, top: 0 }; }
  };
  return elements[id];
}
var windowListeners = {};
global.document = {
  getElementById: function (id) { return el(id); },
  createElement: function () {
    return { width: 0, height: 0, getContext: function () { return makeCtx(); } };
  }
};
global.window = {
  addEventListener: function (ev, fn) {
    (windowListeners[ev] = windowListeners[ev] || []).push(fn);
  },
  location: { search: '', href: 'test' }
};
function fire(target, ev, evt) {
  ((target.listeners && target.listeners[ev]) || []).forEach(function (fn) {
    fn.call(target, evt);
  });
}
function fireWindow(ev, evt) {
  (windowListeners[ev] || []).forEach(function (fn) { fn(evt); });
}
function dragBox(dx, dy, mods) {
  var evt = { clientX: 100, clientY: 100, buttons: 1,
              shiftKey: !!(mods && mods.shift), altKey: !!(mods && mods.alt),
              preventDefault: function () {} };
  fire(el('overview'), 'mousedown', evt);
  fireWindow('mousemove', { clientX: 100 + dx, clientY: 100 + dy, buttons: 1,
    shiftKey: evt.shiftKey, altKey: evt.altKey, preventDefault: function () {} });
  fireWindow('mouseup', {});
}

/* ---- load the page ----------------------------------------------------- */
var fs = require('fs');
global.MV = require('../web/model.js');
var threw = null;
try {
  eval(fs.readFileSync(__dirname + '/../web/rigapp.js', 'utf8'));
} catch (e) { threw = e; }
ok(!threw, 'rigapp loads and renders' + (threw ? ' (' + threw.message + ')' : ''));

/* ---- phase 1: manual bank, then auto-banks via drags ------------------- */
try {
  fire(el('bank'), 'click', {});
  ok(/1 pose/.test(el('bankcount').textContent),
     'manual bank increments the visible counter (' +
     el('bankcount').textContent + ')');
  var moves = [[60, 0], [0, 45], [-40, -30, { shift: 1 }], [30, 25, { shift: 1 }],
               [50, -20], [-60, 10, { shift: 1 }], [20, 40], [45, 15, { shift: 1 }]];
  moves.forEach(function (m) { dragBox(m[0], m[1], m[2]); });
  var n = parseInt(el('bankcount').textContent, 10);
  ok(n >= 4, 'dragging auto-banks novel poses (' + n + ' banked)');
  ok(el('readout').textContent.indexOf('rig solve from') >= 0,
     'readout reports a solve');
  threw = null;
} catch (e) { threw = e; }
ok(!threw, 'phase 1 interactions never throw' +
   (threw ? ' (' + threw.message + ')' : ''));

/* ---- phase 2: accept, move the box, model accumulates ------------------ */
try {
  fire(el('tomeasure'), 'click', {});
  ok(/phase 2/.test(el('phasebadge').textContent),
     'accept switches to phase 2 (' + el('phasebadge').textContent + ')');
  dragBox(25, -20);
  dragBox(-30, 15, { shift: 1 });
  ok(el('readout').textContent.indexOf('model accumulation') >= 0,
     'measure phase reports cloud/TSDF accumulation');
  fire(el('cloudview'), 'click', {});
  ok(/cloud \d+ points/.test(el('readout').textContent),
     'clicking the cloud panel clears and re-estimates the current pose');
  fire(el('clearmodel'), 'click', {});
  fire(el('tosolve'), 'click', {});
  fire(el('resetsolve'), 'click', {});
  ok(/0 poses/.test(el('bankcount').textContent), 'reset returns counter to 0');
  threw = null;
} catch (e) { threw = e; }
ok(!threw, 'phase 2 interactions never throw' +
   (threw ? ' (' + threw.message + ')' : ''));

if (failures) {
  console.log(failures + ' rig-ui test(s) FAILED');
  process.exit(1);
}
console.log('all rig-ui tests passed');
