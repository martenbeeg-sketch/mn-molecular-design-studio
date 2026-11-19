const assert = require('node:assert/strict');
const {test} = require('node:test');
const {readFileSync} = require('node:fs');
const {resolve} = require('node:path');
const {JSDOM} = require('jsdom');

const assets = resolve(__dirname, '../../mn_molecular_design_studio/components/binderscope_inspector');
const source = readFileSync(resolve(assets, 'inspector.js'), 'utf8');
const template = readFileSync(resolve(assets, 'index.html'), 'utf8');

function fixture(confirmed = true) {
  const residue = (side, number) => ({id:`${side}:${number}`, side, number, aa:'A', label:`A${number + 300}`,
    author:`${side}:${number + 300}`, xyz:[number, 0, 0], confidence:side === 'B' && confirmed ? 80 : null,
    burial:20, ss:'H', interface:true});
  return {title:'Saved design', target_label:'Saved target', source_sha256:'verified', limitation:'One conformation',
    target_caption:'Saved numbering', hbond_basis:'computed', extra_chains:0, unmapped_interactions:0,
    confidence_confirmed:confirmed, residues:{A:[residue('A',1),residue('A',2)], B:[residue('B',1),residue('B',2)]},
    contacts:[{target:'A:1',binder:'B:1',kind:'H-bond + salt bridge',distance:3,pairs:3},
      {target:'A:1',binder:'B:2',kind:'Contact',distance:4,pairs:1}],
    interactions:[{kind:'H-bond',start:[0,0,0],end:[3,0,0]}, {kind:'Salt bridge',start:[0,0,0],end:[0,3,0]}],
    hotspots:[{id:'A:2',requested:'X302',status:'missed',minimum_distance_angstrom:6.84}],
    hotspot_status:'assessed', polar:[{id:'B:2',atom:'N'}],
    patches:[{area:30,res:[1]}], liabilities:[{kind:'Test motif',pos:2,reason:'Recorded warning',sev:'warn'}],
    reference:{status:'unavailable',reason:'No fixed original reference'}, pdb:'display-copy'};
}

function inspector(payload = fixture()) {
  const dom = new JSDOM(template, {runScripts:'outside-only'});
  const {window} = dom;
  window.document.getElementById('assessment-data').textContent = JSON.stringify(payload);
  const calls = {cylinders:[], surfaces:[], zooms:[], styles:[], renders:0};
  const model = {setStyle:(selection,style)=>calls.styles.push({selection,style})};
  const viewer = {
    addModel:()=>model, selectedAtoms:()=>[{x:0,y:0,z:0}],
    removeAllShapes:()=>{calls.cylinders=[];}, removeAllLabels:()=>{}, removeAllSurfaces:()=>{},
    removeShape:()=>{}, removeLabel:()=>{}, addSphere:()=>({}), addLabel:()=>({}),
    addCylinder:shape=>{calls.cylinders.push(shape);return shape;},
    addSurface:(type,style,selection)=>{calls.surfaces.push({type,style,selection});return Promise.resolve();},
    render:()=>{calls.renders++;}, zoomTo:selection=>calls.zooms.push(selection), zoom:()=>{},
    setHoverable:()=>{}, setClickable:()=>{}, resize:()=>{}, pngURI:()=>''
  };
  window.$3Dmol = {SurfaceType:{VDW:1},createViewer:()=>viewer};
  window.ResizeObserver = class {observe() {}};
  window.eval(source);
  const document = window.document;
  return {window,document,calls,close:()=>window.close()};
}

test('all analysis views switch locally without changing the measurements', () => {
  const view = inspector();
  for (const name of ['findings','confidence','footprint','interface']) {
    view.document.getElementById('tab-'+name).click();
    assert.equal(view.document.getElementById('tab-'+name).getAttribute('aria-selected'), 'true');
    assert.equal(view.document.getElementById('panel-'+name).hidden, false);
    assert.equal([...view.document.querySelectorAll('[role=tabpanel]')].filter(panel=>!panel.hidden).length, 1);
  }
  assert.match(view.document.getElementById('findings').textContent, /1 intended hotspot not contacted/);
  assert.equal(view.document.querySelectorAll('.network-edge').length, 2);
  assert.equal(view.document.getElementById('viewer').dataset.ready, 'true');
  view.close();
});

test('residue hover and keyboard selection link sequence, network and exact 3D selection', () => {
  const view = inspector();
  const cell = view.document.querySelector('.track-cell[data-residue="B:1"]');
  cell.dispatchEvent(new view.window.Event('pointerenter'));
  assert.match(view.document.getElementById('selection').textContent, /saved B:301/);
  assert(view.document.querySelector('.network-edge[data-binder="B:1"]').classList.contains('picked'));
  cell.click();
  assert.deepEqual(JSON.parse(JSON.stringify(view.calls.zooms.at(-1))), {or:[{chain:'B',resi:1}]});
  const pair = view.document.querySelector('g.edge');
  pair.dispatchEvent(new view.window.KeyboardEvent('keydown', {key:'Enter',bubbles:true}));
  assert.match(view.document.getElementById('selection').textContent, /H-bond \+ salt bridge/);
  assert.deepEqual(JSON.parse(JSON.stringify(view.calls.zooms.at(-1))), {or:[{chain:'A',resi:1},{chain:'B',resi:1}]});
  assert.match(view.document.getElementById('selection').textContent, /not a chemical bond/);
  view.close();
});

test('overlays use supplied atom geometry and proteins stay cartoons', async () => {
  const view = inspector();
  assert.equal(view.calls.cylinders.length, 2);
  assert.deepEqual(JSON.parse(JSON.stringify(view.calls.cylinders[0].end)), {x:3,y:0,z:0});
  const checkbox = view.document.getElementById('show-hbonds');
  checkbox.checked = false;
  checkbox.dispatchEvent(new view.window.Event('change'));
  assert.equal(view.calls.cylinders.length, 1);
  view.document.getElementById('show-footprint').click();
  await Promise.resolve();
  assert.equal(view.document.getElementById('show-surface').checked, true);
  assert.equal(view.calls.surfaces.length, 1);
  assert.deepEqual(JSON.parse(JSON.stringify(view.calls.surfaces[0].selection)), {chain:'A'});
  assert.match(view.calls.surfaces[0].style.colorfunc({chain:'A',resi:1}), /^#[0-9a-f]{6}$/);
  assert(view.calls.styles.every(call=>call.style.cartoon && !call.style.stick));
  view.close();
});

test('unknown confidence is not inferred; confirmed confidence applies only to binder', () => {
  const unknown = inspector(fixture(false));
  assert.equal(unknown.document.querySelector('[data-color="confidence"]').disabled, true);
  assert.match(unknown.document.getElementById('confidence-panel').textContent, /No confirmed per-residue binder pLDDT/);
  unknown.close();
  const confirmed = inspector();
  confirmed.document.querySelector('[data-color="confidence"]').click();
  const color = confirmed.calls.styles.at(-1).style.cartoon.colorfunc;
  assert.equal(color({chain:'A',resi:1}), '#477aa9');
  assert.equal(color({chain:'B',resi:1}), '#70b9d1');
  confirmed.close();
});

test('untrusted names remain text, not markup or script', () => {
  const payload = fixture();
  payload.title = '<img src=x onerror="bad()">';
  payload.target_label = '<script>bad()</script>';
  const view = inspector(payload);
  assert.equal(view.document.getElementById('design-title').textContent, payload.title);
  assert.equal(view.document.getElementById('design-title').children.length, 0);
  assert.equal(view.document.getElementById('target-label').children.length, 0);
  view.close();
});
