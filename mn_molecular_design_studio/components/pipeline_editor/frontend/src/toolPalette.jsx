import React, {useMemo, useState} from 'react';

const artifactNames = {
  target: 'Target structure', molecule_task: 'Molecular inputs', candidates: 'Candidates',
  backbone: 'Backbones', sequence: 'Sequences', prediction: 'Predicted structures',
  compound_set: 'Compound set', ligand_sdf: 'Ligand SDF', docking_result: 'Docking results',
  fingerprint_manifest: 'Fingerprints', centroid_manifest: 'Centroids',
  prepared_complex: 'Prepared complex', md_trajectory: 'Trajectory + topology',
  md_analysis: 'MD measurements', analysis: 'Analysis results', tool_result: 'Saved results',
  active_learning_state: 'Screening state',
};
export const artifactLabel = type => artifactNames[type] || type.replaceAll('_', ' ');
export const portsFor = (spec, side) => spec?.capabilities?.[side] || [];

export function compatibleInlets(catalog, wires, sourceKind, output, types) {
  const options = {};
  const destinations = wires[sourceKind]?.[output] || {};
  for (const [kind, inlets] of Object.entries(destinations)) {
    const spec = catalog[kind];
    if (!spec || spec.hidden) continue;
    const matches = portsFor(spec, 'inputs').filter(port =>
      inlets[port.id] && port.artifact_types.some(type => types.includes(type)));
    if (matches.length) options[kind] = matches;
  }
  return options;
}

const unique = items => [...new Set(items)].sort((a, b) => a.localeCompare(b));
const modalityLabel = value => value === 'Ligand' ? 'Small molecule' : value;
const artifactsFor = (spec, side) => unique(portsFor(spec, side).flatMap(port => port.artifact_types));
const filterFields = [
  ['category', 'Category', spec => [spec.category], value => value],
  ['modality', 'Modality', spec => spec.palette?.modalities || [], modalityLabel],
  ['function', 'Function', spec => spec.palette?.functions || [], value => value],
  ['input', 'Input artifact', spec => artifactsFor(spec, 'inputs'), artifactLabel],
  ['output', 'Output artifact', spec => artifactsFor(spec, 'outputs'), artifactLabel],
];
const groupOrder = ['Input', 'Target inputs', 'Cofolding input', 'Library screening', 'Docking', 'Ligand library',
  'Structure prediction', 'Molecular dynamics', 'Analysis', 'Collection', 'Results'];
const groupRank = group => groupOrder.includes(group) ? groupOrder.indexOf(group) : groupOrder.length;

function ToolCard({spec, inlets, onAdd}) {
  const inputs = portsFor(spec, 'inputs');
  const outputs = portsFor(spec, 'outputs');
  const description = spec.palette?.description || spec.description;
  const preview = ports => ports.map(port => port.label).join(' · ');
  const add = () => onAdd(spec.key, inlets?.[0]?.id);
  const drag = event => event.dataTransfer.setData('application/studio-node', spec.key);
  const content = <>
    <span className="palette-card-header"><strong>{spec.label}</strong>
      {(!inlets || inlets.length === 1) && <span className="palette-add" aria-hidden="true">+</span>}
    </span>
    {description && <span className="palette-card-description" title={description}>{description}</span>}
    <span className="palette-card-ports">
      {inputs.length > 0 && <span><b>IN</b><span>{preview(inputs)}</span></span>}
      {outputs.length > 0 && <span><b>OUT</b><span>{preview(outputs)}</span></span>}
    </span>
  </>;
  return <article className="palette-card" style={{'--tool-color':spec.color}}>
    {inlets?.length > 1 ? <div className="palette-card-content">{content}</div> :
      <button className="palette-card-content" aria-label={'Add ' + spec.label}
        draggable={!inlets} onDragStart={drag} onClick={add}>{content}</button>}
    {inlets?.length > 0 && <div className="palette-inlets">
      {inlets.length === 1 ? <small>Connects to {inlets[0].label}</small> : inlets.map(inlet =>
        <button key={inlet.id} onClick={() => onAdd(spec.key, inlet.id)}>Connect to {inlet.label}</button>)}
    </div>}
  </article>;
}

export function ToolPalette({catalog, branch, usedStages, onConnect, onBranchPort, onClearBranch, onAdd, onClose, onDesignBinder}) {
  const [query, setQuery] = useState('');
  const [filters, setFilters] = useState({});
  const [advanced, setAdvanced] = useState(false);
  const [openGroups, setOpenGroups] = useState({});
  const available = useMemo(() => Object.values(catalog).filter(
    spec => !spec.hidden && (!branch || branch.inlets[spec.key])), [catalog, branch?.inlets]);
  const words = query.trim().toLowerCase().split(/\s+/).filter(Boolean);
  const searched = available.filter(spec => {
    const haystack = [spec.label, spec.key, spec.description, spec.palette?.description, spec.category,
      ...(spec.palette?.functions || []), ...(spec.palette?.modalities || []).map(modalityLabel),
      ...portsFor(spec, 'inputs').map(port => port.label),
      ...portsFor(spec, 'outputs').map(port => port.label)].filter(Boolean).join(' ').toLowerCase();
    return words.every(word => haystack.includes(word));
  });
  const matches = (spec, except) => filterFields.every(([key, , get]) =>
    key === except || !filters[key] || get(spec).includes(filters[key]));
  const matching = searched.filter(spec => matches(spec));
  const shown = matching.filter(spec => !usedStages.has(spec.key));
  const existing = (branch?.existing || []).filter(option => matching.some(spec => spec.key === option.kind));
  const grouped = new Map();
  for (const spec of shown) {
    if (!grouped.has(spec.category)) grouped.set(spec.category, []);
    grouped.get(spec.category).push(spec);
  }
  const groups = [...grouped.keys()].sort((a, b) => groupRank(a) - groupRank(b) || a.localeCompare(b));
  const filterCount = Object.values(filters).filter(Boolean).length;
  const revealGroups = Boolean(words.length || filterCount || branch);
  return <aside className="node-library tool-palette" aria-label="Tool library">
    <div className="palette-controls">
      <div className="panel-heading">
        <span>Tool library</span><button onClick={onClose} title="Hide library" aria-label="Hide library">‹</button>
      </div>
      {branch && <div className="palette-branch">
        <div className="palette-branch-heading"><strong>{branch.node.data.label}</strong>
          <button type="button" onClick={() => {setFilters({});onClearBranch();}}>Show all tools</button></div>
        {onDesignBinder && branch.node.data.kind === 'target' &&
          <button type="button" className="palette-design-binder" onClick={() => onDesignBinder(branch.node.id)}>
            Design a binder…
          </button>}
        <label>Branch from output
          <select aria-label="Branch from output" value={branch.port} onChange={event => {setFilters({});onBranchPort(event.target.value);}}>
            {portsFor(catalog[branch.node.data.kind], 'outputs').map(port =>
              <option key={port.id} value={port.id}>{port.label}</option>)}
          </select>
        </label>
      </div>}
      <input aria-label="Search pipeline tools" placeholder="Search tools…" value={query}
        onChange={event => setQuery(event.target.value)}/>
      <div className="palette-filter-bar">
        <button type="button" aria-expanded={advanced} onClick={() => setAdvanced(!advanced)}>
          Filters{filterCount ? ' (' + filterCount + ')' : ''}
        </button>
        {(filterCount > 0 || words.length > 0) && <button type="button"
          onClick={() => {setQuery('');setFilters({});}}>Clear</button>}
        <span role="status">{shown.length} {shown.length === 1 ? 'tool' : 'tools'}
          {existing.length > 0 && ' · ' + existing.length + ' existing'}</span>
      </div>
      {advanced && <div className="palette-filters">
        {filterFields.map(([key, label, get, format]) => {
          const values = unique(searched.filter(spec => matches(spec, key)).flatMap(get));
          const selected = filters[key] || '';
          if (selected && !values.includes(selected)) values.unshift(selected);
          return <label key={key}>{label}
            <select aria-label={label} value={selected} onChange={event => setFilters({...filters, [key]:event.target.value})}>
              <option value="">All</option>
              {values.map(value => <option key={value} value={value}>{format(value)}</option>)}
            </select>
          </label>;
        })}
      </div>}
    </div>
    {existing.length > 0 && <section className="palette-existing">
      <h3>Connect to an existing step</h3>
      {existing.map(({wire, label, port}) => <button key={wire.target + '/' + wire.targetHandle}
        onClick={() => onConnect(wire)}>Connect to {label} · {port}</button>)}
    </section>}
    {!groups.length && !existing.length && <p className="palette-empty">{branch && !available.length
      ? 'No registered tools accept this output.' : 'No tools match these filters.'}</p>}
    {groups.map(group => {
      const opened = openGroups[group] ?? (revealGroups || groupRank(group) < 3);
      return <section key={group}>
        <button className="palette-group" aria-expanded={opened}
          onClick={() => setOpenGroups({...openGroups, [group]:!opened})}>
          <span>{opened ? '▾' : '▸'} {group}</span><span>{grouped.get(group).length}</span>
        </button>
        {opened && grouped.get(group).sort((a, b) => a.label.localeCompare(b.label)).map(spec =>
          <ToolCard key={spec.key} spec={spec} inlets={branch?.inlets[spec.key]} onAdd={onAdd}/>)}
      </section>;
    })}
  </aside>;
}
