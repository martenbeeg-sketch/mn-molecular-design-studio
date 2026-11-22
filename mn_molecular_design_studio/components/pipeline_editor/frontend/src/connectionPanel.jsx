import React, {useState} from 'react';

const normalizeType = type => type.replaceAll('_', ' ');
const labels = ports => ports.length ? ports.map(item => (
  <li key={item.id}><strong>{item.label}</strong><small>{(item.artifact_types || []).map(normalizeType).join(', ')}</small></li>
)) : <li className="panel-hint">None</li>;

export function ConnectionPanel({node, catalog, matrix, close, configure, connect, getTypes, readOnly}) {
  const [port, setPort] = useState(null);
  const [search, setSearch] = useState('');
  const [handoffsOpen, setHandoffsOpen] = useState(false);
  const spec = node && catalog[node.data.kind];
  const inputs = spec?.capabilities?.inputs || [];
  const outputs = spec?.capabilities?.outputs || [];
  const chosen = outputs.some(row => row.id === port) ? port : outputs[0]?.id;
  const outgoing = chosen ? matrix.wires[node.data.kind]?.[chosen] || {} : {};
  const types = chosen ? getTypes(node, chosen) : [];
  const options = Object.entries(outgoing).flatMap(([kind, matching]) => {
    const target = catalog[kind];
    if (!target || target.hidden) return [];
    return Object.entries(matching).flatMap(([inlet, artifact]) => {
      const input = target.capabilities?.inputs.find(item => item.id === inlet);
      if (!input || !types.some(type => input.artifact_types.includes(type))) return [];
      return [{kind, inlet, artifact: types.find(type => input.artifact_types.includes(type)) || artifact,
        label: target.label, input: input.label}];
    });
  }).filter(item => `${item.label} ${item.input} ${item.artifact}`.toLowerCase().includes(search.trim().toLowerCase()))
    .sort((a, b) => a.label.localeCompare(b.label) || a.input.localeCompare(b.input));
  const handoffs = spec ? [...new Set([spec.workflow, ...(spec.workflow_aliases || [])])]
    .flatMap(workflow => matrix.handoffs[workflow] || []) : [];
  return <aside className="node-inspector connection-panel" aria-label="Connections">
    <div className="panel-heading"><span>Connections</span><button onClick={close} aria-label="Close connections">×</button></div>
    {!spec ? <p className="panel-hint">Select a step to see its inputs, outputs and compatible tools.</p> : <>
      <div className="port-summary">
        <h2>{node.data.label}</h2>
        <button onClick={configure}>Configure step</button>
        <dl><dt>Inputs</dt><dd><ul>{labels(inputs)}</ul></dd><dt>Outputs</dt><dd><ul>{labels(outputs)}</ul></dd></dl>
      </div>
      {readOnly ? <p className="panel-hint">The screening recipe connects these stages.</p> : chosen && <section className="connection-options">
        <h3>Add next step</h3>
        <label>Output<select aria-label="Output socket" value={chosen} onChange={event => setPort(event.target.value)}>
          {outputs.map(item => <option key={item.id} value={item.id}>{item.label}</option>)}
        </select></label>
        <input aria-label="Find a compatible step" value={search} onChange={event => setSearch(event.target.value)} placeholder="Find a compatible step"/>
        {options.length ? <ul>{options.map(item => <li key={`${item.kind}/${item.inlet}`}>
          <button onClick={() => connect(chosen, item.kind, item.inlet)}>
            <span><strong>{item.label}</strong><small>{item.input} · {normalizeType(item.artifact)}</small></span><span aria-hidden="true">+</span>
          </button>
        </li>)}</ul> : <p className="panel-hint">No matching steps for this output.</p>}
      </section>}
      {!!handoffs.length && <section className="standalone-handoffs">
        <button aria-expanded={handoffsOpen} onClick={() => setHandoffsOpen(!handoffsOpen)}>Saved-result handoffs</button>
        {handoffsOpen && <ul>{handoffs.map(item => <li key={`${item.source}/${item.target}`}>
          <strong>{Object.values(catalog).find(row => row.workflow === item.target)?.label || item.target}</strong>
          <small>{normalizeType(item.artifact_type)} via {item.transfer === 'source_job' ? 'saved job' : 'manifest'}</small>
          {!!item.additional_inputs.length && <small>Also needs {item.additional_inputs.map(normalizeType).join(', ')}</small>}
        </li>)}</ul>}
      </section>}
    </>}
  </aside>;
}
