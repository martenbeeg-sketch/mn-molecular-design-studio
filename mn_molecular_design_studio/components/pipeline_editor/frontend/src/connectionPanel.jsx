import React, {useState} from 'react';
import {artifactLabel, portsFor} from './toolPalette.jsx';

const labels = (ports, metadata) => ports.length ? ports.map(item => (
  <li key={item.id}><strong>{item.label}</strong>
    <small>{item.artifact_types.map(artifactLabel).join(', ')}</small>
    {metadata && <small>{item.required ? 'Required connection' : 'Connect or configure in the step'}
      {item.multiple && ' · accepts multiple links'}</small>}
  </li>
)) : <li className="panel-hint">None</li>;

export function ConnectionPanel({node, catalog, matrix, close, configure, branch}) {
  const [handoffsOpen, setHandoffsOpen] = useState(false);
  const spec = node && catalog[node.data.kind];
  const inputs = portsFor(spec, 'inputs');
  const outputs = portsFor(spec, 'outputs');
  const handoffs = spec ? [...new Set([spec.workflow, ...(spec.workflow_aliases || [])])]
    .flatMap(workflow => matrix.handoffs[workflow] || []) : [];
  return <aside className="node-inspector connection-panel" aria-label="Connections">
    <div className="panel-heading"><span>Connections</span><button onClick={close} aria-label="Close connections">×</button></div>
    {!spec ? <p className="panel-hint">Select a step to see its inputs, outputs and compatible tools.</p> : <>
      <div className="port-summary">
        <h2>{node.data.label}</h2>
        <button onClick={configure}>Configure step</button>
        <dl><dt>Inputs</dt><dd><ul>{labels(inputs, true)}</ul></dd><dt>Outputs</dt><dd><ul>{labels(outputs)}</ul></dd></dl>
      </div>
      {outputs.length > 0 && <section className="connection-options">
        <h3>Add next step</h3>
        <ul>{outputs.map(port => <li key={port.id}>
          <button onClick={() => branch(port.id)}>Add from {port.label}</button>
        </li>)}</ul>
      </section>}
      {!!handoffs.length && <section className="standalone-handoffs">
        <button aria-expanded={handoffsOpen} onClick={() => setHandoffsOpen(!handoffsOpen)}>Saved-result handoffs</button>
        {handoffsOpen && <ul>{handoffs.map(item => <li key={`${item.source}/${item.target}`}>
          <strong>{Object.values(catalog).find(row => row.workflow === item.target)?.label || item.target}</strong>
          <small>{artifactLabel(item.artifact_type)} via {item.transfer === 'source_job' ? 'saved job' : 'manifest'}</small>
          {!!item.additional_inputs.length && <small>Also needs {item.additional_inputs.map(artifactLabel).join(', ')}</small>}
        </li>)}</ul>}
      </section>}
    </>}
  </aside>;
}
