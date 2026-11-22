import React, {useEffect} from 'react';
import {Handle, Position, useUpdateNodeInternals} from '@xyflow/react';

export function CollapsedGroup({id, data, selected}) {
  const update = useUpdateNodeInternals();
  const socketLayout = JSON.stringify(data.boundary);
  useEffect(() => update(id), [id, socketLayout, update]);
  return <div className={`studio-node group-card ${selected ? 'is-selected' : ''}`}>
    <div className="node-category">Group · {data.members.length} steps</div>
    <div className="node-title">{data.name}</div>
    <div className="ports">
      {data.boundary.in.map(port => <div key={port.handle} className="port-label input">
        <Handle id={port.handle} type="target" position={Position.Left}/><span>{port.label}</span>
      </div>)}
      {data.boundary.out.map(port => <div key={port.handle} className="port-label output">
        <span>{port.label}</span><Handle id={port.handle} type="source" position={Position.Right}/>
      </div>)}
    </div>
  </div>;
}

export function Frame({data, selected}) {
  return <div className={`canvas-frame ${data.loop ? 'repeat-frame' : ''} ${selected ? 'frame-selected' : ''}`}>
    <div className="frame-caption"><strong>{data.name}</strong><span>{data.subtitle || `${data.members.length} steps`}</span></div>
  </div>;
}

export function RoundOutline({program, configure, close}) {
  if (!program?.screen) return null;
  const {screen, blocks} = program;
  const step = id => <button key={id} className={`round-step ${blocks[id].role}`} onClick={() => configure(id)}>
    <strong>{blocks[id].label}</strong>
    {screen.step_details[id] && <small>{screen.step_details[id]}</small>}
    {(screen.source_refs[id] || []).map(([label, source]) => <small key={source}>{label} <b>{blocks[source].label}</b></small>)}
    {id === screen.selection_id && <>
      <small>{screen.selection.example}</small>
      <small>{screen.selection.predicted_count} model-ranked · {screen.selection.random_count} random</small>
    </>}
  </button>;
  return <aside className="node-inspector round-inspector" aria-label="Screening round">
    <div className="panel-heading"><span>Screening round</span><button onClick={close} aria-label="Close screening round">×</button></div>
    <div className="round-outline">
      <div className="round-stack">{screen.setup.map(step)}</div>
      <section className="round-scope">
        <h3>Repeat docking rounds</h3>
        <div className="round-stack">{screen.rounds.map(step)}</div>
        <section className="round-scope round-decision">
          <button className="round-condition" onClick={() => configure(screen.decision)}>
            <strong>If any stopping condition is met</strong>
            {screen.stop_conditions.map(condition => <small key={condition}>{condition}</small>)}
          </button>
          <div className="round-stack"><div className="round-exit">Exit loop</div></div>
          <h4>Otherwise</h4>
          <div className="round-stack">{screen.next_round.map(step)}</div>
        </section>
        <div className="round-return">↻ Dock the next batch</div>
      </section>
      <div className="round-stack">{screen.finish.map(step)}</div>
    </div>
  </aside>;
}
