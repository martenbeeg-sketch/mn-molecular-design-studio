import React, {memo, useCallback, useEffect, useRef, useState} from 'react';
import {createRoot} from 'react-dom/client';
import {ReactFlow, ReactFlowProvider, Background, Controls, MiniMap, Handle, Position, applyNodeChanges, applyEdgeChanges, addEdge, useReactFlow, useNodesInitialized} from '@xyflow/react';
import '@xyflow/react/dist/style.css';
import './style.css';

const post=(type,data={})=>window.parent.postMessage({isStreamlitMessage:true,type,...data},'*');
const copy=x=>JSON.parse(JSON.stringify(x));
const uid=()=>crypto.randomUUID().replaceAll('-','');
let catalog={};
const authored=(nodes,edges,base,viewport)=>({...base,nodes:nodes.map(n=>({id:n.id,type:'studio',position:n.position,data:{kind:n.data.kind,label:n.data.label,params:n.data.params}})),edges:edges.map(e=>({id:e.id,source:e.source,target:e.target,sourceHandle:e.sourceHandle||'out',targetHandle:e.targetHandle||'in'})),viewport});
const statusLabels={completed:'Complete',cached:'Reused',running:'Running',queued:'Queued',failed:'Failed',blocked:'Blocked',pending:'Pending'};

const StudioNode=memo(({data,selected})=>{
  const spec=catalog[data.kind]||{}; const p=data.params||{};
  const summary=p.num_designs?`${p.num_designs} design${p.num_designs===1?'':'s'} · ${p.binder_length||'variable'} aa`:p.sequences?`${p.sequences} sequence${p.sequences===1?'':'s'} / backbone`:data.kind==='target'?(p.target_pdb?`${p.target_chains||'A'} · target loaded`:'Choose a target structure'):data.kind==='filter'?(p.metric?`${p.metric} ${p.operator} ${p.threshold}`:p.passing_only?'Filter-passing designs':'All candidates'):data.kind==='results'?'Overlay · confidence · RMSD':spec.category;
  return <div className={`studio-node ${selected?'is-selected':''}`} style={{'--node-color':spec.color}}>
    <div className="node-category"><span className="category-dot"/>{spec.category}{data.status&&<span className={`node-status ${data.status}`}>{statusLabels[data.status]||data.status}</span>}</div>
    <div className="node-title">{data.label||spec.label}</div>
    <div className="node-summary">{summary}</div><div className="node-edit-hint">Click to configure <span>↗</span></div>
    {data.count!==undefined&&<div className="node-count">{data.count} candidates</div>}
    <div className="ports">{spec.inputs?.map(p=><div key={p.id} className="port-label input"><Handle id={p.id} type="target" position={Position.Left}/><span>{p.label}</span></div>)}{spec.outputs?.map(p=><div key={p.id} className="port-label output"><span>{p.label}</span><Handle id={p.id} type="source" position={Position.Right}/></div>)}</div>
  </div>;
});
const nodeTypes={studio:StudioNode};

function Editor({args}){
  catalog=args.catalog;
  const initial=args.graph;
  const [nodes,setNodes]=useState(initial.nodes),[edges,setEdges]=useState(initial.edges);
  const [query,setQuery]=useState(''),[notice,setNotice]=useState(''),[selection,setSelection]=useState(null);
  const [palette,setPalette]=useState(true),[focus,setFocus]=useState(false),[historyTick,setHistoryTick]=useState(0),[helpOpen,setHelpOpen]=useState(false),[showMap,setShowMap]=useState(false);
  const base=useRef(initial),viewport=useRef(initial.viewport||{x:0,y:0,zoom:1}),reset=useRef(args.reset_token),lastSent=useRef('');
  const state=useRef({nodes,edges}),history=useRef({past:[],future:[]}),clipboard=useRef(null),dragSnapshot=useRef(null);
  const rf=useReactFlow(), initialized=useNodesInitialized(), fitPending=useRef(true); state.current={nodes,edges};
  const snapshot=useCallback(()=>authored(state.current.nodes,state.current.edges,base.current,viewport.current),[]);
  const emit=useCallback((action=null)=>{const graph=snapshot();lastSent.current=JSON.stringify(graph);post('streamlit:setComponentValue',{value:{graph,...(action?{action:{...(typeof action==='string'?{type:action}:action),id:uid()}}:{})}});},[snapshot]);
  const pushHistory=useCallback((snap=null)=>{const h=history.current;h.past.push(snap||snapshot());if(h.past.length>60)h.past.shift();h.future=[];setHistoryTick(x=>x+1);},[snapshot]);
  const restore=useCallback(g=>{base.current={...base.current,name:g.name};setNodes(copy(g.nodes));setEdges(copy(g.edges));if(g.viewport){viewport.current=g.viewport;rf.setViewport(g.viewport);}},[rf]);
  const undo=useCallback(()=>{const h=history.current;if(!h.past.length)return;h.future.push(snapshot());restore(h.past.pop());setHistoryTick(x=>x+1);},[snapshot,restore]);
  const redo=useCallback(()=>{const h=history.current;if(!h.future.length)return;h.past.push(snapshot());restore(h.future.pop());setHistoryTick(x=>x+1);},[snapshot,restore]);
  useEffect(()=>{
    if(reset.current!==args.reset_token){reset.current=args.reset_token;base.current=args.graph;setNodes(copy(args.graph.nodes));setEdges(copy(args.graph.edges));history.current={past:[],future:[]};setSelection(null);viewport.current=args.graph.viewport;fitPending.current=true;}
    base.current={...base.current,name:args.graph.name};
  },[args.reset_token,args.graph.name,rf]);
  const settingsApplied=useRef(args.settings_update?.id);
  useEffect(()=>{
    const change=args.settings_update;
    if(!change||settingsApplied.current===change.id)return;
    settingsApplied.current=change.id;
    const current=state.current.nodes.find(n=>n.id===change.node_id);
    if(!current||JSON.stringify(current.data.params)===JSON.stringify(change.params)&&current.data.label===change.label)return;
    pushHistory();
    setNodes(prev=>prev.map(n=>n.id===change.node_id?{...n,data:{...n.data,params:copy(change.params),label:change.label}}:n));
  },[args.settings_update?.id,pushHistory]);
  useEffect(()=>{
    if(!args.input_contexts)return;
    setNodes(prev=>prev.map(n=>args.input_contexts[n.id]?{...n,data:{...n.data,params:{...n.data.params,...args.input_contexts[n.id]}}}:n));
  },[JSON.stringify(args.input_contexts)]);
  useEffect(()=>{
    const timer=setTimeout(()=>{const g=snapshot(),serialized=JSON.stringify(g);if(serialized!==lastSent.current){lastSent.current=serialized;post('streamlit:setComponentValue',{value:{graph:g}});}},500);
    return()=>clearTimeout(timer);
  },[nodes,edges,snapshot]);
  useEffect(()=>{if(initialized&&fitPending.current){const id=requestAnimationFrame(()=>{rf.fitView({padding:.4,maxZoom:1});fitPending.current=false;});return()=>cancelAnimationFrame(id);}},[initialized,nodes,rf,args.reset_token]);
  useEffect(()=>{post('streamlit:setFrameHeight',{height:focus?950:760});},[focus]);
  useEffect(()=>{if(notice){const t=setTimeout(()=>setNotice(''),5500);return()=>clearTimeout(t);}},[notice]);
  const addNode=useCallback((kind,position)=>{
    const spec=catalog[kind];pushHistory();const pos=position||rf.screenToFlowPosition({x:window.innerWidth*.49,y:330});
    const offsets=[[0,0],[0,210],[320,0],[-320,0],[0,-210],[320,210],[-320,210],[0,420],[640,0]];const open=offsets.map(([x,y])=>({x:pos.x+x,y:pos.y+y})).find(p=>!state.current.nodes.some(n=>Math.abs(n.position.x-p.x)<280&&Math.abs(n.position.y-p.y)<170))||{x:pos.x,y:pos.y+630};
    const n={id:`${kind}-${uid().slice(0,10)}`,type:'studio',position:open,selected:true,data:{kind,label:spec.label,params:Object.fromEntries(spec.fields.map(f=>[f.name,f.default]))}};
    setNodes(prev=>[...prev.map(n=>({...n,selected:false})),n]);setSelection(n.id);setNotice(`${spec.label} added. Connect its sockets to another step.`);
  },[pushHistory,rf]);
  const compatible=useCallback(connection=>{
    const ns=state.current.nodes,es=state.current.edges;
    const source=ns.find(n=>n.id===connection.source),target=ns.find(n=>n.id===connection.target);
    if(!source||!target||source.id===target.id)return false;
    const a=catalog[source.data.kind],b=catalog[target.data.kind];
    if(!a.outputs.length||!b.inputs.length)return false;
    const infer=(n,seen=new Set())=>{if(!n||seen.has(n.id))return null;seen.add(n.id);if(n.data.kind==='filter'){const input=es.find(e=>e.target===n.id);return input?infer(ns.find(v=>v.id===input.source),seen):'candidates';}return catalog[n.data.kind].outputs[0]?.type;};
    const accepts=b.inputs[0].accepts,output=infer(source);
    if(output==='candidates'?accepts.includes('target'):!accepts.includes(output))return false;
    if(!b.inputs[0].multiple&&es.some(e=>e.target===target.id))return false;
    const visit=(id,seen=new Set())=>{if(id===source.id)return true;if(seen.has(id))return false;seen.add(id);return es.filter(e=>e.source===id).some(e=>visit(e.target,seen));};
    return !visit(target.id);
  },[]);
  const onConnect=useCallback(c=>{if(!compatible(c)){setNotice('These sockets are incompatible, already connected, or would create a loop.');return;}pushHistory();setEdges(prev=>addEdge({...c,id:`edge-${uid()}`},prev));},[compatible,pushHistory]);
  const duplicate=useCallback(()=>{
    const selected=state.current.nodes.filter(n=>n.selected||n.id===selection);if(!selected.length)return;
    pushHistory();const ids=Object.fromEntries(selected.map(n=>[n.id,`${n.data.kind}-${uid().slice(0,10)}`]));
    const clones=selected.map(n=>({...copy(n),id:ids[n.id],position:{x:n.position.x+48,y:n.position.y+48},selected:true}));
    const newEdges=state.current.edges.filter(e=>ids[e.source]&&ids[e.target]).map(e=>({...copy(e),id:`edge-${uid()}`,source:ids[e.source],target:ids[e.target]}));
    setNodes(prev=>[...prev.map(n=>({...n,selected:false})),...clones]);setEdges(prev=>[...prev,...newEdges]);setSelection(clones[0].id);
  },[pushHistory,selection]);
  useEffect(()=>{
    const handler=e=>{
      if(['INPUT','TEXTAREA','SELECT'].includes(e.target.tagName)||e.target.isContentEditable)return;
      const mod=e.ctrlKey||e.metaKey;
      if(mod&&e.key.toLowerCase()==='z'){e.preventDefault();e.shiftKey?redo():undo();}
      if(mod&&e.key.toLowerCase()==='y'){e.preventDefault();redo();}
      if(mod&&e.key.toLowerCase()==='d'){e.preventDefault();duplicate();}
      if(mod&&e.key.toLowerCase()==='c'){const ns=state.current.nodes.filter(n=>n.selected);clipboard.current={nodes:copy(ns),edges:copy(state.current.edges.filter(e=>ns.some(n=>n.id===e.source)&&ns.some(n=>n.id===e.target)))};}
      if(mod&&e.key.toLowerCase()==='v'&&clipboard.current?.nodes.length){e.preventDefault();pushHistory();const c=clipboard.current,ids=Object.fromEntries(c.nodes.map(n=>[n.id,`${n.data.kind}-${uid().slice(0,10)}`]));setNodes(p=>[...p.map(n=>({...n,selected:false})),...c.nodes.map(n=>({...copy(n),id:ids[n.id],position:{x:n.position.x+70,y:n.position.y+70},selected:true}))]);setEdges(p=>[...p,...c.edges.map(e=>({...e,id:`edge-${uid()}`,source:ids[e.source],target:ids[e.target]}))]);}
      if(e.key.toLowerCase()==='f'){e.preventDefault();rf.fitView({padding:.4,duration:250});}
      if(e.key==='Escape')setSelection(null);
    };window.addEventListener('keydown',handler);return()=>window.removeEventListener('keydown',handler);
  },[undo,redo,duplicate,pushHistory,rf]);
  const layout=()=>{pushHistory();const ns=state.current.nodes,es=state.current.edges;const levels=new Map(),visiting=new Set();const level=id=>{if(levels.has(id))return levels.get(id);if(visiting.has(id))return 0;visiting.add(id);const ins=es.filter(e=>e.target===id);const l=ins.length?Math.max(...ins.map(e=>level(e.source)))+1:0;levels.set(id,l);return l;};const rows={};setNodes(ns.map(n=>{const col=level(n.id),row=rows[col]||0;rows[col]=row+1;return {...n,position:{x:70+col*320,y:110+row*210}};}));setTimeout(()=>rf.fitView({padding:.4,duration:250}),30);};
  const displayNodes=nodes.map(n=>({...n,data:{...n.data,status:args.execution?.[n.id]?.status,count:args.execution?.[n.id]?.candidate_count}}));
  const categories=[...new Set(Object.values(catalog).map(s=>s.category))];
  return <div className={`editor-shell ${focus?'focus-mode':''}`}>
    <header className="editor-toolbar"><div className="toolbar-title"><span className="flow-symbol">◇</span><strong>Pipeline canvas</strong><span className="draft-chip">Draft</span></div><div className="toolbar-actions"><button title="Undo (Ctrl+Z)" disabled={!history.current.past.length} onClick={undo}>↶</button><button title="Redo (Ctrl+Shift+Z)" disabled={!history.current.future.length} onClick={redo}>↷</button><span className="toolbar-divider"/><button onClick={layout} title="Arrange the current graph">Arrange</button><button title="Toggle canvas overview" onClick={()=>setShowMap(!showMap)}>Map</button><button title="Canvas shortcuts and help" onClick={()=>{setHelpOpen(!helpOpen);setSelection(null);}}>?</button><button onClick={()=>{setFocus(!focus);setTimeout(()=>rf.fitView({padding:.4,duration:200}),80);}}>{focus?'Exit focus':'Focus canvas'}</button><button onClick={()=>emit('save_version')}>Save version</button><button className="run-button" onClick={()=>emit('run')}>▶ Run pipeline</button></div></header>
    <div className="workspace">
      {palette&&!focus&&<aside className="node-library"><div className="panel-heading"><span>Add a step</span><button onClick={()=>setPalette(false)} title="Hide library">‹</button></div><input aria-label="Search pipeline tools" placeholder="Search tools…" value={query} onChange={e=>setQuery(e.target.value)}/>{categories.map(category=>{const specs=Object.values(catalog).filter(s=>!s.hidden&&s.category===category&&`${s.label} ${s.description}`.toLowerCase().includes(query.toLowerCase()));return specs.length?<section key={category}><h3>{category}</h3>{specs.map(s=><button key={s.key} className="palette-tool" draggable onDragStart={e=>e.dataTransfer.setData('application/studio-node',s.key)} onClick={()=>addNode(s.key)} title={s.description}><span style={{color:s.color}}>◇</span><span>{s.label}</span><span className="add-plus">+</span></button>)}</section>:null;})}<p className="library-help">Drag a tool onto the canvas, or click to add it.</p></aside>}
      <main className="flow-region" onDragOver={e=>{e.preventDefault();e.dataTransfer.dropEffect='move';}} onDrop={e=>{e.preventDefault();const kind=e.dataTransfer.getData('application/studio-node');if(catalog[kind])addNode(kind,rf.screenToFlowPosition({x:e.clientX,y:e.clientY}));}}>
        {!palette&&!focus&&<button className="show-library" onClick={()=>setPalette(true)}>+ Tools</button>}
        <ReactFlow nodes={displayNodes} edges={edges} nodeTypes={nodeTypes} onNodesChange={changes=>{if(changes.some(c=>c.type==='remove'))pushHistory();setNodes(prev=>applyNodeChanges(changes,prev));}} onEdgesChange={changes=>{if(changes.some(c=>c.type==='remove'))pushHistory();setEdges(prev=>applyEdgeChanges(changes,prev));}}
          onNodeDragStart={()=>{dragSnapshot.current=snapshot();}} onNodeDragStop={()=>{if(dragSnapshot.current){pushHistory(dragSnapshot.current);dragSnapshot.current=null;}}}
          onConnect={onConnect} isValidConnection={compatible} onNodeClick={(event,n)=>{if(!event.target.closest('.react-flow__handle')){setSelection(n.id);setHelpOpen(false);if(!event.shiftKey&&!event.ctrlKey&&!event.metaKey)emit({type:'edit_node',node_id:n.id});}}} onPaneClick={()=>setSelection(null)} onEdgeDoubleClick={(_,edge)=>{pushHistory();setEdges(es=>es.filter(e=>e.id!==edge.id));}}
          onMoveEnd={(_,v)=>{viewport.current=v;emit();}} minZoom={.18} maxZoom={2} defaultViewport={initial.viewport} snapToGrid snapGrid={[12,12]} selectionKeyCode="Shift" multiSelectionKeyCode="Shift" panOnDrag={[0,1,2]} deleteKeyCode={['Backspace','Delete']} connectionRadius={28} defaultEdgeOptions={{type:'default',pathOptions:{curvature:.12},style:{stroke:'#8ca99b',strokeWidth:2.2},interactionWidth:22}}>
          <Background gap={24} size={1.2} color="#cad6cf"/><Controls showInteractive={false} fitViewOptions={{padding:.4,maxZoom:1}}/>{showMap&&<MiniMap style={{pointerEvents:'none'}} nodeColor={n=>catalog[n.data.kind]?.color||'#799687'} maskColor="rgba(246,249,246,.72)"/>}
        </ReactFlow>
        {notice&&<div className="canvas-notice" role="status">{notice}</div>}
      </main>
      {!focus&&helpOpen&&<aside className="node-inspector"><div className="panel-heading"><span>Build without losing work</span><button onClick={()=>setHelpOpen(false)} title="Close help">×</button></div><div className="inspector-welcome"><div className="welcome-icon">⌘</div><h2>Connect your experiment</h2><p>Move the blocks freely. Drag between sockets to pass structures and sequences to the next step.</p><p>Click a block to open its familiar input form. Apply changes to save them in this draft.</p><div className="non-destructive-note">Saved versions and run results remain intact when you edit this draft.</div><dl><dt>Pan</dt><dd>Drag the background</dd><dt>Zoom</dt><dd>Scroll</dd><dt>Select several</dt><dd>Shift + drag</dd><dt>Undo / redo</dt><dd>Ctrl + Z / Shift + Z</dd><dt>Duplicate</dt><dd>Ctrl + D</dd><dt>Fit canvas</dt><dd>F</dd><dt>Remove wire</dt><dd>Double-click it</dd></dl></div></aside>}
    </div><footer className="canvas-footer"><span>{nodes.length} steps · {edges.length} connections</span><span>Edits affect this draft. Run creates a separate snapshot.</span><button onClick={()=>{const blob=new Blob([JSON.stringify(snapshot(),null,2)],{type:'application/json'}),url=URL.createObjectURL(blob),a=document.createElement('a');a.href=url;a.download='pipeline.json';a.click();URL.revokeObjectURL(url);}}>Export graph ↗</button></footer>
  </div>;
}

let mounted=false;
const root=createRoot(document.getElementById('root'));
window.addEventListener('message',event=>{if(event.data.type==='streamlit:render'){const args=event.data.args;root.render(<ReactFlowProvider><Editor args={args}/></ReactFlowProvider>);mounted=true;}});
post('streamlit:componentReady',{apiVersion:1});post('streamlit:setFrameHeight',{height:760});
