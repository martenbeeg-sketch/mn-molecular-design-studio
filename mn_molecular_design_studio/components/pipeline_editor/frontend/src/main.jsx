import React, {memo, useCallback, useEffect, useRef, useState} from 'react';
import {createRoot} from 'react-dom/client';
import {ReactFlow, ReactFlowProvider, Background, Controls, MiniMap, Handle, Position, applyNodeChanges, applyEdgeChanges, addEdge, useReactFlow, useNodesInitialized, useUpdateNodeInternals} from '@xyflow/react';
import {frameBounds, originalConnection, projectGroups, pruneGroups} from './canvasGroups.js';
import {ConnectionPanel} from './connectionPanel.jsx';
import {ToolPalette, artifactLabel, compatibleInlets, portsFor} from './toolPalette.jsx';
import {CollapsedGroup, Frame, RoundOutline} from './workspaceBlocks.jsx';
import '@xyflow/react/dist/style.css';
import './style.css';

const post=(type,data={})=>window.parent.postMessage({isStreamlitMessage:true,type,...data},'*');
const copy=x=>JSON.parse(JSON.stringify(x));
const uid=()=>crypto.randomUUID().replaceAll('-','');
let catalog={};
const emittedTypes=(node,handle,edges,nodes,seen=new Set())=>{
  if(!node||seen.has(node.id))return [];
  const spec=catalog[node.data.kind]||{},port=(spec.capabilities?.outputs||spec.outputs||[]).find(item=>item.id===handle);
  if(node.data.kind!=='filter')return port?(port.artifact_types||[port.type].filter(Boolean)):[];
  seen.add(node.id);
  const parents=edges.filter(edge=>edge.target===node.id);
  const inherited=parents.flatMap(edge=>emittedTypes(nodes.find(item=>item.id===edge.source),edge.sourceHandle||'out',edges,nodes,new Set(seen)));
  return [...new Set(inherited.length?inherited:(port?.artifact_types||[port?.type].filter(Boolean)))];
};
const authored=(nodes,edges,base,viewport)=>({...base,nodes:nodes.map(n=>({id:n.id,type:'studio',position:n.position,data:{kind:n.data.kind,label:n.data.label,...(n.data.description?{description:n.data.description}:{}),params:n.data.params,...Object.fromEntries(['enabled','retained_run_dir','execution','evidence_namespace'].filter(k=>k in n.data).map(k=>[k,n.data[k]]))}})),edges:edges.map(e=>({id:e.id,source:e.source,target:e.target,sourceHandle:e.sourceHandle||'out',targetHandle:e.targetHandle||'in',...(e.loop_back?{loop_back:true}:{})})),viewport});
const statusLabels={completed:'Complete',cached:'Reused',running:'Running',queued:'Queued',failed:'Failed',blocked:'Blocked',pending:'Pending'};

function SocketList({ports, side, counts, onBranch}) {
  const input = side === 'input';
  if (!ports.length) return <div aria-hidden="true"/>;
  return <div className={'node-socket-column ' + side}>
    <span className="node-socket-heading">{input ? 'Inputs' : 'Outputs'}</span>
    {ports.map(port => {
      const types = port.artifact_types.map(artifactLabel);
      const detail = input ? (port.multiple ? 'Multiple links' : !port.required ? 'Wire or configure' : '') : '';
      const count = counts?.[port.id] || 0;
      return <div key={port.id} className="node-socket-row" data-port-id={port.id} data-linked={count > 0}>
        <Handle id={port.id} type={input ? 'target' : 'source'} position={input ? Position.Left : Position.Right}
          className="typed-socket" title={port.label + (types.length ? ' · ' + types.join(', ') : '')}/>
        <span className="node-socket-text">
          <strong title={port.label}>{port.label}</strong>
          <small title={[...types, detail].filter(Boolean).join(' · ')}>
            {count > 1 ? count + ' links · ' : ''}{types.length === 1 ? types[0] : types.length + ' artifact types'}
          </small>
        </span>
        {!input && onBranch && <button className="socket-branch nodrag nopan" type="button" title={'Add step after ' + port.label}
          aria-label={'Add step after ' + port.label} onClick={event => {event.stopPropagation();onBranch(port.id);}}>+</button>}
      </div>;
    })}
  </div>;
}

const StudioNode=memo(({id,data,selected})=>{
  const updateNodeInternals=useUpdateNodeInternals();
  useEffect(()=>updateNodeInternals(id),[id,data.kind,updateNodeInternals]);
  const spec=catalog[data.kind]||{}; const p=data.params||{};
  const summary=data.summary||(p.num_designs?`${p.num_designs} design${p.num_designs===1?'':'s'} · ${p.binder_length||'variable'} aa`:p.sequences?`${p.sequences} sequence${p.sequences===1?'':'s'} / backbone`:data.kind==='target'?(p.target_pdb?`${p.target_chains||'A'} · ${p.hotspots?'hotspots set':'no hotspots'} · ${p.coldspots?'coldspots set':'no coldspots'}`:'Choose a target structure'):data.kind==='binder_context'?`${p.hotspots?'hotspot override':p.clear_unselected?'no hotspots':'inherit hotspots'} · ${p.coldspots?'coldspot override':p.clear_unselected?'no coldspots':'inherit coldspots'}`:data.kind==='filter'?(p.metric?`${p.metric} ${p.operator} ${p.threshold}`:p.passing_only?'Filter-passing designs':'All candidates'):data.kind==='results'?'Overlay · confidence · RMSD':spec.category);
  const coldspotLabel={native:'Coldspots supported',mapped:'Coldspots supported · mapped',unsupported:'Coldspots not supported'}[spec.coldspot_support];
  const hotspotLabel={native:'Hotspots supported',conditional:'Hotspots depend on task',unsupported:'Hotspots not supported'}[spec.hotspot_support];
  const targetContextLabel=spec.target_context_support==='structure_only'?'Structure + chains only · interface annotations ignored':null;
  return <div className={`studio-node ${selected?'is-selected':''}`} style={{'--node-color':spec.color}}>
    <div className="node-category"><span className="category-dot"/>{spec.category}{data.status&&<span className={`node-status ${data.status}`}>{statusLabels[data.status]||data.status}</span>}</div>
    <div className="node-title">{data.label||spec.label}</div>
    <div className="node-summary">{summary}</div><div className="node-edit-hint">Select to inspect <span>↗</span></div>
    {coldspotLabel&&<div className={`node-capability ${spec.coldspot_support==='unsupported'?'unsupported':'supported'}`}>{coldspotLabel}</div>}
    {hotspotLabel&&<div className={`node-capability ${spec.hotspot_support==='unsupported'?'unsupported':'supported'}`}>{hotspotLabel}</div>}
    {targetContextLabel&&<div className="node-capability unsupported">{targetContextLabel}</div>}
    {data.count!==undefined&&<div className="node-count">{data.count} candidates</div>}
    <div className="node-sockets">
      <SocketList ports={portsFor(spec, 'inputs')} side="input" counts={data.inputCounts}/>
      <SocketList ports={portsFor(spec, 'outputs')} side="output" counts={data.outputCounts} onBranch={data.onBranch}/>
    </div>
  </div>;
});
const nodeTypes={studio:StudioNode, collapsedGroup:CollapsedGroup, frame:Frame};
const matrixEmpty={wires:{},handoffs:{}};

function Editor({args}){
  catalog=args.catalog;
  const initial=args.graph;
  const readOnly=Boolean(args.template_read_only);
  const matrix=args.connections || matrixEmpty;
  const feedback = (source, output, target, input) => matrix.screening_feedback?.join('/') === [source, output, target, input].join('/');
  const fitPadding=.4;
  const [nodes,setNodes]=useState(initial.nodes),[edges,setEdges]=useState(initial.edges);
  const [groups,setGroups]=useState(initial.layout_groups || []);
  const [panel,setPanel]=useState('');
  const [pipelineName,setPipelineName]=useState(initial.name||'Untitled pipeline');
  const [notice,setNotice]=useState(''),[selection,setSelection]=useState(null),[branch,setBranch]=useState(null);
  const [palette,setPalette]=useState(true),[focus,setFocus]=useState(false),[historyTick,setHistoryTick]=useState(0),[helpOpen,setHelpOpen]=useState(false),[showMap,setShowMap]=useState(false);
  const [fullscreen,setFullscreen]=useState(false),shellRef=useRef(null);
  const base=useRef(initial),viewport=useRef(initial.viewport||{x:0,y:0,zoom:1}),reset=useRef(args.reset_token),lastSent=useRef(''),pipelineNameRef=useRef(pipelineName);
  pipelineNameRef.current=pipelineName;
  const state=useRef({nodes,edges,groups}),history=useRef({past:[],future:[]}),clipboard=useRef(null),dragSnapshot=useRef(null);
  const groupEndpoints=useRef(new Map());
  const rf=useReactFlow(), initialized=useNodesInitialized(), fitPending=useRef(!args.restore_viewport); state.current={nodes,edges,groups};
  const snapshot=useCallback(()=>authored(state.current.nodes,state.current.edges,{...base.current,name:pipelineNameRef.current,layout_groups:pruneGroups(state.current.groups,state.current.nodes)},viewport.current),[]);
  const emit=useCallback((action=null)=>{
    if(['save','save_as','save_template','use_template'].includes(action))viewport.current=rf.getViewport();
    const graph=snapshot();
    lastSent.current=JSON.stringify(graph);
    post('streamlit:setComponentValue',{value:{graph,...(action?{action:{...(typeof action==='string'?{type:action}:action),id:uid()}}:{})}});
  },[rf,snapshot]);
  const pushHistory=useCallback((snap=null)=>{const h=history.current;h.past.push(snap||snapshot());if(h.past.length>60)h.past.shift();h.future=[];setHistoryTick(x=>x+1);},[snapshot]);
  const restore=useCallback(g=>{base.current={...base.current,name:g.name};setPipelineName(g.name||'Untitled pipeline');setNodes(copy(g.nodes));setEdges(copy(g.edges));setGroups(copy(g.layout_groups||[]));if(g.viewport){viewport.current=g.viewport;rf.setViewport(g.viewport);}},[rf]);
  const undo=useCallback(()=>{const h=history.current;if(!h.past.length)return;h.future.push(snapshot());restore(h.past.pop());setHistoryTick(x=>x+1);},[snapshot,restore]);
  const redo=useCallback(()=>{const h=history.current;if(!h.future.length)return;h.past.push(snapshot());restore(h.future.pop());setHistoryTick(x=>x+1);},[snapshot,restore]);
  useEffect(()=>{
    if(reset.current!==args.reset_token){reset.current=args.reset_token;base.current=args.graph;setPipelineName(args.graph.name||'Untitled pipeline');setNodes(copy(args.graph.nodes));setEdges(copy(args.graph.edges));setGroups(copy(args.graph.layout_groups||[]));history.current={past:[],future:[]};setSelection(null);viewport.current=args.graph.viewport;fitPending.current=!args.restore_viewport;if(args.restore_viewport&&viewport.current)rf.setViewport(viewport.current);}
    base.current={...base.current,name:args.graph.name};
    if(args.graph.name!==pipelineNameRef.current)setPipelineName(args.graph.name||'Untitled pipeline');
  },[args.reset_token,args.graph.name,rf]);
  const settingsApplied=useRef(args.settings_update?.id);
  useEffect(()=>{
    const change=args.settings_update;
    if(!change||settingsApplied.current===change.id)return;
    settingsApplied.current=change.id;
    const current=state.current.nodes.find(n=>n.id===change.node_id);
    const data=change.data||Object.fromEntries(['params','label','kind','description','execution','enabled','retained_run_dir','evidence_namespace'].filter(k=>k in change).map(k=>[k,change[k]]));
    if(!current||Object.entries(data).every(([key,value])=>JSON.stringify(current.data[key])===JSON.stringify(value)))return;
    pushHistory();
    setNodes(prev=>prev.map(n=>n.id===change.node_id?{...n,data:{...n.data,...copy(data)}}:n));
  },[args.settings_update?.id,pushHistory]);
  useEffect(()=>{
    if(!args.input_contexts)return;
    setNodes(prev=>prev.map(n=>args.input_contexts[n.id]?{...n,data:{...n.data,params:{...n.data.params,...args.input_contexts[n.id]}}}:n));
  },[JSON.stringify(args.input_contexts)]);
  useEffect(()=>{
    if(readOnly)return;
    const timer=setTimeout(()=>{const g=snapshot(),serialized=JSON.stringify(g);if(serialized!==lastSent.current){lastSent.current=serialized;post('streamlit:setComponentValue',{value:{graph:g}});}},500);
    return()=>clearTimeout(timer);
  },[nodes,edges,groups,pipelineName,readOnly,snapshot]);
  useEffect(()=>{if(initialized&&fitPending.current){const id=requestAnimationFrame(()=>{rf.fitView({padding:fitPadding,maxZoom:1});fitPending.current=false;});return()=>cancelAnimationFrame(id);}},[initialized,nodes,rf,args.reset_token,fitPadding]);
  useEffect(()=>{post('streamlit:setFrameHeight',{height:focus?950:760});},[focus]);
  useEffect(()=>{
    const updateFullscreen=()=>{
      const active=(document.fullscreenElement||document.webkitFullscreenElement)===shellRef.current;
      setFullscreen(active);
      requestAnimationFrame(()=>rf.fitView({padding:fitPadding,duration:180}));
    };
    document.addEventListener('fullscreenchange',updateFullscreen);
    document.addEventListener('webkitfullscreenchange',updateFullscreen);
    return()=>{document.removeEventListener('fullscreenchange',updateFullscreen);document.removeEventListener('webkitfullscreenchange',updateFullscreen);};
  },[rf,fitPadding]);
  useEffect(()=>{if(notice){const t=setTimeout(()=>setNotice(''),5500);return()=>clearTimeout(t);}},[notice]);
  const addNode=useCallback((kind,position,connection=null)=>{
    const spec=catalog[kind];
    if(!spec||spec.hidden)return false;
    if(spec.screening_stage && state.current.nodes.some(node=>node.data.kind===kind)){
      setNotice(spec.label+' is already on the canvas. Connect to its existing input.');
      return false;
    }
    const suggested=args.screening_positions?.[kind];
    pushHistory();const pos=position||(suggested&&{x:suggested.x,y:suggested.y})||rf.screenToFlowPosition({x:window.innerWidth*.49,y:330});
    const offsets=[[0,0],[0,210],[320,0],[-320,0],[0,-210],[320,210],[-320,210],[0,420],[640,0]];const open=offsets.map(([x,y])=>({x:pos.x+x,y:pos.y+y})).find(p=>!state.current.nodes.some(n=>Math.abs(n.position.x-p.x)<280&&Math.abs(n.position.y-p.y)<170))||{x:pos.x,y:pos.y+630};
    const n={id:`${kind}-${uid().slice(0,10)}`,type:'studio',position:open,selected:true,data:{kind,label:spec.label,params:Object.fromEntries(spec.fields.map(f=>[f.name,f.default]))}};
    setNodes(prev=>[...prev.map(n=>({...n,selected:false})),n]);
    if(connection)setEdges(prev=>[...prev,{...connection,id:`edge-${uid()}`,target:n.id}]);
    setSelection(connection?.source||n.id);
    setNotice(connection ? `${spec.label} added and connected.` : `${spec.label} added. Connect its sockets to another step.`);
    return true;
  },[args.screening_positions,pushHistory,rf]);
  const normalizeConnection=connection=>{
    const current=originalConnection(connection,groupEndpoints.current);
    return {...current,sourceHandle:current.sourceHandle||'out',targetHandle:current.targetHandle||'in'};
  };
  const compatible=useCallback(connection=>{
    const ns=state.current.nodes,es=state.current.edges;
    connection=normalizeConnection(connection);
    const source=ns.find(n=>n.id===connection.source),target=ns.find(n=>n.id===connection.target);
    if(!source||!target||source.id===target.id)return false;
    const a=catalog[source.data.kind],b=catalog[target.data.kind];
    const outputs=a.capabilities?.outputs||a.outputs||[],inputs=b.capabilities?.inputs||b.inputs||[];
    const output=outputs.find(port=>port.id===connection.sourceHandle);
    const input=inputs.find(port=>port.id===connection.targetHandle);
    if(!output||!input||!matrix.wires[source.data.kind]?.[output.id]?.[target.data.kind]?.[input.id])return false;
    const produced=emittedTypes(source,connection.sourceHandle,es,ns),accepts=input.artifact_types||input.accepts||[];
    if(!produced.some(type=>accepts.includes(type)))return false;
    const occupied=es.filter(e=>e.target===target.id&&(e.targetHandle||'in')===input.id);
    if(occupied.length&&(!input.multiple||occupied.some(e=>e.source===source.id&&(e.sourceHandle||'out')===output.id)))return false;
    if(feedback(source.data.kind,output.id,target.data.kind,input.id))return true;
    const visit=(id,seen=new Set())=>{if(id===source.id)return true;if(seen.has(id))return false;seen.add(id);return es.filter(e=>e.source===id&&!e.loop_back).some(e=>visit(e.target,seen));};
    return !visit(target.id);
  },[matrix]);
  const onConnect=c=>{
    if(!compatible(c)){setNotice('These sockets are incompatible, occupied, or would create a cycle.');return;}
    pushHistory();
    const wire=normalizeConnection(c);
    const source=nodes.find(node=>node.id===wire.source), target=nodes.find(node=>node.id===wire.target);
    if(feedback(source.data.kind,wire.sourceHandle,target.data.kind,wire.targetHandle))wire.loop_back=true;
    setEdges(prev=>addEdge({...wire,id:`edge-${uid()}`},prev));
    return true;
  };
  const connectNext=(source, output, kind, inlet)=>{
    const match=matrix.wires[source.data.kind]?.[output]?.[kind]?.[inlet];
    const accepted=catalog[kind]?.capabilities?.inputs.find(port=>port.id===inlet)?.artifact_types||[];
    if(!match||!emittedTypes(source,output,edges,nodes).some(type=>accepted.includes(type)))return false;
    const position=catalog[kind].screening_stage ? null : {x:source.position.x+352,y:source.position.y};
    return addNode(kind,position,{source:source.id,sourceHandle:output,targetHandle:inlet,
      ...(feedback(source.data.kind,output,kind,inlet)?{loop_back:true}:{})});
  };
  const groupedMembers=new Set(groups.flatMap(group=>group.members));
  const selectionMembers=nodes.filter(node=>node.selected&&!groupedMembers.has(node.id));
  const selectedGroup=groups.find(group=>group.id===selection);
  const groupSelection=()=>{
    if(selectionMembers.length<2||groups.length>=50)return;
    pushHistory();
    const id='group-'+uid();
    setGroups(items=>[...items,{id,name:'New group',members:selectionMembers.map(node=>node.id),collapsed:false}]);
    setSelection(id);setPanel('group');
  };
  const updateGroup=(id,values)=>{
    pushHistory();setGroups(items=>items.map(item=>item.id===id?{...item,...values}:item));
  };
  const removeGroup=id=>{
    pushHistory();setGroups(items=>items.filter(item=>item.id!==id));setSelection(null);setPanel('');
  };
  const duplicateGroup=group=>{
    const members=new Set(group.members);
    const originals=nodes.filter(node=>members.has(node.id));
    const internal=edges.filter(edge=>members.has(edge.source)&&members.has(edge.target));
    if(groups.length>=50||nodes.length+originals.length>100||edges.length+internal.length>250){
      setNotice('This group would exceed the pipeline size limit.');return;
    }
    const ids=new Map(originals.map(node=>[node.id,node.data.kind+'-'+uid().slice(0,10)]));
    const groupId='group-'+uid();
    const offset=group.collapsed?318:frameBounds(originals).width+35;
    pushHistory();
    setNodes(previous=>[...previous.map(node=>({...node,selected:false})),...originals.map(node=>({
      ...copy(node),id:ids.get(node.id),selected:false,position:{x:node.position.x+offset,y:node.position.y},
    }))]);
    setEdges(previous=>[...previous,...internal.map(edge=>({...edge,id:'edge-'+uid(),source:ids.get(edge.source),target:ids.get(edge.target)}))]);
    setGroups(previous=>[...previous,{...group,id:groupId,members:[...ids.values()]}]);
    setSelection(groupId);
  };
  const changeCanvasNodes=changes=>{
    const real=changes.filter(change=>nodes.some(node=>node.id===change.id));
    const shifts=changes.filter(change=>change.type==='position'&&change.position&&groups.some(group=>group.id===change.id));
    if(real.some(change=>change.type==='remove'))pushHistory();
    if(!real.length&&!shifts.length)return;
    setNodes(previous=>{
      let next=applyNodeChanges(real,previous);
      for(const shift of shifts){
        const group=groups.find(item=>item.id===shift.id);
        const members=next.filter(node=>group.members.includes(node.id));
        if(members.length<2)continue;
        const {position}=frameBounds(members), dx=shift.position.x-position.x, dy=shift.position.y-position.y;
        next=next.map(node=>group.members.includes(node.id)?{...node,position:{x:node.position.x+dx,y:node.position.y+dy}}:node);
      }
      return next;
    });
    if(real.some(change=>change.type==='remove')){
      const remaining=nodes.filter(node=>!real.some(change=>change.id===node.id&&change.type==='remove'));
      setGroups(items=>pruneGroups(items,remaining));
    }
  };
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
      if(e.key.toLowerCase()==='f'){e.preventDefault();rf.fitView({padding:fitPadding,duration:250});}
      if(e.key==='Escape')setSelection(null);
      if(readOnly)return;
      const mod=e.ctrlKey||e.metaKey;
      if(mod&&e.key.toLowerCase()==='z'){e.preventDefault();e.shiftKey?redo():undo();}
      if(mod&&e.key.toLowerCase()==='y'){e.preventDefault();redo();}
      if(mod&&e.key.toLowerCase()==='d'){e.preventDefault();duplicate();}
      if(mod&&e.key.toLowerCase()==='c'){const ns=state.current.nodes.filter(n=>n.selected);clipboard.current={nodes:copy(ns),edges:copy(state.current.edges.filter(e=>ns.some(n=>n.id===e.source)&&ns.some(n=>n.id===e.target)))};}
      if(mod&&e.key.toLowerCase()==='v'&&clipboard.current?.nodes.length){e.preventDefault();pushHistory();const c=clipboard.current,ids=Object.fromEntries(c.nodes.map(n=>[n.id,`${n.data.kind}-${uid().slice(0,10)}`]));setNodes(p=>[...p.map(n=>({...n,selected:false})),...c.nodes.map(n=>({...copy(n),id:ids[n.id],position:{x:n.position.x+70,y:n.position.y+70},selected:true}))]);setEdges(p=>[...p,...c.edges.map(e=>({...e,id:`edge-${uid()}`,source:ids[e.source],target:ids[e.target]}))]);}
    };window.addEventListener('keydown',handler);return()=>window.removeEventListener('keydown',handler);
  },[undo,redo,duplicate,pushHistory,rf,fitPadding,readOnly]);
  const layout=()=>{
    pushHistory();
    const ns=state.current.nodes,es=state.current.edges;
    const levels=new Map(),visiting=new Set();
    const level=id=>{
      if(levels.has(id))return levels.get(id);
      if(visiting.has(id))return 0;
      visiting.add(id);
      const ins=es.filter(e=>e.target===id&&!e.loop_back);
      const value=ins.length?Math.max(...ins.map(e=>level(e.source)))+1:0;
      levels.set(id,value);
      return value;
    };
    const columns=new Map();
    ns.forEach(node=>{
      const column=level(node.id);
      columns.set(column,[...(columns.get(column)||[]),node]);
    });
    const maxRows=Math.max(1,...[...columns.values()].map(group=>group.length));
    const rowGap=220,top=100,nodePositions=new Map();
    for(const [column,group] of columns){
      const ordered=[...group].sort((a,b)=>a.position.y-b.position.y||a.id.localeCompare(b.id));
      const startY=top+(maxRows-ordered.length)*rowGap/2;
      ordered.forEach((node,row)=>nodePositions.set(node.id,{x:70+column*360,y:startY+row*rowGap}));
    }
    setNodes(ns.map(node=>({...node,position:nodePositions.get(node.id)})));
    setTimeout(()=>rf.fitView({padding:.4,duration:250}),30);
  };
  const toggleFullscreen=async()=>{
    const element=shellRef.current;
    if(!element)return;
    try{
      if((document.fullscreenElement||document.webkitFullscreenElement)===element){
        const exit=document.exitFullscreen||document.webkitExitFullscreen;
        if(exit)await exit.call(document);
      }else{
        if(element.requestFullscreen)await element.requestFullscreen({navigationUI:'hide'});
        else if(element.webkitRequestFullscreen)element.webkitRequestFullscreen();
        else setNotice('Fullscreen is not supported by this browser.');
      }
    }catch{
      setNotice('Fullscreen could not be opened. Check that the browser allows fullscreen for this page.');
    }
  };
  // Keep the node library visible while fullscreen. Node forms are rendered by
  // Streamlit outside this iframe, so clicking a node first exits fullscreen.
  const canvasOnly=focus&&!fullscreen;
  const openBranch = (nodeId, port) => {
    setBranch({nodeId, port});
    setPalette(true);
    setFocus(false);
    setSelection(nodeId);
    setPanel('');
    setHelpOpen(false);
  };
  const source = nodes.find(node => node.id === branch?.nodeId);
  let branchData = null;
  if(source && portsFor(catalog[source.data.kind], 'outputs').some(port => port.id === branch.port)){
    const inlets=compatibleInlets(catalog,matrix.wires,source.data.kind,branch.port,
      emittedTypes(source,branch.port,edges,nodes));
    const existing=nodes.flatMap(target => (inlets[target.data.kind] || []).flatMap(port => {
      const wire={source:source.id,sourceHandle:branch.port,target:target.id,targetHandle:port.id};
      return compatible(wire) ? [{wire,label:target.data.label,kind:target.data.kind,port:port.label}] : [];
    }));
    branchData={...branch,node:source,inlets,existing};
  }
  const usedStages = new Set(nodes.filter(node=>catalog[node.data.kind].screening_stage).map(node=>node.data.kind));
  const displayEdges=edges.map(edge => edge.loop_back
    ? {...edge,label:'Next round',style:{stroke:'#ab641f',strokeWidth:2.5,strokeDasharray:'6 4'}}
    : edge);
  const inputCounts = Object.create(null), outputCounts = Object.create(null);
  for (const edge of edges) {
    const input = inputCounts[edge.target] ||= {};
    const output = outputCounts[edge.source] ||= {};
    const targetPort = edge.targetHandle || 'in', sourcePort = edge.sourceHandle || 'out';
    input[targetPort] = (input[targetPort] || 0) + 1;
    output[sourcePort] = (output[sourcePort] || 0) + 1;
  }
  const displayNodes=nodes.map(n=>{
    const incomingKinds=edges.filter(e=>e.target===n.id).map(e=>nodes.find(source=>source.id===e.source)?.data.kind);
    const hasSurfaceAnalysis=incomingKinds.includes('surface_analysis');
    const hasOtherResults=incomingKinds.some(kind=>kind&&kind!=='surface_analysis');
    const hasMDAnalysis=incomingKinds.includes('md_analysis');
    const mdStep=['md_source','md_analysis'].includes(n.data.kind)||(n.data.kind==='results'&&hasMDAnalysis);
    const summary=n.data.kind==='md_source'
      ? n.data.params.option_id?'One recorded segment · read-only':'Choose a recorded MD segment'
      : n.data.kind==='md_analysis'?'RMSD · RMSF · Rg · contacts'
      : n.data.kind==='results'&&hasMDAnalysis?'MD plots · tables · source provenance'
      : n.data.kind==='results'&&hasSurfaceAnalysis
      ? hasOtherResults?'Candidates + surface evidence':'Residue consensus · surface view'
      : undefined;
    return {...n,data:{...n.data,status:args.execution?.[n.id]?.status,
      count:mdStep?undefined:args.execution?.[n.id]?.candidate_count,
      summary:summary || args.round_outline?.screen?.step_details[n.id],
      inputCounts:inputCounts[n.id],outputCounts:outputCounts[n.id],
      onBranch:readOnly?undefined:port=>openBranch(n.id,port)}};
  });
  const projection=projectGroups(displayNodes,displayEdges,edges,groups,catalog);
  groupEndpoints.current=projection.endpoints;
  const focusedNode=nodes.find(node=>node.id===selection);
  const openSettings=async id=>{
    if((document.fullscreenElement||document.webkitFullscreenElement)===shellRef.current){
      const exit=document.exitFullscreen||document.webkitExitFullscreen;
      if(exit)await exit.call(document);
    }
    emit({type:'edit_node',node_id:id});
  };
  const inspectNode=(event,node)=>{
    if(event.target.closest('.react-flow__handle'))return;
    setHelpOpen(false);
    if(node.type==='frame'||node.type==='collapsedGroup'){
      setSelection(node.id);setPanel('group');return;
    }
    const id=node.id;
    setSelection(id);
    if(event.shiftKey||event.ctrlKey||event.metaKey)return;
    if(canvasOnly&&!readOnly)openSettings(id);
    else{
      if(canvasOnly)setFocus(false);
      setPanel('connections');
    }
  };
  return <div ref={shellRef} className={`editor-shell ${focus?'focus-mode':''} ${fullscreen?'fullscreen-mode':''}`}>
    <header className="editor-toolbar">
      <div className="editor-heading">
        <div className="toolbar-title">
          <button className="library-back" onClick={()=>emit('library')}>‹ Pipeline library</button>
          <input className="pipeline-name-input" aria-label={args.shared_template?'Template name':'Pipeline name'} value={pipelineName} readOnly={readOnly}
            onChange={event=>setPipelineName(event.target.value)}
            onBlur={()=>setPipelineName(name=>name.trim()||'Untitled pipeline')} title={args.shared_template?'Template name':'Pipeline name'}/>
          <span className="draft-chip">{args.shared_template?'Shared template':args.draft_saved?'Saved pipeline':'Unsaved pipeline'}</span>
        </div>
        <div className="toolbar-actions primary-actions">
          {args.shared_template ? <>
            <button className={readOnly?'run-button':undefined} disabled={!nodes.length} onClick={()=>emit('use_template')}>Use template</button>
            {!readOnly&&<button className="run-button" disabled={!nodes.length} onClick={()=>emit('save')}>Save template</button>}
            </> : <>
              <button onClick={()=>emit('versions')}>Versions</button>
              <button onClick={()=>emit('runs')}>Runs</button>
              <button className="run-button" disabled={!nodes.length} title={!nodes.length?'Add a step before setting up a run':'Review inputs and execution before submitting'}
                onClick={()=>emit('run')}>Run setup →</button>
            </>}
        </div>
      </div>
      <div className="toolbar-actions canvas-actions">
        {!readOnly && <>
        <button title="Undo (Ctrl+Z)" disabled={!history.current.past.length} onClick={undo}>↶</button>
        <button title="Redo (Ctrl+Shift+Z)" disabled={!history.current.future.length} onClick={redo}>↷</button>
        <span className="toolbar-divider"/>
        <button onClick={layout} title="Arrange the current graph">Arrange</button>
        <button disabled={selectionMembers.length<2} onClick={groupSelection} title="Group selected steps">Group</button>
        </>}
        <button aria-pressed={panel==='connections'} onClick={()=>{setPanel(panel==='connections'?'':'connections');setHelpOpen(false);}}
          title="Inspect inputs, outputs and compatible steps">Connections</button>
        {args.round_outline?.screen&&<button aria-pressed={panel==='round'}
          onClick={()=>{setPanel(panel==='round'?'':'round');setHelpOpen(false);}}>Screening round</button>}
        <button title="Toggle canvas overview" onClick={()=>setShowMap(!showMap)}>Map</button>
        <span className="canvas-spacer"/>
        <button onClick={()=>{setFocus(!focus);setTimeout(()=>rf.fitView({padding:fitPadding,duration:200}),80);}}>{focus?'Exit focus':'Focus canvas'}</button>
        <button onClick={toggleFullscreen} title={fullscreen?'Exit fullscreen':'Enter fullscreen'}>{fullscreen?'Exit fullscreen':'⛶ Fullscreen'}</button>
        {!args.shared_template && <>
          <button title="Configure campaign objective and stopping limits" onClick={()=>emit('campaign')}>Campaign</button>
          <button title="Configure workstation, GPU and batching defaults" onClick={()=>emit('resources')}>Resources</button>
          <button onClick={()=>emit('save')}>Save</button>
          <button onClick={()=>emit('save_as')}>Save as…</button>
          {args.can_edit_templates && <button disabled={!nodes.length} onClick={()=>emit('save_template')}>Save as shared template</button>}
        </>}
        <button title="Canvas shortcuts and help" onClick={()=>{setHelpOpen(!helpOpen);setSelection(null);setPanel('');}}>?</button>
      </div>
    </header>
    <div className="workspace">
      {!readOnly && palette && !canvasOnly && <ToolPalette key={branchData?.nodeId || 'all'}
        catalog={catalog} branch={branchData} usedStages={usedStages} onConnect={wire=>{if(onConnect(wire))setBranch(null);}} onBranchPort={port=>setBranch({...branch,port})}
        onClearBranch={()=>setBranch(null)} onClose={()=>setPalette(false)}
        onAdd={(kind,inlet) => {
          if (!branchData) addNode(kind);
          else if (connectNext(branchData.node, branchData.port, kind, inlet)) setBranch(null);
        }}/>}
      <main className="flow-region" onDragOver={readOnly?undefined:e=>{e.preventDefault();e.dataTransfer.dropEffect='move';}}
        onDrop={readOnly?undefined:e=>{e.preventDefault();const kind=e.dataTransfer.getData('application/studio-node');if(catalog[kind])addNode(kind,rf.screenToFlowPosition({x:e.clientX,y:e.clientY}));}}>
        {!readOnly&&!palette&&!canvasOnly&&<button className="show-library" onClick={()=>setPalette(true)}>+ Tools</button>}
        <ReactFlow nodes={readOnly?projection.nodes.map(node=>({...node,draggable:false})):projection.nodes}
          edges={projection.edges} nodeTypes={nodeTypes} nodesConnectable={!readOnly}
          onNodesChange={changeCanvasNodes}
          onEdgesChange={changes=>{if(changes.some(c=>c.type==='remove'))pushHistory();setEdges(prev=>applyEdgeChanges(changes,prev));}}
          onNodeDragStart={()=>{dragSnapshot.current=snapshot();}} onNodeDragStop={()=>{if(dragSnapshot.current){pushHistory(dragSnapshot.current);dragSnapshot.current=null;}}}
          onConnect={readOnly?undefined:onConnect} isValidConnection={compatible}
          onNodeClick={inspectNode}
          onNodeDoubleClick={(_,node)=>{if(!readOnly&&node.type==='studio'&&!canvasOnly)openSettings(node.id);}}
          onPaneClick={()=>{setSelection(null);if(panel==='group'||panel==='connections')setPanel('');}}
          onEdgeDoubleClick={readOnly?undefined:(_,edge)=>{pushHistory();setEdges(es=>es.filter(e=>e.id!==edge.id));}}
          onMoveEnd={(_,v)=>{viewport.current=v;if(!readOnly)emit();}} minZoom={.18} maxZoom={2} defaultViewport={initial.viewport} snapToGrid snapGrid={[12,12]} selectionKeyCode="Shift" multiSelectionKeyCode="Shift" panOnDrag={[0,1,2]} deleteKeyCode={readOnly?null:['Backspace','Delete']} connectionRadius={28} defaultEdgeOptions={{type:'default',pathOptions:{curvature:.12},style:{stroke:'#8ca99b',strokeWidth:2.2},interactionWidth:22}}>
          <Background gap={24} size={1.2} color="#cad6cf"/><Controls showInteractive={false} fitViewOptions={{padding:fitPadding,maxZoom:1}}/>{showMap&&<MiniMap style={{pointerEvents:'none'}} nodeColor={n=>catalog[n.data.kind]?.color||'#799687'} maskColor="rgba(246,249,246,.72)"/>}
        </ReactFlow>
        {!nodes.length&&<div className="canvas-empty"><h2>Your pipeline starts here</h2>
          <p>Add tools from the library{args.guided_setup?' or start with a connected workflow.':'.'}</p>
          {!palette&&<button onClick={()=>setPalette(true)}>Show tools</button>}
          {args.guided_setup&&<button onClick={()=>emit('guided_setup')}>Guided setup</button>}
        </div>}
        {notice&&<div className="canvas-notice" role="status">{notice}</div>}
      </main>
      {!canvasOnly&&panel==='connections'&&<ConnectionPanel key={focusedNode?.id||'empty'}
        node={focusedNode} catalog={catalog} matrix={matrix}
        close={()=>setPanel('')} configure={readOnly?undefined:()=>openSettings(focusedNode.id)}
        branch={readOnly?undefined:port=>openBranch(focusedNode.id,port)}/>}
      {!canvasOnly&&panel==='round'&&<RoundOutline program={args.round_outline}
        configure={readOnly?id=>{setSelection(id);setPanel('connections');}:openSettings} close={()=>setPanel('')}/>}
      {!canvasOnly&&panel==='group'&&selectedGroup&&<aside className="node-inspector group-panel" aria-label="Group settings">
        <div className="panel-heading"><span>{readOnly?selectedGroup.name:'Group settings'}</span><button aria-label="Close group settings" onClick={()=>setPanel('')}>×</button></div>
        {!readOnly&&<label>Group name<input key={selectedGroup.id+selectedGroup.name} defaultValue={selectedGroup.name} maxLength={80}
          onKeyDown={event=>{if(event.key==='Enter')event.currentTarget.blur();}}
          onBlur={event=>{const name=event.target.value.trim()||'Group';if(name!==selectedGroup.name)updateGroup(selectedGroup.id,{name});}}/></label>}
        <button onClick={()=>{
          const collapsed=!selectedGroup.collapsed;
          if(readOnly)setGroups(items=>items.map(item=>item.id===selectedGroup.id?{...item,collapsed}:item));
          else updateGroup(selectedGroup.id,{collapsed});
        }}>{selectedGroup.collapsed?'Open group':'Collapse group'}</button>
        {!readOnly&&<>
          <button onClick={()=>duplicateGroup(selectedGroup)}>Duplicate group</button>
          <button onClick={()=>removeGroup(selectedGroup.id)}>Ungroup steps</button>
        </>}
      </aside>}
      {!canvasOnly&&helpOpen&&readOnly&&<aside className="node-inspector" aria-label="Template help">
        <div className="panel-heading"><span>Template canvas</span><button onClick={()=>setHelpOpen(false)} title="Close help">×</button></div>
        <p>Select a step to view its settings and connections. Pan by dragging the background and zoom by scrolling.</p>
      </aside>}
      {!canvasOnly&&helpOpen&&!readOnly&&<aside className="node-inspector"><div className="panel-heading"><span>Build without losing work</span><button onClick={()=>setHelpOpen(false)} title="Close help">×</button></div><div className="inspector-welcome"><div className="welcome-icon">⌘</div><h2>Connect your experiment</h2><p>Move the blocks freely. Drag between sockets to pass structures and sequences to the next step.</p><p>Select a step to inspect its connections. Double-click to open its input form.</p>{!args.shared_template&&<div className="non-destructive-note">Versions stores snapshots of this saved pipeline. Runs opens the execution history.</div>}<dl><dt>Pan</dt><dd>Drag the background</dd><dt>Zoom</dt><dd>Scroll</dd><dt>Select several</dt><dd>Shift + drag</dd><dt>Undo / redo</dt><dd>Ctrl + Z / Shift + Z</dd><dt>Duplicate</dt><dd>Ctrl + D</dd><dt>Fit canvas</dt><dd>F</dd><dt>Remove wire</dt><dd>Double-click it</dd></dl></div></aside>}
    </div><footer className="canvas-footer"><span>{nodes.length} steps · {edges.length} connections</span><span>{args.shared_template?'Use template creates a separate pipeline.':'Run setup checks inputs before submission.'}</span><button onClick={()=>{const blob=new Blob([JSON.stringify(snapshot(),null,2)],{type:'application/json'}),url=URL.createObjectURL(blob),a=document.createElement('a');a.href=url;a.download='pipeline.json';a.click();URL.revokeObjectURL(url);}}>Export graph ↗</button></footer>
  </div>;
}

let mounted=false;
const root=createRoot(document.getElementById('root'));
window.addEventListener('message',event=>{if(event.data.type==='streamlit:render'){const args=event.data.args;root.render(<ReactFlowProvider><Editor args={args}/></ReactFlowProvider>);mounted=true;}});
post('streamlit:componentReady',{apiVersion:1});post('streamlit:setFrameHeight',{height:760});
