import React from 'react';
import {BaseEdge, MarkerType, Position, getSmoothStepPath} from '@xyflow/react';

const format=value=>Number(value).toLocaleString('en-US');

export function studyInfo(kind, params, nodes){
  const settings=key=>nodes.find(n=>n.data.kind===key)?.data.params||{};
  const docking=settings('active_docking'), decision=settings('active_decision');
  const batch=docking.batch_size||10000, selection=settings('active_select');
  const best=Math.round(batch*(selection.predicted_fraction_percent??80)/100);
  const info={
    active_seed:{eyebrow:'START · FIRST ROUND ONLY',summary:`Random ${format(params.sample_size||10000)} compounds`,detail:'Sample from the centroid library'},
    active_library:{eyebrow:'FULL SCREENING LIBRARY',summary:'SMILES + IDs + 1024-bit fingerprints',detail:'Used for prediction and random exploration'},
    active_target:{eyebrow:'DOCKING INPUT',summary:params.target_pdb?(params.target_name||'Imported receptor and docking box'):'Choose an imported receptor'},
    active_docking:{eyebrow:'01 · DOCK',summary:params.engine||'UniDock Pro',detail:`${format(batch)} compounds per round · ${docking.worker_nodes?.length?`${docking.worker_nodes.length} computers`:'choose computers'}`},
    active_accumulate:{eyebrow:'02 · ACCUMULATE & RANK',summary:'All rounds · best score per source and ID',detail:`Reorder all scores; compare cumulative top ${format(decision.top_n||100)}`},
    active_decision:{eyebrow:'03 · DECISION',summary:`Compare top ${format(params.top_n||100)} scores at each rank`,detail:`Tolerance ${params.score_tolerance??0.05} ${docking.score_units||'score units'} · ${params.stable_rounds||3} stable rounds`},
    active_binder_rule:{eyebrow:'RELABEL EVERY ROUND',summary:`Best ${params.binder_percent_start??10}% → ${params.binder_percent_target??0.01}%`,detail:`Log schedule · ${params.tightening_rounds||10} rounds`},
    active_train:{eyebrow:'04 · LEARN',summary:`1024 bits → 3 × ${format(params.hidden_nodes||3000)} dense nodes`,detail:'Train on accumulated docking results'},
    active_predict:{eyebrow:'05 · PREDICT',summary:'Predict binder probability',detail:'Across the full compound library'},
    active_select:{eyebrow:'06 · SELECT',summary:`${format(best)} best + ${format(batch-best)} random`,detail:'Both from the same library · exclude already docked IDs'},
    active_results:{eyebrow:'FINISH',summary:'Ranked hits · models · iteration history',detail:'Collect the retained poses and final tables'},
  };
  return info[kind]||{};
}

export const studyRole=kind=>({active_library:'input',active_seed:'input',active_target:'file',active_decision:'decision',active_results:'output',active_accumulate:'math',active_binder_rule:'math',active_select:'math'})[kind]||'tool';

// Ports follow the drawing: down the left column, up the right column,
// then a horizontal return from selection to docking.
const routes={
  'active_library:active_seed':{hidden:true},
  'active_seed:active_docking':{from:'bottom',to:'top',label:'Random initial sample'},
  'active_target:active_docking':{from:'left',to:'right',toOffset:20,source:true},
  'active_docking:active_accumulate':{from:'bottom',to:'top',label:'Docking scores'},
  'active_accumulate:active_decision':{from:'bottom',to:'top',label:'Evaluate convergence'},
  'active_decision:active_binder_rule':{from:'right',to:'left',label:'No'},
  'active_binder_rule:active_train':{from:'right',to:'left'},
  'active_train:active_predict':{from:'top',to:'bottom',label:'Trained model'},
  'active_predict:active_select':{from:'top',to:'bottom',label:'Top predictions'},
  'active_library:active_predict':{from:'right',to:'left',reference:true,source:true},
  'active_library:active_select':{from:'bottom',to:'top',label:'Random exploration',source:true},
  'active_decision:active_results':{from:'bottom',to:'top',label:'Yes · stop',stop:true},
  'active_select:active_docking':{from:'left',to:'right',label:'Best + random · next docking round',loop:true},
};

export function studyConnections(nodes,edges){
  const byId=new Map(nodes.map(n=>[n.id,n]));
  const handles=new Map(nodes.map(n=>[n.id,{incoming:[],outgoing:[]}]));
  const references=[];
  const rightmost=Math.max(...nodes.map(n=>n.position.x+(n.data.studySize?.width||320)));
  const topmost=Math.min(...nodes.map(n=>n.position.y));
  const displayEdges=edges.map(edge=>{
    const route=routes[`${byId.get(edge.source)?.data.kind}:${byId.get(edge.target)?.data.kind}`]||{from:'right',to:'left'};
    if(route.hidden)return {...edge,hidden:true};
    let source=edge.source;
    if(route.reference){
      const original=byId.get(edge.source),target=byId.get(edge.target);
      source=`input-reference-${edge.id}`;
      references.push({id:source,type:'studio',position:{x:target.position.x-305,y:target.position.y+30},draggable:false,deletable:false,
        data:{...original.data,reference_to:original.id,studySize:{width:230,height:120},studyInfo:{eyebrow:'INPUT · SAME LIBRARY',summary:'Full fingerprint dataset',detail:'Shared library selection'}}});
      handles.set(source,{incoming:[],outgoing:[]});
    }
    const sourceHandle=`study-out-${edge.id}`,targetHandle=`study-in-${edge.id}`;
    handles.get(source)?.outgoing.push({id:sourceHandle,side:route.from,offset:route.fromOffset??50});
    handles.get(edge.target)?.incoming.push({id:targetHandle,side:route.to,offset:route.toOffset??50});
    const color=route.loop?'#b76526':route.source?'#8b9da9':route.stop?'#4a8268':'#47748e';
    return {...edge,source,type:'study',sourceHandle,targetHandle,
      label:route.label,style:{stroke:color,strokeWidth:route.loop?3:route.source?1.8:2.5,...(route.source?{strokeDasharray:'5 5'}:{})},
      markerEnd:{type:MarkerType.ArrowClosed,color,width:18,height:18},
      data:{lane:route.lane,laneX:rightmost+78,laneY:topmost-52,loop:route.loop,source:route.source},
      ariaLabel:route.label||`Flow from ${byId.get(edge.source)?.data.label} to ${byId.get(edge.target)?.data.label}`};
  });
  return {handles,edges:displayEdges,references};
}

function roundedPath(points){
  let path=`M ${points[0][0]} ${points[0][1]}`;
  for(let i=1;i<points.length-1;i++){
    const [x,y]=points[i],prev=points[i-1],next=points[i+1];
    const before=Math.hypot(x-prev[0],y-prev[1]),after=Math.hypot(next[0]-x,next[1]-y);
    if(!before||!after)continue;
    const r=Math.min(15,before/2,after/2);
    path+=` L ${x-(x-prev[0])*r/before} ${y-(y-prev[1])*r/before} Q ${x} ${y} ${x+(next[0]-x)*r/after} ${y+(next[1]-y)*r/after}`;
  }
  return path+` L ${points.at(-1)[0]} ${points.at(-1)[1]}`;
}

export function StudyEdge({id,sourceX,sourceY,targetX,targetY,sourcePosition,targetPosition,markerEnd,style,label,data={}}){
  let [path,labelX,labelY]=getSmoothStepPath({sourceX,sourceY,targetX,targetY,sourcePosition,targetPosition,borderRadius:16,offset:28});
  if(data.lane==='top'){
    path=roundedPath([[sourceX,sourceY],[sourceX,data.laneY],[targetX,data.laneY],[targetX,targetY]]);
    labelX=(sourceX+targetX)/2;labelY=data.laneY;
  }else if(data.lane==='right'){
    path=roundedPath([[sourceX,sourceY],[data.laneX,sourceY],[data.laneX,targetY],[targetX,targetY]]);
    labelX=data.laneX-40;labelY=(sourceY+targetY)/2;
  }
  return <BaseEdge id={id} path={path} markerEnd={markerEnd} style={style} interactionWidth={22}
    label={label} labelX={labelX} labelY={labelY}
    labelStyle={{fontSize:13,fontWeight:data.loop?650:550,fill:data.loop?'#96501d':'#4b6574'}}
    labelBgStyle={{fill:data.loop?'#fff4e9':'#f7fafc',fillOpacity:1}} labelBgPadding={[7,5]} labelBgBorderRadius={5}/>;
}

export function studyPortStyle(handle){
  return ['top','bottom'].includes(handle.side)?{left:`${handle.offset}%`}:{top:`${handle.offset}%`};
}

export const studyPortPosition=side=>({top:Position.Top,bottom:Position.Bottom,left:Position.Left,right:Position.Right})[side];
