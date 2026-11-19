(() => {
  'use strict';
  const data = JSON.parse(document.getElementById('assessment-data').textContent);
  const lookup = new Map(Object.values(data.residues).flat().map(residue => [residue.id, residue]));
  const palette = {target:'#477aa9',binder:'#318c70',hbond:'#148e91',salt:'#935bb5',contact:'#a7b9b4',hotspot:'#d89425'};
  const state = {color:'role', selected:[], title:'', viewer:null, model:null, hoverShapes:[], hoverLabels:[], surfaceEpoch:0};
  const byId = identifier => document.getElementById(identifier);
  const element = (tag, className, text) => {
    const node = document.createElement(tag);
    if (className) node.className = className;
    if (text !== undefined) node.textContent = text;
    return node;
  };
  const svgElement = (tag, attributes, text) => {
    const node = document.createElementNS('http://www.w3.org/2000/svg', tag);
    Object.entries(attributes || {}).forEach(([name,value]) => node.setAttribute(name,String(value)));
    if (text !== undefined) node.textContent = text;
    return node;
  };
  const label = identifier => lookup.get(identifier)?.label || 'Unmapped';
  const descriptive = identifier => {
    const residue = lookup.get(identifier);
    return residue ? `${residue.side === 'A' ? 'Target' : 'Binder'} ${residue.label} · saved ${residue.author}` : 'Unmapped residue';
  };
  const selection = identifiers => ({or:identifiers.map(identifier => {
    const residue = lookup.get(identifier);
    return residue ? {chain:residue.side,resi:residue.number} : {chain:'!unmapped'};
  })});
  const point = coordinates => ({x:coordinates[0],y:coordinates[1],z:coordinates[2]});
  const finite = value => typeof value === 'number' && Number.isFinite(value);
  const rounded = (value,precision=1) => finite(value) ? value.toFixed(precision) : '—';
  const confidenceColor = value => !finite(value) ? '#c7d0d3' : value >= 90 ? '#2d64af' : value >= 70 ? '#70b9d1' : value >= 50 ? '#e2c44d' : '#dc8c49';
  const burialColor = (value, maximum) => {
    const proportion = finite(value) ? Math.min(1,value / Math.max(1,maximum)) : 0;
    return '#' + [225-183*proportion,236-99*proportion,233-128*proportion].map(channel=>Math.round(channel).toString(16).padStart(2,'0')).join('');
  };
  byId('design-title').textContent = data.title;
  byId('target-label').textContent = data.target_label + ` · ${data.residues.B.length}-residue binder`;
  byId('context-note').textContent = 'Saved model only · author residue numbering · original admission unchanged. Hover across the diagram, sequence and 3D view; click to focus.';
  byId('scope-note').textContent = data.limitation;
  byId('source-digest').textContent = 'Source SHA256: ' + data.source_sha256;
  byId('numbering-note').textContent = data.target_caption;
  byId('technical-note').textContent = 'Contacts ≤4.5 Å heavy atoms · H-bonds use recorded atom geometry · salt bridges use charged-group centroids ≤4 Å. Hydrogen basis: ' + data.hbond_basis + '.';
  if(data.unmapped_interactions)byId('technical-note').append(` ${data.unmapped_interactions} atom-level overlays unavailable; no endpoints were guessed.`);
  if(data.extra_chains)byId('technical-note').append(` ${data.extra_chains} additional saved chains are shown in grey; analysis concerns the selected target and binder.`);
  if (!data.confidence_confirmed) {
    document.querySelector('[data-color="confidence"]').disabled = true;
    document.querySelector('[data-color="confidence"]').title = 'No confirmed per-residue binder pLDDT for this exact model.';
  }

  function clearHover() {
    if (!state.viewer) return;
    state.hoverShapes.forEach(shape => state.viewer.removeShape(shape));
    state.hoverLabels.forEach(annotation => state.viewer.removeLabel(annotation));
    state.hoverShapes = []; state.hoverLabels = [];
  }
  function highlight(identifiers, text, focus=false) {
    identifiers = identifiers.filter(identifier => lookup.has(identifier));
    document.querySelectorAll('.picked').forEach(node => node.classList.remove('picked'));
    document.querySelectorAll('[data-residue]').forEach(node => {
      if (identifiers.includes(node.dataset.residue)) node.classList.add('picked');
    });
    document.querySelectorAll('.network-edge').forEach(node => {
      if (identifiers.includes(node.dataset.target) || identifiers.includes(node.dataset.binder)) node.classList.add('picked');
    });
    byId('selection').textContent = text || 'Hover a residue or contact; click to focus it in 3D.';
    if (!state.viewer) return;
    clearHover();
    identifiers.forEach(identifier => {
      const residue = lookup.get(identifier), center = point(residue.xyz);
      state.hoverShapes.push(state.viewer.addSphere({center,radius:1.05,color:'#e2a526',alpha:.75}));
      state.hoverLabels.push(state.viewer.addLabel(residue.label,{position:{x:center.x,y:center.y+2,z:center.z},fontSize:11,fontColor:'#263d42',backgroundColor:'#ffedb3',backgroundOpacity:.9,borderThickness:0,inFront:true}));
    });
    if (identifiers.length === 2) {
      state.hoverShapes.push(state.viewer.addCylinder({start:point(lookup.get(identifiers[0]).xyz),end:point(lookup.get(identifiers[1]).xyz),radius:.08,color:'#e2a526',dashed:true}));
    }
    if (focus && identifiers.length) {
      state.viewer.zoomTo(selection(identifiers),450);
      state.viewer.zoom(.72);
    }
    state.viewer.render();
  }
  function bindResidues(node, identifiers, text) {
    node.addEventListener('pointerenter',() => highlight(identifiers,text));
    node.addEventListener('pointerleave',() => highlight(state.selected,state.title));
    node.addEventListener('click',() => {
      state.selected = identifiers.filter(identifier => lookup.has(identifier)); state.title = text;
      highlight(state.selected,state.title,true);
    });
    if (node.namespaceURI?.includes('svg')) {
      node.setAttribute('tabindex','0'); node.setAttribute('role','button'); node.setAttribute('aria-label',text);
      node.addEventListener('focus',() => highlight(identifiers,text));
      node.addEventListener('keydown',event => { if (event.key === 'Enter' || event.key === ' ') {event.preventDefault();node.dispatchEvent(new MouseEvent('click'));} });
    }
  }
  function contactText(contact) {
    return `${descriptive(contact.target)} ↔ ${descriptive(contact.binder)} · ${contact.kind} · minimum heavy-atom distance ${rounded(contact.distance,2)} Å · ${contact.pairs} atom pairs. Orange connector is a selection guide, not a chemical bond.`;
  }
  function interactionColor(kind) {
    return kind === 'Contact' ? palette.contact : kind === 'H-bond' ? palette.hbond : kind === 'Salt bridge' ? palette.salt : '#be7826';
  }

  function buildNetwork() {
    const target = [...new Set(data.contacts.map(row => row.target))];
    const binder = [...new Set(data.contacts.map(row => row.binder))];
    const neighbors = new Map([...target,...binder].map(identifier => [identifier,[]]));
    data.contacts.forEach(row => {neighbors.get(row.target).push(row.binder);neighbors.get(row.binder).push(row.target);});
    for (let sweep=0;sweep<6;sweep++) {
      const binderRank = new Map(binder.map((identifier,index) => [identifier,index]));
      target.sort((first,second) => {
        const score = identifier => neighbors.get(identifier).reduce((total,value)=>total+binderRank.get(value),0)/neighbors.get(identifier).length;
        return score(first)-score(second) || first.localeCompare(second,undefined,{numeric:true});
      });
      const targetRank = new Map(target.map((identifier,index) => [identifier,index]));
      binder.sort((first,second) => {
        const score = identifier => neighbors.get(identifier).reduce((total,value)=>total+targetRank.get(value),0)/neighbors.get(identifier).length;
        return score(first)-score(second) || first.localeCompare(second,undefined,{numeric:true});
      });
    }
    const diagram = byId('contact-diagram'), height = Math.max(230,Math.max(target.length,binder.length)*18+34);
    diagram.setAttribute('viewBox',`0 0 320 ${height}`);
    const positions = new Map();
    [target,binder].forEach(group => group.forEach((identifier,index) => positions.set(identifier,25+index*(height-45)/Math.max(1,group.length-1))));
    diagram.append(svgElement('text',{x:74,y:11,'text-anchor':'end',fill:'#758d96','font-size':10},'TARGET'));
    diagram.append(svgElement('text',{x:246,y:11,fill:'#758d96','font-size':10},'BINDER'));
    data.contacts.forEach(contact => {
      const start = positions.get(contact.target), end = positions.get(contact.binder);
      const path = `M84 ${start} C142 ${start},178 ${end},236 ${end}`;
      const group = svgElement('g',{class:'edge'});
      const edge = svgElement('path',{d:path,stroke:interactionColor(contact.kind),'stroke-width':contact.kind==='Contact'?1:1.7,opacity:contact.kind==='Contact'?.45:.86,fill:'none',class:'network-edge','data-target':contact.target,'data-binder':contact.binder});
      group.append(edge,svgElement('path',{d:path,class:'edge-hit'}));
      group.append(svgElement('title',{},contactText(contact)));
      bindResidues(group,[contact.target,contact.binder],contactText(contact)); diagram.append(group);
    });
    const hotspotIds = new Set(data.hotspots.map(row=>row.id).filter(Boolean));
    [[target,84,74,'end',palette.target],[binder,236,246,'start',palette.binder]].forEach(([numbers,horizontal,textX,anchor,color]) => {
      numbers.forEach(identifier => {
        const group = svgElement('g',{'data-residue':identifier});
        group.append(svgElement('circle',{cx:horizontal,cy:positions.get(identifier),r:3,fill:hotspotIds.has(identifier)?palette.hotspot:color,class:'node-dot'}));
        group.append(svgElement('text',{x:textX,y:positions.get(identifier)+3.5,'text-anchor':anchor,class:'node-label',fill:hotspotIds.has(identifier)?palette.hotspot:color},label(identifier)));
        bindResidues(group,[identifier],descriptive(identifier)); diagram.append(group);
      });
    });
    if (!data.contacts.length) diagram.append(svgElement('text',{x:20,y:50,fill:'#768a90'},'No recorded contacts'));
  }
  function buildBars(container, side, maximumRows=8) {
    const heading = element('div','bar-heading',side==='A'?'Target':'Binder'); container.append(heading);
    const rows = data.residues[side].filter(row=>finite(row.burial)&&row.burial>0).sort((first,second)=>second.burial-first.burial).slice(0,maximumRows);
    const maximum = Math.max(1,...rows.map(row=>row.burial));
    rows.forEach(residue => {
      const button = element('button','bar-row'); button.type='button'; button.dataset.residue=residue.id;
      const background = element('span','bar-bg'), fill = element('span','bar-fill');
      fill.style.width=`${100*residue.burial/maximum}%`; background.append(fill);
      button.append(element('span','',residue.label),background,element('span','bar-value',rounded(residue.burial)));
      bindResidues(button,[residue.id],`${descriptive(residue.id)} · ${rounded(residue.burial)} Å² buried surface area, not binding energy`);
      container.append(button);
    });
    if (!rows.length) container.append(element('p','muted','No per-residue burial recorded.'));
  }
  function buildContactTable() {
    const table = element('table','contact-table'), head = element('thead'), heading = element('tr');
    ['Target','Binder','Type','Å','Pairs'].forEach(text=>heading.append(element('th','',text)));head.append(heading);table.append(head);
    const body = element('tbody');
    [...data.contacts].sort((first,second)=>second.pairs-first.pairs).slice(0,15).forEach(contact=>{
      const row = element('tr'); row.dataset.pair=`${contact.target}|${contact.binder}`; row.tabIndex=0;
      [label(contact.target),label(contact.binder),contact.kind,rounded(contact.distance),String(contact.pairs)].forEach(text=>row.append(element('td','',text)));
      bindResidues(row,[contact.target,contact.binder],contactText(contact));
      row.addEventListener('keydown',event=>{if(event.key==='Enter')row.click();});body.append(row);
    });table.append(body);byId('contact-table').append(table);
    byId('contact-table').append(element('p','muted',`${data.contacts.length} total residue pairs; largest 15 shown here. Complete measurements remain in the CSV table.`));
  }

  function trackCell(residue, text, extraClass='') {
    const cell = element('button','track-cell '+extraClass,text);cell.type='button';cell.dataset.residue=residue.id;
    const message = `${descriptive(residue.id)} · ${residue.interface?'interface':'not a recorded contact'} · pLDDT ${rounded(residue.confidence)} · burial ${rounded(residue.burial)} Å²`;
    cell.title=message;cell.setAttribute('aria-label',message);bindResidues(cell,[residue.id],message);return cell;
  }
  function buildTracks() {
    ['B','A'].forEach(side=>{
      const group = element('div','track-group');group.append(element('div','track-heading',side==='B'?'BINDER':'TARGET · SAVED NUMBERING'));
      const fields = side==='B'?['Sequence','Structure','pLDDT','Burial','Flags']:['Sequence','Burial','Hotspots'];
      const hotspots = new Set(data.hotspots.map(row=>row.id));
      const patch = new Set(data.patches.flatMap(row=>row.res.map(number=>`B:${number}`)));
      const liabilities = new Set(data.liabilities.filter(row=>row.pos).map(row=>`B:${row.pos}`));
      const maximum = Math.max(1,...data.residues[side].map(row=>row.burial||0));
      fields.forEach(field=>{
        const track=element('div','track');track.append(element('span','track-name',field));
        data.residues[side].forEach(residue=>{
          const text = field==='Sequence'?residue.aa:field==='Structure'?(residue.ss||'·'):field==='Flags'?(liabilities.has(residue.id)?'▲':patch.has(residue.id)?'~':''):field==='Hotspots'?(hotspots.has(residue.id)?'●':''):'';
          const cell=trackCell(residue,text,field==='Sequence'&&residue.interface?'iface':'');
          if(field==='pLDDT'){const shade=element('span','shade');shade.style.background=confidenceColor(residue.confidence);cell.append(shade);}
          if(field==='Burial'&&finite(residue.burial)){const bar=element('span','bar');bar.style.height=`${Math.max(1,12*residue.burial/maximum)}px`;cell.append(bar);}
          track.append(cell);
        });group.append(track);
      });byId('sequence-tracks').append(group);
    });
  }
  function finding(title,description,identifiers,warning=false) {
    const button=element('button','finding'+(warning?' warn':''));button.type='button';
    button.append(element('strong','',title),element('span','',description));
    bindResidues(button,identifiers,`${title} · ${description}`);byId('findings').append(button);
  }
  function buildFindings() {
    const missed=data.hotspots.filter(row=>row.status==='missed'),unmapped=data.hotspots.filter(row=>row.status==='unmapped');
    if(missed.length)finding(`${missed.length} intended hotspot${missed.length===1?'':'s'} not contacted`,missed.map(row=>`${row.requested}: ${rounded(row.minimum_distance_angstrom)} Å`).join(' · ')+'. Recorded 4.5 Å heavy-atom contact rule; original admission is unchanged.',missed.map(row=>row.id).filter(Boolean),true);
    if(unmapped.length)finding(`${unmapped.length} hotspot identities cannot be mapped`,'No residue-number offset is guessed. Unmapped is not passed.',unmapped.map(row=>row.id).filter(Boolean),true);
    if(!data.hotspots.length)finding('No intended hotspot assessment',`Declaration status: ${data.hotspot_status.replaceAll('_',' ')}. A missing declaration is not an all-clear.`,[]);
    if(data.polar.length)finding(`${data.polar.length} buried unpaired polar atom${data.polar.length===1?'':'s'}`,data.polar.map(row=>`${label(row.id)} ${row.atom}`).join(' · ')+'. Geometric screen, not an experimental measurement.',data.polar.map(row=>row.id),true);
    data.liabilities.forEach(row=>finding(row.kind,`${row.motif||''}${row.pos?' at '+label('B:'+row.pos):''}. ${row.reason||''}`,row.pos?['B:'+row.pos]:[],row.sev==='warn'));
    data.patches.forEach((row,index)=>finding(`Exposed hydrophobic patch ${index+1}: ${rounded(row.area)} Å²`,'Exposed patch, separate from buried interface surface. No universal pass threshold is applied.',row.res.map(number=>'B:'+number)));
    if(data.reference.status!=='assessed')finding('Reference context is limited',data.reference.reason||'No verified original target conformation is available.',[]);
    else finding('Original target comparison',`Mapped target Cα RMSD ${rounded(data.reference.target_ca_rmsd_angstrom,3)} Å; not a full-protein environment or specificity assessment.`,[]);
  }
  function buildConfidence() {
    const container=byId('confidence-panel'),values=data.residues.B.filter(row=>finite(row.confidence));
    if(!data.confidence_confirmed||!values.length){container.append(element('p','muted','No confirmed per-residue binder pLDDT. B-factors are not reinterpreted from their numeric range. Geometry-based measurements remain available.'));return;}
    const mean=values.reduce((total,row)=>total+row.confidence,0)/values.length,minimum=Math.min(...values.map(row=>row.confidence));
    const stats=element('div','confidence-stats');
    [[mean,'Mean pLDDT'],[minimum,'Lowest pLDDT']].forEach(([value,text])=>{const stat=element('div','confidence-stat');stat.append(element('strong','',rounded(value)),element('span','',text));stats.append(stat);});container.append(stats);
    const control=element('label','threshold','Highlight binder residues below '),input=element('input');input.type='number';input.min=0;input.max=100;input.value=75;input.setAttribute('aria-label','Exploratory confidence threshold');control.append(input);container.append(control);
    const plot=svgElement('svg',{viewBox:'0 0 330 145',class:'mini-plot',role:'img','aria-label':'Per-residue binder pLDDT'});
    [50,70,90].forEach(value=>{const vertical=130-value*1.15;plot.append(svgElement('line',{x1:25,x2:325,y1:vertical,y2:vertical,stroke:'#e8efec'}),svgElement('text',{x:0,y:vertical+3,fill:'#809298','font-size':9},String(value)));});
    const points=values.map((row,index)=>`${25+index*300/Math.max(1,values.length-1)},${130-row.confidence*1.15}`).join(' ');
    plot.append(svgElement('polyline',{points,fill:'none',stroke:palette.target,'stroke-width':1.7}));
    values.forEach((row,index)=>{const dot=svgElement('circle',{cx:25+index*300/Math.max(1,values.length-1),cy:130-row.confidence*1.15,r:3,fill:confidenceColor(row.confidence)});bindResidues(dot,[row.id],descriptive(row.id)+` · pLDDT ${rounded(row.confidence)}`);plot.append(dot);});container.append(plot);
    const list=element('div','footprint-residues'),note=element('p','muted');container.append(list,note);
    function update(){const threshold=Number(input.value);list.replaceChildren();values.filter(row=>row.confidence<threshold).forEach(row=>{const button=element('button','',row.label);bindResidues(button,[row.id],descriptive(row.id)+` · pLDDT ${rounded(row.confidence)}`);list.append(button);});note.textContent=`${values.filter(row=>row.confidence<threshold).length} binder residues below ${threshold}. Exploratory display threshold, not the campaign's admission filter. Target B-factors remain uninterpreted.`;}
    input.addEventListener('input',update);update();
  }
  function buildFootprint() {
    buildBars(byId('target-footprint'),'A',15);
    const list=element('div','footprint-residues');data.residues.A.filter(row=>row.interface).forEach(row=>{const button=element('button','',row.label);button.dataset.residue=row.id;bindResidues(button,[row.id],descriptive(row.id)+` · ${rounded(row.burial)} Å² buried`);list.append(button);});byId('target-footprint').append(list);
    byId('show-footprint').addEventListener('click',()=>{state.color='burial';byId('show-surface').checked=true;document.querySelectorAll('[data-color]').forEach(button=>button.setAttribute('aria-pressed',String(button.dataset.color===state.color)));drawStructure();});
  }

  function drawStructure() {
    if(!state.viewer||!state.model)return;
    clearHover();state.viewer.removeAllShapes();state.viewer.removeAllLabels();state.viewer.removeAllSurfaces();state.surfaceEpoch++;
    const maximum=Math.max(1,...Object.values(data.residues).flat().map(row=>row.burial||0));
    const colorAtom=atom=>{
      const residue=lookup.get(`${atom.chain}:${atom.resi}`);
      if(!residue)return '#b1bdc2';
      if(state.color==='confidence')return residue.side==='B'?confidenceColor(residue.confidence):palette.target;
      if(state.color==='interface')return residue.interface?'#8068ae':residue.side==='A'?'#b4c6d4':'#acd1c1';
      if(state.color==='burial')return burialColor(residue.burial,maximum);
      return residue.side==='A'?palette.target:palette.binder;
    };
    state.model.setStyle({},{cartoon:{colorfunc:colorAtom,arrows:false,thickness:.3}});
    if(data.extra_chains)state.model.setStyle({hetflag:true,not:{chain:['A','B']}},{stick:{radius:.15}});
    data.interactions.forEach(interaction=>{
      if(interaction.kind==='H-bond'?!byId('show-hbonds').checked:!byId('show-salt').checked)return;
      state.viewer.addCylinder({start:point(interaction.start),end:point(interaction.end),radius:.085,dashed:true,color:interaction.kind==='H-bond'?palette.hbond:palette.salt,fromCap:1,toCap:1});
    });
    if(byId('show-hotspots').checked)data.hotspots.filter(row=>row.id&&lookup.has(row.id)).forEach(row=>{const residue=lookup.get(row.id),center=point(residue.xyz);state.viewer.addSphere({center,radius:.8,color:palette.hotspot,alpha:.75});state.viewer.addLabel(residue.label,{position:{x:center.x,y:center.y+2,z:center.z},fontSize:10,fontColor:'#965a22',backgroundOpacity:.8,backgroundColor:'#fff8e9',borderThickness:0,inFront:true});});
    if(byId('show-polar').checked)data.polar.forEach(row=>{const atom=state.viewer.selectedAtoms({chain:row.id.split(':')[0],resi:Number(row.id.split(':')[1]),atom:row.atom})[0];if(atom)state.viewer.addSphere({center:atom,radius:.5,color:'#c67368',alpha:.85});});
    if(byId('show-surface').checked){const epoch=state.surfaceEpoch;Promise.resolve(state.viewer.addSurface($3Dmol.SurfaceType.VDW,{opacity:.52,colorfunc:colorAtom},{chain:'A'})).then(()=>{if(epoch===state.surfaceEpoch)state.viewer.render();}).catch(()=>{byId('viewer-error').textContent='Target surface unavailable; cartoon and interface diagram remain available.';});}
    highlight(state.selected,state.title);state.viewer.render();
  }
  function startViewer() {
    try {
      state.viewer=$3Dmol.createViewer(byId('viewer'),{backgroundColor:'#fbfcfd',antialias:true});
      state.model=state.viewer.addModel(data.pdb,'pdb');
      if(!state.viewer.selectedAtoms({atom:'CA'}).length)throw new Error('No protein Cα coordinates were rendered.');
      drawStructure();state.viewer.zoomTo({chain:['A','B']});state.viewer.render();
      state.viewer.setHoverable({},true,atom=>{const identifier=`${atom.chain}:${atom.resi}`;if(lookup.has(identifier))highlight([identifier],descriptive(identifier));},()=>highlight(state.selected,state.title));
      state.viewer.setClickable({},true,atom=>{const identifier=`${atom.chain}:${atom.resi}`;if(lookup.has(identifier)){state.selected=[identifier];state.title=descriptive(identifier);highlight(state.selected,state.title,true);}});
      new ResizeObserver(()=>{state.viewer.resize();state.viewer.render();}).observe(byId('viewer'));
      byId('viewer').dataset.ready='true';
    } catch(error) {byId('viewer-error').textContent='3D rendering unavailable: '+error.message+'. The recorded 2D diagrams and measurements remain usable.';}
  }
  document.querySelectorAll('[data-color]').forEach(button=>button.addEventListener('click',()=>{state.color=button.dataset.color;document.querySelectorAll('[data-color]').forEach(other=>other.setAttribute('aria-pressed',String(other===button)));drawStructure();}));
  document.querySelectorAll('.layer-controls input').forEach(input=>input.addEventListener('change',drawStructure));
  document.querySelectorAll('[data-tab]').forEach(button=>button.addEventListener('click',()=>{document.querySelectorAll('[data-tab]').forEach(other=>{const active=other===button;other.setAttribute('aria-selected',String(active));byId('panel-'+other.dataset.tab).hidden=!active;});}));
  byId('reset-view').addEventListener('click',()=>{state.selected=[];state.title='';highlight([],state.title);if(state.viewer){state.viewer.zoomTo({chain:['A','B']},450);state.viewer.render();}});
  function download(name,href){const anchor=element('a');anchor.download=name;anchor.href=href;anchor.click();}
  byId('save-image').addEventListener('click',()=>{if(state.viewer)download('binderscope-complex.png',state.viewer.pngURI());});
  byId('save-diagram').addEventListener('click',()=>{
    const copy=byId('contact-diagram').cloneNode(true);
    copy.setAttribute('xmlns','http://www.w3.org/2000/svg');copy.querySelectorAll('.edge-hit').forEach(node=>node.remove());
    copy.querySelectorAll('text').forEach(node=>{node.setAttribute('font-family','monospace');node.setAttribute('font-size','10');});
    const url=URL.createObjectURL(new Blob([new XMLSerializer().serializeToString(copy)],{type:'image/svg+xml'}));download('binderscope-residue-contacts.svg',url);setTimeout(()=>URL.revokeObjectURL(url),1000);
  });
  buildNetwork();buildBars(byId('burial-bars'),'B');buildBars(byId('burial-bars'),'A');buildContactTable();buildTracks();buildFindings();buildConfidence();buildFootprint();startViewer();
})();
