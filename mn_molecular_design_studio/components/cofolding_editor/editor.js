(() => {
  const root=document.getElementById('root');
  let state=null,cap={},previews={},yamlCatalog=[],yamlError='',examples={},menu=null,dragged=null,timer=null,height=0;
  const ids='ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz';
  const labels={protein:'Protein',rna:'RNA',dna:'DNA',ligand:'Ligand'};
  const fields={samples:'Number of Samples',optimization_mode:'Claude-inspired inference optimization',use_yaml:'Use YAML file input',num_recycles:'Number of Recycles',step_scale:'Step Scale',seed_count:'Number of Seeds',seed:'Seed (0 for random)',output_format:'Output Format',boltz_version:'Boltz Version',method:'Method',template_distance_threshold:'Template Distance Threshold',use_potentials:'Use Potentials',use_msa:'Use Multiple Sequence Alignment',max_msa_seqs:'Max MSA Sequences',generate_ipsae:'Generate ipSAE metrics',sampling_steps:'Diffusion Steps',max_parallel_samples:'Parallel Samples',prioritize:'Prioritize in scheduler queue'};
  const fieldName=path=>(cap.prediction?.fields||[]).find(field=>field.key===path.split('.').at(-1))?.label||fields[path.split('.').at(-1)]||path.split('.').at(-1).replaceAll('_',' ');
  const esc=x=>String(x??'').replace(/[&<>"']/g,ch=>({'&':'&amp;','<':'&lt;','>':'&gt;','"':'&quot;',"'":'&#39;'}[ch]));
  const icon=name=>({trash:'<svg viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="1.7"><path d="M4 7h16M9 7V4h6v3M7 7l1 14h8l1-14M10 10v7m4-7v7"/></svg>',copy:'<svg viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="1.6"><path d="M8 8h12v13H8zM16 8V3H3v13h5"/></svg>',file:'<svg viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="1.7"><path d="M6 2h8l4 4v16H6zM14 2v5h4M12 18v-7m-3 3 3-3 3 3"/></svg>'}[name]);
  const post=(type,fields={})=>parent.postMessage({isStreamlitMessage:true,type,...fields},'*');
  const fit=()=>requestAnimationFrame(()=>{const next=Math.ceil(root.getBoundingClientRect().height+4);if(next!==height){height=next;post('streamlit:setFrameHeight',{height:next});}});
  const notify=(immediate=false)=>{clearTimeout(timer);const send=()=>post('streamlit:setComponentValue',{value:state});if(immediate)send();else timer=setTimeout(send,180);};
  const clone=x=>JSON.parse(JSON.stringify(x));
  const normalizePrediction=()=>{for(const field of cap.prediction?.fields||[]){const value=state.settings[field.key];if(value==null||value===''&&['int','float'].includes(field.type||'int'))state.settings[field.key]=field.default;}};
  const fresh=()=>({id:Date.now()+Math.floor(Math.random()*1000),kind:'protein',value:'',copies:1,modifications:[]});
  const kindOf=card=>card.kind==='ligand'?(card.mode||(cap.ccd?'ccd':'smiles')):card.kind;
  const catalog=()=>{if(state.settings?.use_yaml&&yamlCatalog.length)return yamlCatalog;let rows=[];for(const kind of ['protein','dna','rna','ccd','smiles'])state.cards.forEach((card,index)=>{if(kindOf(card)!==kind)return;for(let copy=0;copy<card.copies;copy++)rows.push({id:ids[rows.length],kind,card:index,value:card.value||''});});return rows.filter(row=>row.id);};
  const hasLigand=()=>state.settings?.use_yaml&&cap.boltz?yamlCatalog.some(row=>row.kind==='ligand'):catalog().some(row=>['ccd','smiles'].includes(row.kind));
  const compoundIdentity=()=>`<section class="compound-identity"><h3>Compound identification</h3><div class="identity-fields"><label>Compound ID<input type="text" aria-label="Compound ID" data-path="compound_id" value="${esc(state.compound_id||'')}" placeholder="Optional; defaults to compound name"></label><label>Compound name *<input type="text" aria-label="Compound name" data-path="compound_name" value="${esc(state.compound_name||'')}" placeholder="Name shown in result tables"></label></div><p>The name is required and is saved with the ligand input. The ID defaults to the name when left blank.</p></section>`;
  const chainOptions=(chosen,filter=()=>true)=>catalog().filter(filter).map(row=>`<option value="${row.id}" ${chosen===row.id?'selected':''}>Chain ${row.id}</option>`).join('');
  const control=(path,value,type='number',extra='')=>`<input type="${type}" aria-label="${esc(fieldName(path))}" data-path="${path}" value="${esc(value)}" ${extra}>`;
  const select=(path,value,choices)=>`<select aria-label="${esc(fieldName(path))}" data-path="${path}">${choices.map(([id,label])=>`<option value="${esc(id)}" ${value===id?'selected':''}>${esc(label)}</option>`).join('')}</select>`;
  const check=(path,value,label)=>`<label class="check"><input type="checkbox" data-path="${path}" ${value?'checked':''}>${label}</label>`;
  const handle=(index,restraint=false)=>`<div class="handle"><span class="drag" draggable="true" data-drag="${restraint?'r':'c'}-${index}" title="Drag to reorder" aria-label="Drag to reorder">⠿</span><button class="duplicate" data-action="${restraint?'duplicate-restraint':'duplicate'}" data-index="${index}" title="Duplicate" aria-label="Duplicate ${restraint?'restraint':'chain'}">${icon('copy')}</button></div>`;
  const getPath=path=>path.split('.').reduce((value,key)=>value?.[key],state);
  const setPath=(path,value)=>{const keys=path.split('.');let target=state;for(const key of keys.slice(0,-1))target=target[key];target[keys.at(-1)]=value;};
  const upload=(kind,description,accept,multiple=true)=>`<label class="upload" role="button" tabindex="0" aria-label="Upload ${kind.replaceAll('_',' ')}" data-upload="${kind}"><span>Drag & Drop or <strong>Browse</strong> ${description}</span><input class="hidden" type="file" data-file="${kind}" accept="${accept}" ${multiple?'multiple':''}></label>`;
  function files(kind,items){return items.map((item,index)=>`<div class="file-row"><span class="filename">${esc(item.name)}</span>${kind==='msa_files'?`<select aria-label="Protein chain for ${esc(item.name)}" data-path="msa_files.${index}.chain_id"><option value="">Assign protein chain</option>${chainOptions(item.chain_id,row=>row.kind==='protein')}</select>`:''}${kind==='templates'&&cap.templates?`<select aria-label="Input chain for template ${esc(item.name)}" data-path="templates.${index}.query_chain_id"><option value="">Match protein chains automatically</option>${chainOptions(item.query_chain_id,row=>row.kind==='protein')}</select>`:''}<button class="icon" data-action="remove-file" data-kind="${kind}" data-index="${index}" aria-label="Remove ${esc(item.name)}">${icon('trash')}</button></div>`).join('');}
  function molecule(card,index){
    const chainRows=catalog().filter(row=>row.card===index),chainLabel=chainRows.length>1?`Chains ${chainRows[0]?.id}–${chainRows.at(-1)?.id}`:`Chain ${chainRows[0]?.id||'—'}`;
    const polymer=card.kind!=='ligand',sequence=String(card.value||'').replace(/\s/g,'');
    const type=`<div class="type-picker"><button class="type-button" data-action="type-menu" data-index="${index}" aria-haspopup="listbox" aria-expanded="${menu===index}"><span class="dot ${card.kind}"></span>${labels[card.kind]}<span class="chevron">▾</span></button>${menu===index?`<div class="type-menu" role="listbox">${cap.types.map(kind=>`<button role="option" aria-selected="${card.kind===kind}" class="${card.kind===kind?'selected':''}" data-action="set-type" data-index="${index}" data-kind="${kind}"><span class="dot ${kind}"></span>${labels[kind]}</button>`).join('')}</div>`:''}</div>`;
    const toolbar=`<div class="toolbar">${type}<div class="copies">Copies <div class="stepper"><button data-action="copies" data-index="${index}" data-delta="-1" aria-label="Fewer copies" ${card.copies<=1?'disabled':''}>−</button><input type="number" aria-label="Copies of chain ${chainRows[0]?.id||''}" min="1" max="52" data-path="cards.${index}.copies" value="${card.copies}"><button data-action="copies" data-index="${index}" data-delta="1" aria-label="More copies" ${card.copies>=52?'disabled':''}>+</button></div></div><span class="spacer"></span>${polymer?`<button class="pill" data-action="plain" data-index="${index}">${card.plain?'View Sequence':'View Plain Text'}</button>`:''}<span class="chain">${chainLabel}</span><button class="icon" data-action="delete" data-index="${index}" title="Remove chain" aria-label="Remove chain ${chainRows[0]?.id||''}">${icon('trash')}</button></div>`;
    let body='';
    if(polymer){
      let contents='';
      if(!sequence||card.plain||card.editing)contents=`<textarea aria-label="${labels[card.kind]} sequence for chain ${chainRows[0]?.id||''}" placeholder="Paste or type your ${labels[card.kind].toLowerCase()} sequence here…" data-sequence="${index}" data-path="cards.${index}.value">${esc(card.value)}</textarea>`;
      else {let blocks=[];for(let start=0;start<sequence.length;start+=20)blocks.push(`<span class="sequence-block">${esc(sequence.slice(start,start+20))}<sup>${Math.min(start+20,sequence.length)}</sup></span>`);contents=`<div class="styled-sequence" tabindex="0" role="textbox" aria-label="Edit ${labels[card.kind].toLowerCase()} sequence for chain ${chainRows[0]?.id||''}" data-action="edit" data-index="${index}">${blocks.join('')}</div>`;}
      body=`<div class="sequence-area">${contents}<span class="count">${sequence.length?`${sequence.length} ${card.kind==='protein'?'residues':'nucleotides'}`:''}</span></div>`;
      if(card.kind==='protein'&&(cap.cyclic||cap.modifications))body+=`<div class="features">${cap.cyclic?check(`cards.${index}.cyclic`,card.cyclic,'Cyclic'):''}${cap.modifications?`<button class="text-action" data-action="add-modification" data-index="${index}">＋ &nbsp;Add modification</button>`:''}</div>`;
      body+=`<div class="modifications">${(card.modifications||[]).map((mod,j)=>`<div class="modification"><label>Residue position${control(`cards.${index}.modifications.${j}.position`,mod.position,'number','min="1"')}</label><label>CCD code${control(`cards.${index}.modifications.${j}.code`,mod.code,'text','placeholder="MSE"')}</label><button class="icon" data-action="remove-modification" data-index="${index}" data-mod="${j}" aria-label="Remove modification">${icon('trash')}</button></div>`).join('')}</div>`;
    }else{
      const mode=card.mode||(cap.ccd?'ccd':'smiles');
      body=`<div class="ligand-help">Ligand SMILES or CCD code (use CCD code for ions)</div><div class="tabs" role="tablist">${(cap.ccd?['ccd']:[]).concat(cap.smiles?['smiles']:[]).map(m=>`<button role="tab" aria-selected="${mode===m}" class="${mode===m?'active':''}" data-action="ligand-mode" data-index="${index}" data-mode="${m}">${m==='ccd'?'Enter CCD Code':'SMILES'}</button>`).join('')}</div><input class="ligand-input" aria-label="${mode==='ccd'?'CCD code':'Ligand SMILES'}" placeholder="${mode==='ccd'?'Enter CCD component code, e.g. ATP':'Paste or type ligand SMILES here…'}" data-path="cards.${index}.value" value="${esc(card.value)}"><div data-preview-slot="${card.id}">${mode==='smiles'&&previews[String(card.id)]?`<img class="ligand-preview" alt="Ligand chemical structure" src="${esc(previews[String(card.id)])}">`:''}</div>`;
    }
    return `<section class="entity" data-card="${index}">${handle(index)}${toolbar}${body}</section>`;
  }
  function restraint(row,index){
    const kind=Object.keys(row)[0],d=row[kind];
    const descriptions={contact:'A contact restraint constrains the distance between a single residue/atom on one chain and a single residue/atom on another chain.',pocket:'A pocket restraint constrains the distance between selected pocket residues and a binder chain.',bond:'A covalent bond connects a specific atom on one chain to a specific atom on another chain.'};
    const endpoint=(which,title,atom=false)=>`<div class="endpoint"><select aria-label="${title}" data-path="restraints.${index}.${kind}.${which}.0">${chainOptions(d[which][0])}</select><input aria-label="${title} residue or atom" placeholder="${atom?'Residue index':'Residue index / atom'}" data-path="restraints.${index}.${kind}.${which}.1" value="${esc(d[which][1])}">${atom?`<input aria-label="${title} atom name" placeholder="Atom name" data-path="restraints.${index}.${kind}.${which}.2" value="${esc(d[which][2])}">`:''}</div>`;
    let body='';
    if(kind==='contact'||kind==='bond')body=`<div class="endpoint-row">${endpoint(kind==='contact'?'token1':'atom1','Chain 1',kind==='bond')}${endpoint(kind==='contact'?'token2':'atom2','Chain 2',kind==='bond')}</div>`;
    else{
      const pocketChain=d.pocket_chain||catalog().find(r=>r.kind==='protein')?.id||'A',pocketSeq=catalog().find(r=>r.id===pocketChain)?.value.replace(/\s/g,'')||'';
      const positions=new Set((d.pocket_residues||'').split(/[\s,]+/).filter(Boolean).map(Number));
      body=`<div class="endpoint-row pocket-fields"><label>Binder chain<select data-path="restraints.${index}.pocket.binder">${chainOptions(d.binder)}</select></label><label>Pocket chain<select data-path="restraints.${index}.pocket.pocket_chain">${chainOptions(pocketChain,r=>!['ccd','smiles'].includes(r.kind))}</select></label><label style="flex:1">Pocket residues (e.g. 45 46 47)<input style="width:100%;display:block" data-path="restraints.${index}.pocket.pocket_residues" value="${esc(d.pocket_residues||'')}"></label></div>${pocketSeq?`<div class="small-note">Click or drag to select pocket residues on chain ${pocketChain}.</div><div class="residue-picker" data-picker="${index}">${[...pocketSeq].map((aa,j)=>`<button title="${pocketChain}${j+1}" aria-label="${pocketChain}${j+1}" data-residue="${j+1}" data-restraint="${index}" class="${positions.has(j+1)?'selected':''}">${esc(aa)}</button>`).join('')}</div><button class="text-action" data-action="select-residues" data-index="${index}" data-all="1">Select all</button> &nbsp;<button class="text-action" data-action="select-residues" data-index="${index}" data-all="0">Clear</button>`:''}`;
    }
    if(kind!=='bond')body+=`<label class="distance">${kind==='pocket'?'Max distance':'Distance'} (Å)${control(`restraints.${index}.${kind}.max_distance`,d.max_distance,'number','min="4" max="20" step="0.5" style="width:166px;display:block"')}</label>${check(`restraints.${index}.${kind}.force`,d.force,'Force restraint with potentials')}`;
    return `<section class="restraint" data-restraint-card="${index}" style="${restraintError(row)?'border-color:#ef5555':''}">${handle(index,true)}<div class="toolbar"><select class="restraint-type" aria-label="Restraint type" data-restraint-type="${index}">${[['contact','Contact'],['pocket','Pocket'],['bond','Covalent Bond']].map(([id,label])=>`<option value="${id}" ${kind===id?'selected':''}>${label}</option>`).join('')}</select><span class="spacer"></span><button class="icon" data-action="delete-restraint" data-index="${index}" aria-label="Remove restraint">${icon('trash')}</button></div><p class="restraint-description">${descriptions[kind]}</p>${body}<div class="error" data-restraint-error="${index}">${esc(restraintError(row))}</div></section>`;
  }
  function setting(title,description,content,toggle=false){return `<section class="setting ${toggle?'toggle-setting':''}"><h3 class="setting-title">${title}</h3><p class="setting-description">${description}</p>${content}</section>`;}
  const toggle=(key)=>`<input class="toggle" type="checkbox" role="switch" aria-label="${fields[key]||key}" data-path="settings.${key}" ${state.settings[key]?'checked':''}>`;
  function boltzOptions(){
    const s=state.settings;
    return `<label class="samples">Number of Samples *${control('settings.samples',s.samples,'number','min="1" max="5"')}</label><details class="optional" ${state.optional_open?'open':''}><summary>Optional settings</summary>
    ${setting('Claude-inspired inference optimization','Choose the optimized inference profile for the local Boltz image.',select('settings.optimization_mode',s.optimization_mode,(cap.modes||[]).map(mode=>[mode,mode==='off'?'None':mode[0].toUpperCase()+mode.slice(1)])))}
    ${setting('Use YAML file input','Upload a Boltz YAML configuration file instead of filling out the molecule cards. Missing protein MSAs are generated automatically when MSA is enabled.',toggle('use_yaml')+(s.use_yaml?upload('yaml','your yaml file','.yaml,.yml',false)+files('yaml_files',state.yaml_text?[{name:state.yaml_name||'Boltz YAML input'}]:[])+upload('yaml_files','files referenced by the YAML','.a3m,.txt,.pdb,.cif,.mmcif,.ent')+files('yaml_files',state.yaml_files):''),true)}
    ${setting('Number of Recycles','Number of refinement rounds. Higher values may improve accuracy but take longer.',control('settings.num_recycles',s.num_recycles,'number','min="1" max="100"'))}
    ${setting('Step Scale','Lower values increase diversity; recommended between 1 and 2.',control('settings.step_scale',s.step_scale,'number','min="0.1" max="3" step="0.1"'))}
    ${setting('Number of Seeds','Number of prediction runs. Seed 0 chooses a random starting seed.',control('settings.seed_count',s.seed_count,'number','min="1" max="16"'))}
    ${setting('Seed (0 for random)','Set a nonzero seed for reproducibility. The default, 0, chooses a fresh random seed.',`<div class="seed-controls">${control('settings.seed',s.seed,'number','min="0" max="2147483647"')}<button class="text-action" data-action="randomize">Randomize</button></div>`)}
    ${setting('Output Format','File format used for predicted structures.',select('settings.output_format',s.output_format,[['pdb','PDB'],['mmcif','mmCIF']]))}
    ${setting('Boltz Version','Model version installed in the selected Docker image.',select('settings.boltz_version',s.boltz_version,[['2.2.1','2.2.1']]))}
    ${setting('Method','Experimental method bias for the prediction (optional).',select('settings.method',s.method,['','md','x-ray diffraction','electron microscopy','solution nmr','solid-state nmr','neutron diffraction','theoretical model','afdb','boltz-1','other'].map(value=>[value,value||'None'])))}
    ${!cap.settings_only?setting('Templates','Upload template CIF or PDB files to guide the structure prediction.',upload('templates','your cif/pdb files','.cif,.mmcif,.pdb,.ent')+files('templates',state.templates)):''}
    ${!cap.settings_only?setting('Template Distance Threshold','Maximum deviation from a template in Å. A value enables template restraints with potentials; leave blank for an unforced template.',control('settings.template_distance_threshold',s.template_distance_threshold,'number','min="0.1" step="0.5"')):''}
    ${setting('Use Potentials','Apply inference time potentials to improve the physical quality of predicted poses.',toggle('use_potentials'),true)}
    ${setting('Use Multiple Sequence Alignment','Uses uploaded A3Ms first, then matching A3Ms in the shared reference files; remaining protein chains use the ColabFold server.',toggle('use_msa'),true)}
    ${!cap.settings_only?setting('Custom A3M Files','Upload custom A3M alignments and assign each file to its protein chain. MSA generation is skipped for those chains.',upload('msa_files','your a3m files','.a3m,.txt')+files('msa_files',state.msa_files)):''}
    ${setting('Max MSA Sequences','Limit the sequences used for each MSA. Leave blank to use the engine default (8192).',control('settings.max_msa_seqs',s.max_msa_seqs,'number','min="1" max="8192"'))}
    ${setting('Generate ipSAE metrics','Calculate ipSAE for protein–protein interfaces after structure prediction.',toggle('generate_ipsae'),true)}
    ${setting('Diffusion Steps','Number of diffusion sampling steps.',control('settings.sampling_steps',s.sampling_steps,'number','min="1" max="1000"'))}
    ${setting('Parallel Samples','Maximum samples evaluated concurrently. Lower values reduce GPU memory use.',control('settings.max_parallel_samples',s.max_parallel_samples,'number','min="1" max="5"'))}
    ${setting('Prioritize in scheduler queue','Move this job ahead of normal queued work when resources are available.',toggle('prioritize'),true)}
    </details>`;
  }
  function restraintError(row){
    const kind=Object.keys(row)[0],d=row[kind],known=new Set(catalog().map(chain=>chain.id));
    if(kind==='pocket')return !known.has(d.binder)?'Choose a binder chain.':(!d.contacts_input?.trim()?'Select at least one pocket residue.':'');
    for(const key of kind==='contact'?['token1','token2']:['atom1','atom2']){const endpoint=d[key];if(!known.has(endpoint[0]))return 'Choose a chain for each endpoint.';if(!String(endpoint[1]).trim())return `Chain ${endpoint[0]} is missing a residue/atom.`;if(kind==='bond'&&!String(endpoint[2]).trim())return `Chain ${endpoint[0]} is missing an atom name.`;}
    return '';
  }
  function predictionOptions(){
    const config=cap.prediction,s=state.settings;
    const min=config.sample_min_modes?.[s.optimization_mode]||1;
    const max=config.sample_max||5;
    const renderField=field=>{
      const path=`settings.${field.key}`,value=s[field.key]??field.default;
      let input;
      if(field.type==='bool')input=`<input class="toggle" type="checkbox" role="switch" aria-label="${esc(field.label)}" data-path="${path}" ${value?'checked':''}>`;
      else if(field.type==='select')input=select(path,value,field.choices.map(choice=>Array.isArray(choice)?choice:[choice,choice]));
      else if(field.type==='text')input=control(path,value,'text');
      else input=control(path,value,'number',`min="${field.min}" max="${field.max}" step="${field.step||(field.type==='float'?0.1:1)}"`);
      return setting(esc(field.label),esc(field.help||''),input,field.type==='bool');
    };
    const fields=config.fields.filter(field=>!field.required&&!field.hidden);
    const templates=cap.templates?templateControls():'';
    const options=fields.map(field=>(field.key==='use_msa'?templates:'')+renderField(field)).join('')+(!fields.some(field=>field.key==='use_msa')?templates:'');
    return `<label class="samples">Number of Samples *${control('settings.samples',s.samples,'number',`min="${min}" max="${max}"`)}</label><details class="optional" ${state.optional_open?'open':''}><summary>Optional settings</summary>${options}${(config.fixed_settings||[]).map(text=>`<p class="small-note">${esc(text)}</p>`).join('')}</details>${config.fields.filter(field=>field.required&&!field.hidden).map(renderField).join('')}`;
  }
  function templateControls(){return setting('Templates','Upload PDB or mmCIF structures. Protein chains are aligned by sequence; select a chain to restrict the template.',upload('templates','your cif/pdb files','.cif,.mmcif,.pdb,.ent')+files('templates',state.templates));}
  function error(){
    if(cap.prediction){
      const s=state.settings,min=cap.prediction.sample_min_modes?.[s.optimization_mode]||1;
      const max=cap.prediction.sample_max||5;
      if(!Number.isInteger(Number(s.samples))||s.samples<min||s.samples>max)return `Request between ${min} and ${max} samples for this inference mode.`;
      for(const field of cap.prediction.fields){const value=s[field.key]??field.default;if(field.required&&!value)return 'Review and accept the AlphaFold3 weight terms.';if(!['bool','select','text'].includes(field.type)&&(!Number.isFinite(Number(value))||value===''||Number(value)<field.min||Number(value)>field.max))return `Enter a valid ${field.label.toLowerCase()}.`;}
    }
    if(cap.settings_only)return '';
    if(state.settings?.use_yaml){if(!state.yaml_text?.trim())return 'Upload a Boltz YAML file.';if(yamlError)return yamlError;if(hasLigand()&&!state.compound_name?.trim())return 'Enter a compound name for this ligand complex.';return '';}
    if(!state.cards.length||state.cards.some(card=>!String(card.value||'').trim()))return 'All sequences must be non empty.';
    if(state.cards.reduce((sum,card)=>sum+Number(card.copies||1),0)>52)return 'Enter at most 52 chains and ligand copies.';
    if(state.cards.some(card=>(card.modifications||[]).some(mod=>!mod.code?.trim())))return 'Every modification needs a CCD code.';
    if((state.msa_files||[]).some(file=>!file.chain_id)&&state.settings.use_msa)return 'Assign each A3M file to its protein chain.';
    if(hasLigand()&&!state.compound_name?.trim())return 'Enter a compound name for this ligand complex.';
    if(state.restraints.some(row=>restraintError(row)))return 'Complete all restraint endpoints.';
    return '';
  }
  function render(){
    if(!state)return;
    let html='';
    if(cap.examples){const examples=[['7QF4','7QF4 - MiniSOG monomer'],['8JYR','8JYR - Antibody complex']];if(cap.ligand_examples)examples.push(['ligand',`Protein-Ligand Complex${cap.affinity?' + Affinity':''}`]);html=`<div class="caption">EXAMPLES</div><div class="examples">${examples.map(([name,label])=>`<button data-action="example" data-name="${name}">${label}</button>`).join('')}</div>`;}
    if(state.settings.use_yaml&&cap.boltz)html+=`<div class="yaml-notice">YAML input is enabled. ${state.yaml_text?`Using ${esc(state.yaml_name||'uploaded Boltz YAML')}.`:'Upload your file under Optional settings.'}</div>`+(hasLigand()?compoundIdentity():'');
    else if(!cap.settings_only){
      html+=state.cards.map(molecule).join('')+`<div class="actions"><button class="outline" data-action="add" ${state.cards.length>=16?'disabled':''}>＋ &nbsp;Add chain</button><button class="outline" data-action="fasta">${icon('file')} Import FASTA</button><button class="clear" data-action="clear">Clear all</button><input class="hidden" type="file" data-file="fasta" accept=".fa,.fasta,.fna,.txt"></div>`;
      if(hasLigand())html+=compoundIdentity();
      if(cap.restraints)html+=(state.restraints.length?'<div class="caption">RESTRAINTS</div>':'')+state.restraints.map(restraint).join('')+`<div class="restraint-help"><button class="outline" data-action="add-restraint">＋ &nbsp;Add restraint</button><span>Restraints constrain how different chains are arranged in your complex, and can incorporate prior knowledge of structural contacts.</span></div>`;
      const ligands=catalog().filter(row=>['ccd','smiles'].includes(row.kind));
      if(cap.affinity_control&&ligands.length)html+=`<div class="affinity">${check('settings.predict_affinity',state.settings.predict_affinity,'Predict affinity (small molecule ligand)')}${state.settings.predict_affinity?`<label>Binder chain <select data-path="settings.affinity_binder">${chainOptions(state.settings.affinity_binder,row=>['ccd','smiles'].includes(row.kind))}</select></label>`:''}</div>`;
    }
    if(!cap.boltz&&!cap.prediction&&(cap.templates||cap.msa))html+=`<details class="optional" ${state.optional_open?'open':''}><summary>Optional settings</summary>${cap.templates?templateControls():''}${cap.msa?setting('Custom A3M Files','Assign an alignment to its protein chain.',upload('msa_files','your a3m files','.a3m,.txt')+files('msa_files',state.msa_files)):''}</details>`;
    if(cap.boltz||cap.prediction)html+=(cap.boltz?boltzOptions():predictionOptions())+`<div class="submit-area"><button class="submit" data-action="submit" ${error()||!cap.can_submit?'disabled':''}>Submit</button><div class="validation">${esc(error()||(!cap.can_submit?'Select a project to submit.':''))}</div></div>`;
    root.innerHTML=html;
    root.querySelector('.optional')?.addEventListener('toggle',event=>{state.optional_open=event.target.open;notify();fit();});
    root.querySelectorAll('[data-sequence]').forEach(input=>input.addEventListener('blur',event=>{const card=state.cards[Number(input.dataset.sequence)];card.value=card.value.replace(/\s/g,'').toUpperCase();card.editing=false;notify(true);if(!card.plain&&!root.contains(event.relatedTarget))render();}));
    root.querySelectorAll('[data-upload]').forEach(zone=>{zone.addEventListener('dragover',event=>{event.preventDefault();zone.classList.add('drag-over');});zone.addEventListener('dragleave',()=>zone.classList.remove('drag-over'));zone.addEventListener('drop',event=>{event.preventDefault();zone.classList.remove('drag-over');readFiles(zone.dataset.upload,[...event.dataTransfer.files]);});});
    root.querySelectorAll('[data-drag]').forEach(node=>node.addEventListener('dragstart',event=>{dragged=node.dataset.drag;event.dataTransfer.setData('text/plain',dragged);}));
    root.querySelectorAll('[data-card],[data-restraint-card]').forEach(node=>{node.addEventListener('dragover',event=>{event.preventDefault();node.classList.add('drag-over');});node.addEventListener('dragleave',()=>node.classList.remove('drag-over'));node.addEventListener('drop',event=>{event.preventDefault();if(!dragged)return;const [kind,from]=dragged.split('-');const isCard=node.dataset.card!==undefined;if((kind==='c')!==isCard)return;const rows=isCard?state.cards:state.restraints;const to=Number(isCard?node.dataset.card:node.dataset.restraintCard);rows.splice(to,0,rows.splice(Number(from),1)[0]);dragged=null;render();notify(true);});});
    fit();
  }
  function newRestraint(kind='contact'){
    const chainRows=catalog(),a=chainRows[0]?.id||'A',b=chainRows[1]?.id||a;
    return kind==='contact'?{contact:{token1:[a,''],token2:[b,''],max_distance:5,force:false}}:kind==='bond'?{bond:{atom1:[a,'',''],atom2:[b,'','']}}:{pocket:{binder:b,pocket_chain:a,pocket_residues:'',contacts_input:'',max_distance:5,force:false}};
  }
  const updatePocket=index=>{const d=state.restraints[index]?.pocket;if(d)d.contacts_input=(d.pocket_residues||'').split(/[\s,]+/).filter(Boolean).map(residue=>`${d.pocket_chain}, ${residue}`).join('\n');};
  function updateInput(node){
    if(!node.dataset.path)return;
    let value=node.type==='checkbox'?node.checked:node.value;
    if(node.type==='number')value=node.value===''?'':Number(node.value);
    if(node.dataset.path.endsWith('.copies'))value=Math.max(1,Math.min(52,Number(value)||1));
    setPath(node.dataset.path,value);
    const presets=cap.prediction?.parameter_presets;
    if(presets){
      const presetKey=presets.selector_key, presetPath=`settings.${presetKey}`, presetValues=presets.values||{};
      const customLabel=presets.custom_label||'Custom';
      const presetFields=new Set(Object.keys(Object.values(presetValues)[0]||{}));
      if(node.dataset.path===presetPath){
        const values=presetValues[value];
        if(values)for(const [key,presetValue] of Object.entries(values))state.settings[key]=presetValue;
      }else if(presetFields.has(node.dataset.path.split('.').at(-1))){
        const match=Object.entries(presetValues).find(([,values])=>Object.entries(values).every(([key,presetValue])=>Number(state.settings[key])===Number(presetValue)));
        state.settings[presetKey]=match?.[0]||customLabel;
      }
    }
    if(node.dataset.path==='settings.predict_affinity'&&value&&!state.settings.affinity_binder)state.settings.affinity_binder=catalog().find(row=>['ccd','smiles'].includes(row.kind))?.id||'';
    const match=node.dataset.path.match(/^restraints\.(\d+)\.pocket\./);if(match)updatePocket(Number(match[1]));
    root.querySelectorAll('[data-restraint-card]').forEach(section=>{const index=Number(section.dataset.restraintCard),message=restraintError(state.restraints[index]);section.style.borderColor=message?'#ef5555':'';section.querySelector('[data-restraint-error]').textContent=message;});
    const validation=root.querySelector('.validation'),submit=root.querySelector('.submit');
    if(validation)validation.textContent=error()||(!cap.can_submit?'Select a project to submit.':'');if(submit)submit.disabled=Boolean(error()||!cap.can_submit);
    notify();
  }
  root.addEventListener('input',event=>updateInput(event.target));
  root.addEventListener('change',event=>{
    const node=event.target;
    if(node.dataset.file){readFiles(node.dataset.file,[...node.files]);return;}
    if(node.dataset.restraintType!==undefined){state.restraints[Number(node.dataset.restraintType)]=newRestraint(node.value);render();notify(true);return;}
    updateInput(node);
    const ligandPath=node.dataset.path?.match(/^cards\.(\d+)\.value$/);if(ligandPath&&state.cards[Number(ligandPath[1])].kind==='ligand'){state.cards[Number(ligandPath[1])].preview_value=node.value;notify(true);}
    if(node.type==='checkbox'||node.tagName==='SELECT'||node.dataset.path?.endsWith('.copies')||node.dataset.path?.endsWith('.pocket_residues')){render();notify(true);}
  });
  root.addEventListener('click',async event=>{
    const node=event.target.closest('[data-action]');if(!node)return;
    const action=node.dataset.action,index=Number(node.dataset.index),card=state.cards[index];
    if(action==='type-menu'){menu=menu===index?null:index;render();return;}
    if(action==='set-type'){card.kind=node.dataset.kind;card.value='';card.modifications=[];card.cyclic=false;card.mode=cap.ccd?'ccd':'smiles';menu=null;}
    if(action==='add')state.cards.push(fresh());
    if(action==='delete'){state.cards.splice(index,1);if(!state.cards.length)state.cards.push(fresh());}
    if(action==='duplicate'&&state.cards.length<16)state.cards.splice(index+1,0,{...clone(card),id:Date.now()});
    if(action==='copies')card.copies=Math.max(1,Math.min(52,card.copies+Number(node.dataset.delta)));
    if(action==='plain')card.plain=!card.plain;
    if(action==='edit'){card.editing=true;render();root.querySelector(`[data-sequence="${index}"]`)?.focus();return;}
    if(action==='add-modification')(card.modifications||=([])).push({position:1,code:''});
    if(action==='remove-modification')card.modifications.splice(Number(node.dataset.mod),1);
    if(action==='ligand-mode'){card.mode_values||={};card.mode_values[card.mode||(cap.ccd?'ccd':'smiles')]=card.value;card.mode=node.dataset.mode;card.value=card.mode_values[card.mode]||'';}
    if(action==='fasta'){root.querySelector('[data-file="fasta"]').click();return;}
    if(action==='clear'){state.cards=[fresh()];state.restraints=[];state.msa_files=[];state.settings.predict_affinity=false;}
    if(action==='add-restraint')state.restraints.push(newRestraint());
    if(action==='delete-restraint')state.restraints.splice(index,1);
    if(action==='duplicate-restraint')state.restraints.splice(index+1,0,clone(state.restraints[index]));
    if(action==='select-residues'){const d=state.restraints[index].pocket;const sequence=catalog().find(row=>row.id===d.pocket_chain)?.value||'';d.pocket_residues=node.dataset.all==='1'?[...sequence.replace(/\s/g,'')].map((_,i)=>i+1).join(' '):'';updatePocket(index);}
    if(action==='randomize'){const values=new Uint32Array(1);crypto.getRandomValues(values);state.settings.seed=values[0]%2147483646+1;}
    if(action==='remove-file'){if(node.dataset.kind==='yaml_files'&&node.closest('.setting')?.querySelector('[data-file="yaml"]')&&node.closest('.file-row').querySelector('.filename').textContent===state.yaml_name){state.yaml_text='';state.yaml_name='';}else state[node.dataset.kind].splice(index,1);}
    if(action==='example'){try{if(!examples[node.dataset.name])examples=await(await fetch('examples.json')).json();const example=examples[node.dataset.name];if(example){state.cards=clone(example.cards).filter(card=>{const kind=card.kind==='ligand'?(card.mode||'smiles'):card.kind;return kind!=='ligand'||cap.ligand_examples;});state.cards.filter(card=>card.kind==='ligand').forEach(card=>card.preview_value=card.value);state.restraints=cap.restraints?clone(example.restraints||[]):[];state.msa_files=[];state.settings.use_yaml=false;state.settings.predict_affinity=Boolean(cap.affinity_control&&example.affinity);state.settings.affinity_binder=state.settings.predict_affinity?(catalog().find(row=>['ccd','smiles'].includes(row.kind))?.id||''):'';}}catch{state.message='Could not load the example.';}}
    if(action==='submit'){if(error()||!cap.can_submit)return;state.submit_token=String(Date.now())+'-'+Math.random().toString(36).slice(2);notify(true);return;}
    render();notify(true);
  });
  let selection=null;
  root.addEventListener('pointerdown',event=>{const node=event.target.closest('[data-residue]');if(!node)return;event.preventDefault();selection={index:Number(node.dataset.restraint),start:Number(node.dataset.residue),last:Number(node.dataset.residue)};});
  root.addEventListener('pointerover',event=>{const node=event.target.closest('[data-residue]');if(selection&&node&&Number(node.dataset.restraint)===selection.index)selection.last=Number(node.dataset.residue);});
  window.addEventListener('pointerup',()=>{if(!selection)return;const d=state.restraints[selection.index].pocket,selected=new Set((d.pocket_residues||'').split(/[\s,]+/).filter(Boolean).map(Number));for(let i=Math.min(selection.start,selection.last);i<=Math.max(selection.start,selection.last);i++){if(selection.start===selection.last&&selected.has(i))selected.delete(i);else selected.add(i);}d.pocket_residues=[...selected].sort((a,b)=>a-b).join(' ');updatePocket(selection.index);selection=null;render();notify(true);});
  root.addEventListener('keydown',event=>{if((event.key==='Enter'||event.key===' ')&&event.target.dataset.upload){event.preventDefault();event.target.querySelector('input[type=file]').click();}if(event.target.dataset.action==='edit'&&(event.key==='Enter'||event.key===' ')){event.preventDefault();event.target.click();}if(event.target.dataset.residue&&(event.key==='Enter'||event.key===' ')){event.preventDefault();const index=Number(event.target.dataset.restraint),position=Number(event.target.dataset.residue),d=state.restraints[index].pocket,positions=new Set(d.pocket_residues.split(/[\s,]+/).filter(Boolean).map(Number));if(positions.has(position))positions.delete(position);else positions.add(position);d.pocket_residues=[...positions].sort((a,b)=>a-b).join(' ');updatePocket(index);render();notify(true);}if(event.key==='Escape'&&menu!==null){menu=null;render();}});
  async function readFiles(kind,fileList){
    for(const file of fileList){
      if(file.size>20*1024*1024){state.message='Each input file must be smaller than 20 MB.';continue;}
      const text=await file.text();
      if(kind==='fasta'){const records=text.split(/^>/m).map(record=>{const lines=record.trim().split(/\r?\n/);return record.startsWith('>')?'':lines.slice(1).join('').replace(/\s/g,'');}).filter(Boolean);if(!text.trim().startsWith('>')){state.message='Import a FASTA file beginning with a > header.';continue;}if(!records.length||records.length>16){state.message='Import 1 to 16 FASTA records.';continue;}state.cards=records.map(value=>({...fresh(),value:value.toUpperCase()}));state.restraints=[];state.msa_files=[];}
      else if(kind==='yaml'){state.yaml_text=text;state.yaml_name=file.name;}
      else{const item={name:file.name,text};if(kind==='msa_files'){const query=text.split(/^>/m)[1]?.split(/\r?\n/).slice(1).join('').replace(/[a-z.\-\s]/g,'')||'';item.chain_id=catalog().find(row=>row.kind==='protein'&&row.value.replace(/\s/g,'').toUpperCase()===query)?.id||'';}state[kind]=state[kind].filter(existing=>existing.name!==file.name);state[kind].push(item);}
    }
    render();if(state.message){const p=document.createElement('p');p.className='error';p.textContent=state.message;root.append(p);state.message='';fit();}notify(true);
  }
  window.addEventListener('message',event=>{if(event.data.type!=='streamlit:render')return;const capabilitiesChanged=JSON.stringify(cap)!==JSON.stringify(event.data.args.capabilities);cap=event.data.args.capabilities;previews=event.data.args.previews||{};const changed=JSON.stringify(yamlCatalog)!==JSON.stringify(event.data.args.yaml_chains||[]);yamlCatalog=event.data.args.yaml_chains||[];yamlError=event.data.args.yaml_error||'';if(!state){state=clone(event.data.args.initial);state.settings||={};state.templates||=[];state.msa_files||=[];state.yaml_files||=[];state.compound_name||='';state.compound_id||='';normalizePrediction();render();}else {state.compound_name||='';state.compound_id||='';if(changed||capabilitiesChanged){normalizePrediction();render();return;}const submit=root.querySelector('.submit'),validation=root.querySelector('.validation');if(submit)submit.disabled=Boolean(error()||!cap.can_submit);if(validation)validation.textContent=error()||(!cap.can_submit?'Select a project to submit.':'');root.querySelectorAll('[data-preview-slot]').forEach(slot=>{const source=previews[slot.dataset.previewSlot];slot.innerHTML=source?`<img class="ligand-preview" alt="Ligand chemical structure" src="${esc(source)}">`:'';});fit();}});
  new ResizeObserver(fit).observe(root);
  post('streamlit:componentReady',{apiVersion:1});
})();
