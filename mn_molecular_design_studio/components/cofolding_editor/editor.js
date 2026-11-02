(() => {
  const root=document.getElementById('root');
  let state=null,cap={},previews={},yamlCatalog=[],yamlError='',examples={},menu=null,dragged=null,timer=null,height=0,serverCapSignature='',lastSent='';
  const ids='ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz';
  const labels={protein:'Protein',rna:'RNA',dna:'DNA',ligand:'Ligand',ccd:'Ligand (CCD)',smiles:'Ligand (SMILES)',modification:'Modification',bond:'Bond',crosslink:'Crosslink',contact:'Contact',pocket:'Pocket'};
  const decoration=card=>['modification','bond','crosslink','contact','pocket'].includes(card.kind);
  function applyLocalfoldProfile(){
    const profile=cap.localfold_profiles?.[cap.prediction?state?.settings?.model:cap.localfold_model];
    if(profile){const hasSettings=Boolean(cap.prediction);Object.assign(cap,clone(profile));if(!hasSettings)cap.prediction=null;}
  }
  const fields={samples:'Number of Samples',optimization_mode:'Claude-inspired inference optimization',use_yaml:'Use YAML file input',num_recycles:'Number of Recycles',step_scale:'Step Scale',seed_count:'Number of Seeds',seed:'Seed (0 for random)',output_format:'Output Format',boltz_version:'Boltz Version',method:'Method',template_distance_threshold:'Template Distance Threshold',use_potentials:'Use Potentials',use_msa:'Use Multiple Sequence Alignment',max_msa_seqs:'Max MSA Sequences',generate_ipsae:'Generate ipSAE metrics',sampling_steps:'Diffusion Steps',max_parallel_samples:'Parallel Samples',prioritize:'Prioritize in scheduler queue'};
  const fieldName=path=>(cap.prediction?.fields||[]).find(field=>field.key===path.split('.').at(-1))?.label||fields[path.split('.').at(-1)]||path.split('.').at(-1).replaceAll('_',' ');
  const esc=x=>String(x??'').replace(/[&<>"']/g,ch=>({'&':'&amp;','<':'&lt;','>':'&gt;','"':'&quot;',"'":'&#39;'}[ch]));
  const icon=name=>({trash:'<svg viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="1.7"><path d="M4 7h16M9 7V4h6v3M7 7l1 14h8l1-14M10 10v7m4-7v7"/></svg>',copy:'<svg viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="1.6"><path d="M8 8h12v13H8zM16 8V3H3v13h5"/></svg>',file:'<svg viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="1.7"><path d="M6 2h8l4 4v16H6zM14 2v5h4M12 18v-7m-3 3 3-3 3 3"/></svg>'}[name]);
  const post=(type,fields={})=>parent.postMessage({isStreamlitMessage:true,type,...fields},'*');
  const fit=()=>requestAnimationFrame(()=>{const next=Math.ceil(root.getBoundingClientRect().height+4);if(next!==height){height=next;post('streamlit:setFrameHeight',{height:next});}});
  const notify=(immediate=false)=>{clearTimeout(timer);const send=()=>{const serialized=JSON.stringify(state);if(serialized!==lastSent){lastSent=serialized;post('streamlit:setComponentValue',{value:state});}};if(immediate)send();else timer=setTimeout(send,180);};
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
  function entityTypePicker(card,index){
    const selected=card.kind==='ligand'?(card.mode||'ccd'):card.kind;
    const kinds=cap.localfold?[...cap.types.filter(kind=>kind!=='ligand'),...(cap.ccd?['ccd']:[]),...(cap.smiles?['smiles']:[]),...(cap.modifications?['modification']:[]),...(cap.bonds?['bond']:[]),...(cap.crosslinks?['crosslink']:[]),...(cap.contacts?['contact']:[]),...(cap.pockets?['pocket']:[])]:cap.types;
    return `<div class="type-picker"><button class="type-button" data-action="type-menu" data-index="${index}" aria-haspopup="listbox" aria-expanded="${menu===index}"><span class="dot ${card.kind}"></span>${labels[cap.localfold?selected:card.kind]||card.kind}<span class="chevron">▾</span></button>${menu===index?`<div class="type-menu" role="listbox">${kinds.map(kind=>`<button role="option" aria-selected="${selected===kind}" class="${selected===kind?'selected':''}" data-action="set-type" data-index="${index}" data-kind="${kind}"><span class="dot ${kind}"></span>${labels[kind]}</button>`).join('')}</div>`:''}</div>`;
  }
  function localfoldDecoration(card,index){
    const path=`cards.${index}`,type=entityTypePicker(card,index);
    const polymer=row=>['protein','dna','rna'].includes(row.kind);
    const end=(key,atoms=false)=>{
      const value=card[key]||['','',''];
      return `<div class="localfold-endpoint"><select aria-label="${key} chain" data-path="${path}.${key}.0"><option value="">Chain</option>${chainOptions(value[0],atoms||card.kind==='contact'?()=>true:row=>row.kind==='protein')}</select><input aria-label="${key} residue" type="number" min="1" placeholder="res" data-path="${path}.${key}.1" value="${esc(value[1])}">${atoms?`<input aria-label="${key} atom" placeholder="atom" data-path="${path}.${key}.2" value="${esc(value[2])}">`:''}</div>`;
    };
    let fields;
    if(card.kind==='modification')fields=`<select aria-label="Modified chain" data-path="${path}.chain_id"><option value="">Chain</option>${chainOptions(card.chain_id,polymer)}</select><input type="number" min="1" aria-label="Modified residue" placeholder="res" data-path="${path}.position" value="${esc(card.position)}"><input aria-label="Modification CCD code" placeholder="CCD, e.g. SEP" data-path="${path}.code" value="${esc(card.code)}">`;
    else if(card.kind==='bond')fields=`${end('atom1',true)}<span>–</span>${end('atom2',true)}`;
    else if(card.kind==='contact')fields=`${end('token1',true)}<span>–</span>${end('token2',true)}<label>Max Å<input type="number" min="4" max="20" aria-label="Contact maximum distance" data-path="${path}.max_distance" value="${esc(card.max_distance)}"></label>`;
    else if(card.kind==='pocket')fields=`<select aria-label="Pocket binder" data-path="${path}.binder"><option value="">Binder chain</option>${chainOptions(card.binder)}</select><input aria-label="Pocket residues" placeholder="A, 12; A, 16" data-path="${path}.contacts_input" value="${esc(card.contacts_input)}"><label>Max Å<input type="number" min="4" max="20" aria-label="Pocket maximum distance" data-path="${path}.max_distance" value="${esc(card.max_distance)}"></label>`;
    else fields=`<select aria-label="Crosslinker" data-path="${path}.name">${(cap.localfold_crosslinkers||[]).map(name=>`<option ${card.name===name?'selected':''}>${esc(name)}</option>`).join('')}</select>${end('residue1')}<span>–</span>${end('residue2')}`;
    return `<section class="entity localfold-entity localfold-decoration" data-card="${index}">${type}<div class="localfold-decoration-fields">${fields}</div><button class="localfold-remove" data-action="delete" data-index="${index}" aria-label="Remove ${card.kind}">×</button></section>`;
  }
  function localfoldEntityError(card){
    const chains=catalog();
    if(!decoration(card))return '';
    if(!cap[card.kind==='modification'?'modifications':card.kind==='bond'?'bonds':card.kind==='contact'?'contacts':card.kind==='pocket'?'pockets':'crosslinks'])return `${labels[card.kind]} is unavailable for the selected model.`;
    const valid=(end,protein=false)=>{const chain=chains.find(row=>row.id===end?.[0]);return chain&&(!protein||chain.kind==='protein')&&Number.isInteger(Number(end[1]))&&Number(end[1])>=1&&Number(end[1])<=(['protein','dna','rna'].includes(chain.kind)?chain.value.replace(/\s/g,'').length:1);};
    if(card.kind==='modification')return !valid([card.chain_id,card.position])||!chains.some(row=>row.id===card.chain_id&&['protein','dna','rna'].includes(row.kind))?'Choose a polymer chain and residue for the modification.':!String(card.code||'').trim()?'Enter a modification CCD code.':'';
    if(card.kind==='contact')return !valid(card.token1)||!valid(card.token2)||!Number.isFinite(Number(card.max_distance))||card.max_distance<4||card.max_distance>20?'Complete both contact endpoints and a distance within 4–20 Å.':'';
    if(card.kind==='pocket')return !chains.some(row=>row.id===card.binder)||!card.contacts_input?.trim()?'Choose a binder and pocket residues.':'';
    if(card.kind==='bond')return !valid(card.atom1)||!valid(card.atom2)||!card.atom1[2]?.trim()||!card.atom2[2]?.trim()?'Complete both bond chain, residue and atom references.':'';
    return !valid(card.residue1,true)||!valid(card.residue2,true)?'Choose two protein residues for the crosslink.':'';
  }
  function molecule(card,index){
    if(cap.localfold&&decoration(card))return localfoldDecoration(card,index);
    const chainRows=catalog().filter(row=>row.card===index),chainLabel=chainRows.length>1?`Chains ${chainRows[0]?.id}–${chainRows.at(-1)?.id}`:`Chain ${chainRows[0]?.id||'—'}`;
    const polymer=card.kind!=='ligand',sequence=String(card.value||'').replace(/\s/g,'');
    const type=entityTypePicker(card,index);
    const toolbar=`<div class="toolbar">${type}<div class="copies">Copies <div class="stepper"><button data-action="copies" data-index="${index}" data-delta="-1" aria-label="Fewer copies" ${card.copies<=1?'disabled':''}>−</button><input type="number" aria-label="Copies of chain ${chainRows[0]?.id||''}" min="1" max="52" data-path="cards.${index}.copies" value="${card.copies}"><button data-action="copies" data-index="${index}" data-delta="1" aria-label="More copies" ${card.copies>=52?'disabled':''}>+</button></div></div><span class="spacer"></span>${polymer?`<button class="pill" data-action="plain" data-index="${index}">${card.plain?'View Sequence':'View Plain Text'}</button>`:''}<span class="chain">${chainLabel}</span><button class="icon" data-action="delete" data-index="${index}" title="Remove chain" aria-label="Remove chain ${chainRows[0]?.id||''}">${icon('trash')}</button></div>`;
    let body='';
    if(polymer){
      if(cap.localfold){
        const copies=`<input class="localfold-copies" type="number" aria-label="Copies" min="1" max="52" data-path="cards.${index}.copies" value="${card.copies}">`;
        const optionsOpen=Boolean(card.localfold_options_open);
        const proteinCard=card.kind==='protein';
        const options=proteinCard&&(cap.msa||cap.templates)?`<div class="localfold-options-wrap"><button class="localfold-options-button" data-action="localfold-options" data-index="${index}" aria-label="Options for chain ${chainRows[0]?.id||''}" aria-expanded="${optionsOpen}">⋮</button>${optionsOpen?`<div class="localfold-options"><div class="localfold-options-head"><strong>Options · ${chainLabel}</strong><button data-action="localfold-options-close" data-index="${index}" aria-label="Close options">×</button></div>${cap.msa?`<label>Alignment<select data-path="cards.${index}.localfold_alignment"><option value="job" ${card.localfold_alignment!=='single_sequence'?'selected':''}>The job’s alignment</option><option value="single_sequence" ${card.localfold_alignment==='single_sequence'?'selected':''}>Single sequence</option></select><small>This chain uses the selected job alignment.</small></label>`:''}${cap.templates?`<label>Template<select data-path="cards.${index}.localfold_use_template"><option value="true" ${card.localfold_use_template!==false?'selected':''}>Job templates</option><option value="false" ${card.localfold_use_template===false?'selected':''}>No template</option></select><small>Use uploaded templates for this protein chain.</small></label>`:''}</div>`:''}</div>`:'';
        body=`<div class="localfold-main-row"><div class="localfold-field"><span class="localfold-field-label">Entity type</span>${type}</div><label class="localfold-field localfold-copy-label"><span class="localfold-field-label">Copies</span>${copies}</label><div class="localfold-field localfold-sequence"><span class="localfold-field-label">Sequence</span><textarea aria-label="${labels[card.kind]} sequence for chain ${chainRows[0]?.id||''}" placeholder="Paste or type your ${labels[card.kind].toLowerCase()} sequence here…" data-sequence="${index}" data-path="cards.${index}.value">${esc(card.value)}</textarea><span class="count">${sequence.length?`${sequence.length} ${card.kind==='protein'?'residues':'nucleotides'}`:''}</span></div>${options}<button class="localfold-remove" data-action="delete" data-index="${index}" title="Remove chain" aria-label="Remove chain ${chainRows[0]?.id||''}">×</button></div>`;
      }else{
      let contents='';
      if(!sequence||card.plain||card.editing)contents=`<textarea aria-label="${labels[card.kind]} sequence for chain ${chainRows[0]?.id||''}" placeholder="Paste or type your ${labels[card.kind].toLowerCase()} sequence here…" data-sequence="${index}" data-path="cards.${index}.value">${esc(card.value)}</textarea>`;
      else {let blocks=[];for(let start=0;start<sequence.length;start+=20)blocks.push(`<span class="sequence-block">${esc(sequence.slice(start,start+20))}<sup>${Math.min(start+20,sequence.length)}</sup></span>`);contents=`<div class="styled-sequence" tabindex="0" role="textbox" aria-label="Edit ${labels[card.kind].toLowerCase()} sequence for chain ${chainRows[0]?.id||''}" data-action="edit" data-index="${index}">${blocks.join('')}</div>`;}
      body=`<div class="sequence-area">${contents}<span class="count">${sequence.length?`${sequence.length} ${card.kind==='protein'?'residues':'nucleotides'}`:''}</span></div>`;
      }
      if(card.kind==='protein'&&(cap.cyclic||cap.modifications))body+=`<div class="features">${cap.cyclic?check(`cards.${index}.cyclic`,card.cyclic,'Cyclic'):''}${cap.modifications?`<button class="text-action" data-action="add-modification" data-index="${index}">＋ &nbsp;Add modification</button>`:''}</div>`;
      body+=`<div class="modifications">${(card.modifications||[]).map((mod,j)=>`<div class="modification"><label>Residue position${control(`cards.${index}.modifications.${j}.position`,mod.position,'number','min="1"')}</label><label>CCD code${control(`cards.${index}.modifications.${j}.code`,mod.code,'text','placeholder="MSE"')}</label><button class="icon" data-action="remove-modification" data-index="${index}" data-mod="${j}" aria-label="Remove modification">${icon('trash')}</button></div>`).join('')}</div>`;
    }else{
      const mode=card.mode||(cap.ccd?'ccd':'smiles');
      body=`<div class="ligand-help">Ligand SMILES or CCD code (use CCD code for ions)</div><div class="tabs" role="tablist">${(cap.ccd?['ccd']:[]).concat(cap.smiles?['smiles']:[]).map(m=>`<button role="tab" aria-selected="${mode===m}" class="${mode===m?'active':''}" data-action="ligand-mode" data-index="${index}" data-mode="${m}">${m==='ccd'?'Enter CCD Code':'SMILES'}</button>`).join('')}</div><input class="ligand-input" aria-label="${mode==='ccd'?'CCD code':'Ligand SMILES'}" placeholder="${mode==='ccd'?'Enter CCD component code, e.g. ATP':'Paste or type ligand SMILES here…'}" data-path="cards.${index}.value" value="${esc(card.value)}"><div data-preview-slot="${card.id}">${mode==='smiles'&&previews[String(card.id)]?`<img class="ligand-preview" alt="Ligand chemical structure" src="${esc(previews[String(card.id)])}">`:''}</div>`;
    }
    const compactLocalfold=cap.localfold&&polymer;
    return `<section class="entity ${compactLocalfold?'localfold-entity':''}" data-card="${index}">${compactLocalfold?'':handle(index)}${compactLocalfold?'':toolbar}${body}</section>`;
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
      else if(field.type==='text')input=field.multiline?`<textarea aria-label="${esc(field.label)}" data-path="${path}" rows="5">${esc(value)}</textarea>`:control(path,value,'text');
      else input=control(path,value,'number',`min="${field.min}" max="${field.max}" step="${field.step||(field.type==='float'?0.1:1)}"`);
      return setting(esc(field.label),esc(field.help||''),input,field.type==='bool');
    };
    const fields=config.fields.filter(field=>!field.required&&!field.hidden);
    const templates=cap.templates?templateControls():'';
    if(cap.localfold){
      const model=String(s.model||'alphafold3'),meta=cap.localfold_models?.[model]||{};
      const isEsm=model.startsWith('esmfold2'),isAf2=model.startsWith('af2_'),isBioemu=model==='bioemu2';
      const family=isEsm?'esmfold2':model;
      const globalMsa=s.msa_mode!=='single_sequence'&&!meta.no_msa;
      const primaryKeys=new Set(['num_recycles','seed','num_sampling_steps','af2_num_models','bioemu2_mode','msa_mode','max_msa_seqs','max_extra_msa_seqs']);
      const primaryOrder=['af2_num_models','seed','num_recycles','num_sampling_steps','bioemu2_mode','msa_mode','max_msa_seqs','max_extra_msa_seqs'];
      const primary=fields.filter(field=>primaryKeys.has(field.key)
        &&(field.key!=='num_sampling_steps'||!isAf2)
        &&(field.key!=='af2_num_models'||isAf2)
        &&(field.key!=='bioemu2_mode'||isBioemu)
        &&(field.key!=='msa_mode'||!meta.no_msa)
        &&(!['max_msa_seqs','max_extra_msa_seqs'].includes(field.key)||globalMsa))
        .sort((left,right)=>primaryOrder.indexOf(left.key)-primaryOrder.indexOf(right.key));
      const advanced=fields.filter(field=>!primaryKeys.has(field.key)&&field.key!=='model'&&(!['bioemu2_lag_ns','bioemu2_propagation_steps'].includes(field.key)||Boolean(s.bioemu2_structure_text)));
      const compactField=field=>{
        const path=`settings.${field.key}`,value=s[field.key]??field.default;
        let input;
        if(field.key==='af2_num_models')input=select(path,value,[[1,'1'],[2,'2'],[3,'3'],[4,'4'],[5,'All 5']]);
        else if(field.type==='select')input=select(path,value,field.choices.map(choice=>Array.isArray(choice)?choice:[choice,choice]));
        else input=control(path,value,'number',`min="${field.min}" max="${field.max}" step="${field.step||(field.type==='float'?0.1:1)}"`);
        let label=field.label.replace('Number of ','');
        if(field.key==='num_sampling_steps')label=isBioemu?'Steps':'Diffusion';
        if(field.key==='af2_num_models')label='AF2 models';
        return `<label class="localfold-control" data-setting="${esc(field.key)}" title="${esc(field.help||field.label)}"><span>${esc(label)}</span>${input}</label>`;
      };
      const modelGroups=cap.localfold_model_groups||[];
      const modelOptions=modelGroups.map(([group,rows])=>`<optgroup label="${esc(group)}">${rows.map(([id,label])=>`<option value="${esc(id)}" ${family===id?'selected':''}>${esc(label)}</option>`).join('')}</optgroup>`).join('');
      const modelControl=`<label class="localfold-control localfold-model-control" data-setting="model"><span>Model</span><select aria-label="LocalFold model" data-path="settings.model">${modelOptions}</select></label>`;
      const variants=cap.localfold_esmfold2_variants||[];
      const variantControl=isEsm?`<label class="localfold-control" data-setting="model-variant"><span>Variant</span><select aria-label="ESMFold2 variant" data-path="settings.model" data-localfold-model-variant>${variants.map(([id,label])=>`<option value="${esc(id)}" ${model===id?'selected':''}>${esc(label)}</option>`).join('')}</select></label>`:'';
      const msaValues=isAf2?['512:1024','256:512','128:256','64:128','32:64','16:32']:
        [...new Set([String(meta.msa_depth||s.max_msa_seqs||1024),'512','256','128','64','32','16'])];
      const currentMsa=isAf2?`${s.max_msa_seqs||128}:${s.max_extra_msa_seqs||256}`:String(s.max_msa_seqs??meta.msa_depth??1024);
      const ownMsaCap=isBioemu||(isEsm&&!meta.no_msa);
      const msaChoices=isAf2?msaValues.map(value=>[value,value]):
        msaValues.map(value=>[value,ownMsaCap&&value===String(meta.msa_depth)?`Default (${value})`:value]);
      if(!isAf2&&isBioemu)msaChoices[0]=[String(meta.msa_depth||4092),`Default (${meta.msa_depth||4092})`];
      const msaControl=globalMsa?`<label class="localfold-control" data-setting="max_msa"><span>Max MSA</span><select aria-label="Maximum MSA sequences" data-localfold-msa>${msaChoices.map(([value,label])=>`<option value="${esc(value)}" ${currentMsa===value?'selected':''}>${esc(label)}</option>`).join('')}</select></label>`:'';
      const primaryControls=primary.filter(field=>!['msa_mode','max_msa_seqs','max_extra_msa_seqs'].includes(field.key)).map(compactField).join('');
      const alignmentField=primary.find(field=>field.key==='msa_mode');
      const alignmentControl=alignmentField?compactField(alignmentField):'';
      const bioemuSamples=isBioemu?`<label class="localfold-control" data-setting="samples"><span>Samples</span><select aria-label="BioEmu2 samples" data-path="settings.samples">${[1,2,3,5,8,10,25,50,100].map(value=>`<option value="${value}" ${Number(s.samples||1)===value?'selected':''}>${value}</option>`).join('')}</select></label>`:'';
      const msaFiles=cap.msa&&globalMsa&&!cap.settings_only?setting('Custom A3M files','Upload an alignment and assign it to a protein chain. This skips MSA search for that chain.',upload('msa_files','your a3m files','.a3m,.txt')+files('msa_files',state.msa_files)):'';
      const propagationUpload=isBioemu?setting('Initial structure for propagation (optional)','Upload an mmCIF of this input to propagate it in equilibrium mode.',upload('bioemu_structure','an initial mmCIF','.cif,.mmcif',false)+(s.bioemu2_structure_text?`<div class="file-row"><span>${esc(state.bioemu_structure_name||'Initial mmCIF')}</span><button data-action="remove-bioemu-structure">×</button></div>`:'')):'';
      const advancedOptions=advanced.map(field=>(field.key==='use_msa'?templates:'')+renderField(field)).join('')+(!advanced.some(field=>field.key==='use_msa')?templates:'');
      const outputSamples=!isAf2&&!isBioemu?setting('Output samples','Number of structures to generate for this model and seed.',control('settings.samples',s.samples,'number',`min="${min}" max="${max}"`)):'';
      const requiredOptions=config.fields.filter(field=>field.required&&!field.hidden).map(renderField).join('');
      return `<div class="localfold-primary-settings">${modelControl}${primaryControls}${variantControl}${bioemuSamples}${alignmentControl}${msaControl}</div>${requiredOptions}<details class="optional" ${state.optional_open?'open':''}><summary>Advanced settings</summary>${outputSamples}${msaFiles}${propagationUpload}${advancedOptions}${(config.fixed_settings||[]).map(text=>`<p class="small-note">${esc(text)}</p>`).join('')}</details>`;
    }
    const options=fields.map(field=>(field.key==='use_msa'?templates:'')+renderField(field)).join('')+(!fields.some(field=>field.key==='use_msa')?templates:'');
    return `<label class="samples">Number of Samples *${control('settings.samples',s.samples,'number',`min="${min}" max="${max}"`)}</label><details class="optional" ${state.optional_open?'open':''}><summary>Optional settings</summary>${options}${(config.fixed_settings||[]).map(text=>`<p class="small-note">${esc(text)}</p>`).join('')}</details>${config.fields.filter(field=>field.required&&!field.hidden).map(renderField).join('')}`;
  }
  function templateControls(){return setting('Templates','Upload PDB or mmCIF structures. Protein chains are aligned by sequence; select a chain to restrict the template.',upload('templates','your cif/pdb files','.cif,.mmcif,.pdb,.ent')+files('templates',state.templates));}
  function error(){
    if(cap.prediction){
      const s=state.settings,min=cap.prediction.sample_min_modes?.[s.optimization_mode]||1;
      const max=cap.prediction.sample_max||5;
      if(!Number.isInteger(Number(s.samples))||s.samples<min||s.samples>max)return `Request between ${min} and ${max} samples for this inference mode.`;
      for(const field of cap.prediction.fields){const value=s[field.key]??field.default;if(field.required&&!value)return field.key==='accept_helixfold3_terms'?'Review and accept HelixFold3’s published non-commercial weight terms.':'Review and accept AlphaFold3’s published weight terms.';if(!['bool','select','text'].includes(field.type)&&(!Number.isFinite(Number(value))||value===''||Number(value)<field.min||Number(value)>field.max))return `Enter a valid ${field.label.toLowerCase()}.`;}
    }
    if(cap.settings_only)return '';
    if(state.settings?.use_yaml){if(!state.yaml_text?.trim())return 'Upload a Boltz YAML file.';if(yamlError)return yamlError;if(hasLigand()&&!state.compound_name?.trim())return 'Enter a compound name for this ligand complex.';return '';}
    if(!state.cards.some(card=>!decoration(card))||state.cards.some(card=>!decoration(card)&&!String(card.value||'').trim()))return 'All sequences must be non empty.';
    if(state.cards.filter(card=>!decoration(card)).reduce((sum,card)=>sum+Number(card.copies||1),0)>52)return 'Enter at most 52 chains and ligand copies.';
    if(cap.localfold){const problem=state.cards.map(localfoldEntityError).find(Boolean);if(problem)return problem;}
    if(state.cards.some(card=>(card.modifications||[]).some(mod=>!mod.code?.trim())))return 'Every modification needs a CCD code.';
    if((state.msa_files||[]).some(file=>!file.chain_id)&&(cap.localfold?state.settings.msa_mode!=='single_sequence':state.settings.use_msa))return 'Assign each A3M file to its protein chain.';
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
      if(cap.localfold)html+=`<div class="localfold-entity-header"><strong>Entities</strong><button class="outline" data-action="add" ${state.cards.length>=52?'disabled':''}>＋ Add entity</button></div>`;
      html+=state.cards.map(molecule).join('')+`<div class="actions">${cap.localfold?'':`<button class="outline" data-action="add" ${state.cards.length>=16?'disabled':''}>＋ &nbsp;Add chain</button>`}<button class="outline" data-action="fasta">${icon('file')} Import FASTA</button><button class="clear" data-action="clear">Clear all</button><input class="hidden" type="file" data-file="fasta" accept=".fa,.fasta,.fna,.txt"></div>`;
      if(hasLigand())html+=compoundIdentity();
      if(cap.restraints)html+=(state.restraints.length?'<div class="caption">RESTRAINTS</div>':'')+state.restraints.map(restraint).join('')+`<div class="restraint-help"><button class="outline" data-action="add-restraint">＋ &nbsp;Add restraint</button><span>Restraints constrain how different chains are arranged in your complex, and can incorporate prior knowledge of structural contacts.</span></div>`;
      const ligands=catalog().filter(row=>['ccd','smiles'].includes(row.kind));
      if(cap.affinity_control&&ligands.length)html+=`<div class="affinity">${check('settings.predict_affinity',state.settings.predict_affinity,'Predict affinity (small molecule ligand)')}${state.settings.predict_affinity?`<label>Binder chain <select data-path="settings.affinity_binder">${chainOptions(state.settings.affinity_binder,row=>['ccd','smiles'].includes(row.kind))}</select></label>`:''}</div>`;
    }
    if(!cap.boltz&&!cap.prediction&&(cap.templates||cap.msa))html+=`<details class="optional" ${state.optional_open?'open':''}><summary>Optional settings</summary>${cap.templates?templateControls():''}${cap.msa?setting('Custom A3M Files','Assign an alignment to its protein chain.',upload('msa_files','your a3m files','.a3m,.txt')+files('msa_files',state.msa_files)):''}</details>`;
    if(cap.boltz||cap.prediction)html+=(cap.boltz?boltzOptions():predictionOptions())+`<div class="submit-area"><button class="submit" data-action="submit" ${error()||!cap.can_submit?'disabled':''}>Submit</button><div class="validation">${esc(error()||(!cap.can_submit?'Select a project to submit.':''))}</div></div>`;
    root.innerHTML=html;
    root.querySelector('.optional')?.addEventListener('toggle',event=>{state.optional_open=event.target.open;notify();fit();});
    root.querySelectorAll('[data-sequence]').forEach(input=>input.addEventListener('blur',event=>{const card=state.cards[Number(input.dataset.sequence)];card.value=card.value.replace(/\s/g,'').toUpperCase();card.editing=false;notify();if(!cap.localfold&&!card.plain&&!root.contains(event.relatedTarget))render();}));
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
    const path=node.dataset.path||'';
    if(!path&&node.dataset.localfoldMsa===undefined)return;
    let value=node.type==='checkbox'?node.checked:node.value;
    if(node.type==='number')value=node.value===''?'':Number(node.value);
    if(cap.localfold&&['settings.af2_num_models','settings.num_recycles','settings.num_sampling_steps','settings.samples'].includes(path))value=Number(value);
    if(cap.localfold&&/cards\.\d+\.(?:code|(?:atom[12]|token[12])\.2)$/.test(path))value=String(value).trim().toUpperCase();
    if(path.endsWith('.copies'))value=Math.max(1,Math.min(52,Number(value)||1));
    const modelPath=path==='settings.model'&&cap.localfold;
    if(modelPath){
      const previous=String(state.settings.model||'alphafold3');
      const requested=String(value);
      value=requested==='esmfold2'&&node.dataset.localfoldModelVariant===undefined
        ?(previous.startsWith('esmfold2')?previous:(cap.localfold_default_esmfold2_variant||'esmfold2_lm600m'))
        :requested;
      const previousMeta=cap.localfold_models?.[previous]||{};
      const selected=cap.localfold_models?.[value];
      setPath(path,value);
      if(selected){
        if(!previousMeta.no_msa&&selected.no_msa){state.localfold_saved_msa_mode=state.settings.msa_mode||'automatic';state.settings.msa_mode='single_sequence';}
        else if(previousMeta.no_msa&&!selected.no_msa){state.settings.msa_mode=state.localfold_saved_msa_mode||'automatic';delete state.localfold_saved_msa_mode;}
        state.settings.samples=selected.samples;
        state.settings.num_recycles=selected.recycles;
        state.settings.num_sampling_steps=selected.diffusion_steps??0;
        state.settings.max_msa_seqs=selected.msa_depth;
        state.settings.max_extra_msa_seqs=selected.extra_msa_depth??256;
      }
      if(previous==='bioemu2'&&value!=='bioemu2')for(const field of cap.localfold_profiles?.bioemu2?.prediction.fields||[])if(field.key.startsWith('bioemu2_'))state.settings[field.key]=field.default;
      applyLocalfoldProfile();normalizePrediction();
    }else if(node.dataset.localfoldMsa!==undefined){
      if(String(value).includes(':')){
        const [clusters,extra]=String(value).split(':').map(Number);
        state.settings.max_msa_seqs=clusters;state.settings.max_extra_msa_seqs=extra;
      }else state.settings.max_msa_seqs=Number(value);
    }else if(path)setPath(path,value);
    const presets=cap.prediction?.parameter_presets;
    if(presets&&path){
      const presetKey=presets.selector_key, presetPath=`settings.${presetKey}`, presetValues=presets.values||{};
      const customLabel=presets.custom_label||'Custom';
      const presetFields=new Set(Object.keys(Object.values(presetValues)[0]||{}));
      if(path===presetPath){
        const values=presetValues[value];
        if(values)for(const [key,presetValue] of Object.entries(values))state.settings[key]=presetValue;
      }else if(presetFields.has(path.split('.').at(-1))){
        const match=Object.entries(presetValues).find(([,values])=>Object.entries(values).every(([key,presetValue])=>Number(state.settings[key])===Number(presetValue)));
        state.settings[presetKey]=match?.[0]||customLabel;
      }
    }
    if(path==='settings.predict_affinity'&&value&&!state.settings.affinity_binder)state.settings.affinity_binder=catalog().find(row=>['ccd','smiles'].includes(row.kind))?.id||'';
    const match=path.match(/^restraints\.(\d+)\.pocket\./);if(match)updatePocket(Number(match[1]));
    root.querySelectorAll('[data-restraint-card]').forEach(section=>{const index=Number(section.dataset.restraintCard),message=restraintError(state.restraints[index]);section.style.borderColor=message?'#ef5555':'';section.querySelector('[data-restraint-error]').textContent=message;});
    const validation=root.querySelector('.validation'),submit=root.querySelector('.submit');
    if(validation)validation.textContent=error()||(!cap.can_submit?'Select a project to submit.':'');if(submit)submit.disabled=Boolean(error()||!cap.can_submit);
    notify();
  }
  root.addEventListener('input',event=>{if(event.target.tagName!=='SELECT')updateInput(event.target);});
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
    if(action==='localfold-options'){card.localfold_options_open=!card.localfold_options_open;render();return;}
    if(action==='localfold-options-close'){card.localfold_options_open=false;render();return;}
    if(action==='set-type'){
      const kind=node.dataset.kind,chainRows=catalog(),a=chainRows[0]?.id||'',b=chainRows[1]?.id||a;
      Object.assign(card,{kind:['ccd','smiles'].includes(kind)?'ligand':kind,value:'',modifications:[],cyclic:false,copies:cap.localfold?1:card.copies,mode:['ccd','smiles'].includes(kind)?kind:(cap.ccd?'ccd':'smiles')});
      if(kind==='modification')Object.assign(card,{chain_id:chainRows.find(row=>['protein','dna','rna'].includes(row.kind))?.id||'',position:'',code:''});
      if(kind==='bond')Object.assign(card,{atom1:[a,'',''],atom2:[b,'','']});
      if(kind==='contact')Object.assign(card,{token1:[a,'',''],token2:[b,'',''],max_distance:6});
      if(kind==='pocket')Object.assign(card,{binder:b,contacts_input:'',max_distance:6});
      if(kind==='crosslink')Object.assign(card,{name:'DSSO',residue1:[a,''],residue2:[b,'']});
      menu=null;
    }
    if(action==='add'){state.cards.push(fresh());if(cap.localfold)menu=state.cards.length-1;}
    if(action==='delete'){state.cards.splice(index,1);if(!state.cards.length)state.cards.push(fresh());}
    if(action==='duplicate'&&state.cards.length<16)state.cards.splice(index+1,0,{...clone(card),id:Date.now()});
    if(action==='copies')card.copies=Math.max(1,Math.min(52,card.copies+Number(node.dataset.delta)));
    if(action==='plain')card.plain=!card.plain;
    if(action==='edit'){card.editing=true;render();root.querySelector(`[data-sequence="${index}"]`)?.focus();return;}
    if(action==='add-modification'){if(cap.localfold){const chain=catalog().find(row=>row.card===index);state.cards.push({...fresh(),kind:'modification',chain_id:chain?.id||'',position:'',code:''});}else(card.modifications||=([])).push({position:1,code:''});}
    if(action==='remove-modification')card.modifications.splice(Number(node.dataset.mod),1);
    if(action==='ligand-mode'){card.mode_values||={};card.mode_values[card.mode||(cap.ccd?'ccd':'smiles')]=card.value;card.mode=node.dataset.mode;card.value=card.mode_values[card.mode]||'';}
    if(action==='fasta'){root.querySelector('[data-file="fasta"]').click();return;}
    if(action==='clear'){state.cards=[fresh()];state.restraints=[];state.msa_files=[];state.settings.predict_affinity=false;}
    if(action==='add-restraint')state.restraints.push(newRestraint());
    if(action==='delete-restraint')state.restraints.splice(index,1);
    if(action==='duplicate-restraint')state.restraints.splice(index+1,0,clone(state.restraints[index]));
    if(action==='select-residues'){const d=state.restraints[index].pocket;const sequence=catalog().find(row=>row.id===d.pocket_chain)?.value||'';d.pocket_residues=node.dataset.all==='1'?[...sequence.replace(/\s/g,'')].map((_,i)=>i+1).join(' '):'';updatePocket(index);}
    if(action==='randomize'){const values=new Uint32Array(1);crypto.getRandomValues(values);state.settings.seed=values[0]%2147483646+1;}
    if(action==='remove-bioemu-structure'){state.settings.bioemu2_structure_text='';state.bioemu_structure_name='';}
    if(action==='remove-file'){if(node.dataset.kind==='yaml_files'&&node.closest('.setting')?.querySelector('[data-file="yaml"]')&&node.closest('.file-row').querySelector('.filename').textContent===state.yaml_name){state.yaml_text='';state.yaml_name='';}else state[node.dataset.kind].splice(index,1);}
    if(action==='example'){try{if(!examples[node.dataset.name])examples=await(await fetch('examples.json')).json();const example=examples[node.dataset.name];if(example){state.cards=clone(example.cards).filter(card=>{const kind=card.kind==='ligand'?(card.mode||'smiles'):card.kind;return kind!=='ligand'||cap.ligand_examples;});state.cards.filter(card=>card.kind==='ligand').forEach(card=>card.preview_value=card.value);state.restraints=cap.restraints?clone(example.restraints||[]):[];state.msa_files=[];state.settings.use_yaml=false;state.settings.predict_affinity=Boolean(cap.affinity_control&&example.affinity);state.settings.affinity_binder=state.settings.predict_affinity?(catalog().find(row=>['ccd','smiles'].includes(row.kind))?.id||''):'';}}catch{state.message='Could not load the example.';}}
    if(action==='submit'){if(error()||!cap.can_submit)return;state.submit_token=String(Date.now())+'-'+Math.random().toString(36).slice(2);notify(true);return;}
    render();notify(true);
  });
  let selection=null;
  root.addEventListener('pointerdown',event=>{const node=event.target.closest('[data-residue]');if(!node)return;event.preventDefault();selection={index:Number(node.dataset.restraint),start:Number(node.dataset.residue),last:Number(node.dataset.residue)};});
  root.addEventListener('pointerover',event=>{const node=event.target.closest('[data-residue]');if(selection&&node&&Number(node.dataset.restraint)===selection.index)selection.last=Number(node.dataset.residue);});
  window.addEventListener('pointerup',()=>{if(!selection)return;const d=state.restraints[selection.index].pocket,selected=new Set((d.pocket_residues||'').split(/[\s,]+/).filter(Boolean).map(Number));for(let i=Math.min(selection.start,selection.last);i<=Math.max(selection.start,selection.last);i++){if(selection.start===selection.last&&selected.has(i))selected.delete(i);else selected.add(i);}d.pocket_residues=[...selected].sort((a,b)=>a-b).join(' ');updatePocket(selection.index);selection=null;render();notify(true);});
  root.addEventListener('keydown',event=>{if((event.key==='Enter'||event.key===' ')&&event.target.dataset.upload){event.preventDefault();event.target.querySelector('input[type=file]').click();}if(event.target.dataset.action==='edit'&&(event.key==='Enter'||event.key===' ')){event.preventDefault();event.target.click();}if(event.target.dataset.residue&&(event.key==='Enter'||event.key===' ')){event.preventDefault();const index=Number(event.target.dataset.restraint),position=Number(event.target.dataset.residue),d=state.restraints[index].pocket,positions=new Set(d.pocket_residues.split(/[\s,]+/).filter(Boolean).map(Number));if(positions.has(position))positions.delete(position);else positions.add(position);d.pocket_residues=[...positions].sort((a,b)=>a-b).join(' ');updatePocket(index);render();notify(true);}if(event.key==='Escape'){let changed=menu!==null;menu=null;state.cards.forEach(card=>{if(card.localfold_options_open){card.localfold_options_open=false;changed=true;}});if(changed)render();}});
  async function readFiles(kind,fileList){
    for(const file of fileList){
      if(file.size>20*1024*1024){state.message='Each input file must be smaller than 20 MB.';continue;}
      const text=await file.text();
      if(kind==='fasta'){const records=text.split(/^>/m).map(record=>{const lines=record.trim().split(/\r?\n/);return record.startsWith('>')?'':lines.slice(1).join('').replace(/\s/g,'');}).filter(Boolean);if(!text.trim().startsWith('>')){state.message='Import a FASTA file beginning with a > header.';continue;}if(!records.length||records.length>16){state.message='Import 1 to 16 FASTA records.';continue;}state.cards=records.map(value=>({...fresh(),value:value.toUpperCase()}));state.restraints=[];state.msa_files=[];}
      else if(kind==='bioemu_structure'){state.settings.bioemu2_structure_text=text;state.settings.bioemu2_mode='equilibrium';state.bioemu_structure_name=file.name;}
      else if(kind==='yaml'){state.yaml_text=text;state.yaml_name=file.name;}
      else{const item={name:file.name,text};if(kind==='msa_files'){const query=text.split(/^>/m)[1]?.split(/\r?\n/).slice(1).join('').replace(/[a-z.\-\s]/g,'')||'';item.chain_id=catalog().find(row=>row.kind==='protein'&&row.value.replace(/\s/g,'').toUpperCase()===query)?.id||'';}state[kind]=state[kind].filter(existing=>existing.name!==file.name);state[kind].push(item);}
    }
    render();if(state.message){const p=document.createElement('p');p.className='error';p.textContent=state.message;root.append(p);state.message='';fit();}notify(true);
  }
  window.addEventListener('message',event=>{if(event.data.type!=='streamlit:render')return;const signature=JSON.stringify(event.data.args.capabilities),capabilitiesChanged=serverCapSignature!==signature;serverCapSignature=signature;cap=clone(event.data.args.capabilities);if(state)applyLocalfoldProfile();previews=event.data.args.previews||{};const changed=JSON.stringify(yamlCatalog)!==JSON.stringify(event.data.args.yaml_chains||[]);yamlCatalog=event.data.args.yaml_chains||[];yamlError=event.data.args.yaml_error||'';if(!state){state=clone(event.data.args.initial);state.settings||={};state.templates||=[];state.msa_files||=[];state.yaml_files||=[];state.compound_name||='';state.compound_id||='';applyLocalfoldProfile();normalizePrediction();render();}else {state.compound_name||='';state.compound_id||='';if(changed||capabilitiesChanged){applyLocalfoldProfile();normalizePrediction();render();return;}const submit=root.querySelector('.submit'),validation=root.querySelector('.validation');if(submit)submit.disabled=Boolean(error()||!cap.can_submit);if(validation)validation.textContent=error()||(!cap.can_submit?'Select a project to submit.':'');root.querySelectorAll('[data-preview-slot]').forEach(slot=>{const source=previews[slot.dataset.previewSlot];slot.innerHTML=source?`<img class="ligand-preview" alt="Ligand chemical structure" src="${esc(source)}">`:'';});fit();}});
  new ResizeObserver(fit).observe(root);
  post('streamlit:componentReady',{apiVersion:1});
})();
