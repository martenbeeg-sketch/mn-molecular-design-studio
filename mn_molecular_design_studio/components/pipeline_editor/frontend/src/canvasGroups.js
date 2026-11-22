// Groups change the drawing. Wires always resolve to the original tool ports.
const socket = (node, port) => `${node}/${port}`;

export function pruneGroups(groups, nodes) {
  const ids = new Set(nodes.map(node => node.id));
  return groups.map(group => ({...group, members: group.members.filter(id => ids.has(id))}))
    .filter(group => group.members.length > 1);
}

export function frameBounds(nodes) {
  const left = Math.min(...nodes.map(node => node.position.x));
  const top = Math.min(...nodes.map(node => node.position.y));
  const right = Math.max(...nodes.map(node => node.position.x + (node.measured?.width || node.style?.width || 286)));
  const bottom = Math.max(...nodes.map(node => node.position.y + (node.measured?.height || node.style?.height || 185)));
  return {position: {x: left - 24, y: top - 50}, width: right - left + 48, height: bottom - top + 75};
}

export function projectGroups(nodes, edges, authoredEdges, groups, catalog) {
  const byId = new Map(nodes.map(node => [node.id, node]));
  const realEdges = new Map(authoredEdges.map(edge => [edge.id, edge]));
  const hiddenBy = new Map(), endpoints = new Map(), frames = [];
  for (const group of pruneGroups(groups, nodes)) {
    const members = new Set(group.members);
    const children = group.members.map(id => byId.get(id));
    const {position, width, height} = frameBounds(children);
    const boundary = {in: [], out: []};
    if (group.collapsed) {
      for (const node of children) {
        hiddenBy.set(node.id, group.id);
        const spec = catalog[node.data.kind];
        for (const [side, ports] of [['in', spec.inputs], ['out', spec.outputs]]) {
          for (const port of ports) {
            const onInput = side === 'in';
            const linked = authoredEdges.filter(edge => onInput
              ? edge.target === node.id && (edge.targetHandle || 'in') === port.id
              : edge.source === node.id && (edge.sourceHandle || 'out') === port.id);
            const external = linked.some(edge => !members.has(onInput ? edge.source : edge.target));
            if (external || !linked.length) {
              const handle = socket(node.id, port.id);
              boundary[side].push({handle, label: `${node.data.label} · ${port.label}`});
              endpoints.set(`${group.id}:${side}:${handle}`, {id: node.id, port: port.id});
            }
          }
        }
      }
    }
    frames.push({id: group.id, type: group.collapsed ? 'collapsedGroup' : 'frame',
      position, draggable: group.collapsed, deletable: false, zIndex: group.collapsed ? 0 : -1,
      style: group.collapsed ? {width: 282} : {width, height},
      data: {...group, boundary}});
  }
  const projected = edges.flatMap(edge => {
    const source = hiddenBy.get(edge.source), target = hiddenBy.get(edge.target);
    if (source && source === target) return [];
    const real = realEdges.get(edge.id) || edge;
    return [{...edge,
      ...(source ? {source, sourceHandle: socket(real.source, real.sourceHandle || 'out')} : {}),
      ...(target ? {target, targetHandle: socket(real.target, real.targetHandle || 'in')} : {}),
    }];
  });
  return {nodes: [...frames, ...nodes.filter(node => !hiddenBy.has(node.id))], edges: projected, endpoints};
}

export function originalConnection(connection, endpoints) {
  const source = endpoints.get(`${connection.source}:out:${connection.sourceHandle}`);
  const target = endpoints.get(`${connection.target}:in:${connection.targetHandle}`);
  return {...connection, source: source?.id || connection.source, sourceHandle: source?.port || connection.sourceHandle,
    target: target?.id || connection.target, targetHandle: target?.port || connection.targetHandle};
}
