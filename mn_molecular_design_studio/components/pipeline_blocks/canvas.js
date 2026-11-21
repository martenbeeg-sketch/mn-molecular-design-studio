let flow = null, pendingOutput = null, flowOrigin = null;

function flowPoint(dom, side) {
  const rect = dom.getBoundingClientRect();
  const frame = flow.viewport.getBoundingClientRect();
  return [rect.left + rect.width / 2 - frame.left - 1, rect[side] - frame.top - 1];
}

function coordinates(left, top) {
  return {
    x: Math.round((Math.max(24, top) - flow.origin.top) / .55),
    y: Math.round((Math.max(24, left) - flow.origin.left) / 1.45),
  };
}

function displayPosition({x, y}, managed) {
  return managed ? {left: x * .76, top: y * .59} : {left: y * 1.45, top: x * .55};
}

function drawWires(pointer = null) {
  const {viewport, canvas, ports, edges} = flow;
  const width = viewport.clientWidth, height = viewport.clientHeight;
  if (!width || !height) return;
  const ratio = window.devicePixelRatio || 1;
  if (canvas.width !== Math.round(width * ratio) || canvas.height !== Math.round(height * ratio)) {
    canvas.width = Math.round(width * ratio);
    canvas.height = Math.round(height * ratio);
    canvas.style.width = width + "px";
    canvas.style.height = height + "px";
  }
  const ctx = canvas.getContext("2d");
  ctx.setTransform(ratio, 0, 0, ratio, 0, 0);
  ctx.clearRect(0, 0, width, height);
  const colors = getComputedStyle(viewport);
  function wire(source, end, emphasis = false, loop = false, label = "") {
    const start = flowPoint(source, "bottom");
    const lift = Math.max(25, Math.min(100, Math.abs(end[1] - start[1]) / 2));
    ctx.strokeStyle = colors.getPropertyValue(emphasis ? "--active-wire" : "--wire").trim();
    ctx.lineWidth = emphasis ? 2.5 : 1.6;
    ctx.setLineDash(loop ? [6, 4] : []);
    ctx.beginPath();
    ctx.moveTo(...start);
    ctx.bezierCurveTo(start[0], start[1] + lift, end[0], end[1] - lift, ...end);
    ctx.stroke();
    ctx.setLineDash([]);
    if (label) {
      ctx.font = "600 10px system-ui";
      const x = (start[0] + end[0]) / 2, y = (start[1] + end[1]) / 2;
      ctx.fillStyle = colors.getPropertyValue("--block").trim() || "#fbfdfb";
      ctx.fillRect(x - 4, y - 11, ctx.measureText(label).width + 8, 15);
      ctx.fillStyle = colors.getPropertyValue("--text").trim() || "#253b34";
      ctx.fillText(label, x, y);
    }
  }
  for (const edge of edges) {
    const source = ports.get(portKey(edge.source, "out", edge.sourceHandle || "out"));
    const target = ports.get(portKey(edge.target, "in", edge.targetHandle || "in"));
    if (!source || !target) continue;
    const emphasis = [edge.source, edge.target].includes(selected);
    const decision = flow.model.nodes[edge.source].kind === "active_decision";
    let label = "";
    if (edge.loop_back) label = "Next round";
    else if (decision) label = flow.model.nodes[edge.target].kind === "active_results" ? "Stop" : "Continue";
    wire(source, flowPoint(target, "top"), emphasis, edge.loop_back, label);
  }
  if (pointer && pendingOutput) {
    const {left, top} = viewport.getBoundingClientRect();
    wire(pendingOutput, [pointer.clientX - left, pointer.clientY - top], true);
  }
}

function portKey(id, direction, port) {
  return [id, direction, port].join("/");
}

function canvasMatches(output, kind) {
  const {id, port} = output.dataset;
  const source = flow.model.nodes[id];
  return flow.model.compatible[source.kind]?.[port]?.filter(([candidate]) => candidate === kind) || [];
}

function canWire(output, input) {
  const target = input.dataset.id;
  if (output.dataset.id === target || input.dataset.direction !== "in") return false;
  const kind = flow.model.nodes[target].kind;
  if (!canvasMatches(output, kind).some(([, port]) => port === input.dataset.port)) return false;
  const occupancy = flow.edges.filter(edge => edge.target === target && (edge.targetHandle || "in") === input.dataset.port);
  if (occupancy.length && (input.dataset.multiple !== "true" || occupancy.some(edge =>
    edge.source === output.dataset.id && (edge.sourceHandle || "out") === output.dataset.port))) return false;
  return !flow.ancestors[output.dataset.id]?.has(target);
}

function refreshFlow() {
  if (!flow) return;
  const block = dragging || picked;
  flow.placeButton.hidden = !picked;
  flow.placeButton.textContent = picked ? "Place " + args.cards[picked].label : "";
  for (const node of flow.nodes.values()) {
    const active = node.dataset.id === selected;
    node.classList.toggle("selected", active);
    node.querySelector(".flow-step").setAttribute("aria-pressed", String(active));
  }
  for (const port of flow.ports.values()) {
    const outgoing = port.dataset.direction === "out";
    const allowed = outgoing ? block && canvasMatches(port, block).length === 1
      : pendingOutput && canWire(pendingOutput, port);
    port.classList.toggle("available", Boolean(allowed));
    port.classList.toggle("mismatch", outgoing ? Boolean(block && !allowed) : Boolean(pendingOutput && !allowed));
    if (outgoing && !flow.model.managed) port.setAttribute("aria-pressed", String(port === pendingOutput));
    port.classList.remove("over");
  }
}

function commitFlow(payload) {
  if (submitting) return;
  pendingOutput = null;
  picked = dragging = "";
  submitting = true;
  updateSockets();
  report(payload);
}

function snapPlacement(port) {
  const parent = port.closest(".flow-node");
  const top = parent.offsetTop + parent.offsetHeight + 46;
  let left = parent.offsetLeft;
  const blocked = () => [...flow.nodes.values()].some(node =>
    node.offsetTop < top + 125 && node.offsetTop + node.offsetHeight + 10 > top
    && node.offsetLeft < left + 222 && node.offsetLeft + node.offsetWidth + 10 > left);
  while (blocked()) left += parent.offsetWidth + 42;
  return coordinates(left, top);
}

function flowPlacement(event, card) {
  if (!flow || flow.model.managed) return null;
  const element = document.elementFromPoint(event.clientX, event.clientY);
  if (!flow.viewport.contains(element)) return null;
  const port = element.closest('.flow-port[data-direction="out"]');
  const snap = port && canvasMatches(port, card).length === 1;
  if (!snap && element.closest(".flow-node")) return null;
  const {left, top} = flow.plane.getBoundingClientRect();
  const position = snap ? snapPlacement(port)
    : coordinates(event.clientX - left - 12, event.clientY - top - 12);
  return {action: "place", card, position, ...(snap && {source: port.dataset.id, output: port.dataset.port})};
}

function highlightPaletteDrop(event) {
  if (!flow) return;
  const output = document.elementFromPoint(event.clientX, event.clientY)?.closest(".flow-port.available");
  for (const port of flow.ports.values()) port.classList.toggle("over", port === output);
}

function pointerGesture(button, start, move, end) {
  let motion = null, moved = false;
  button.onpointerdown = event => {
    if (submitting || !event.isPrimary || event.button !== 0) return;
    motion = {x: event.clientX, y: event.clientY, id: event.pointerId};
    moved = false;
    button.setPointerCapture(event.pointerId);
  };
  button.onpointermove = event => {
    if (motion?.id !== event.pointerId) return;
    const dx = event.clientX - motion.x, dy = event.clientY - motion.y;
    if (!moved) {
      if (Math.hypot(dx, dy) < 7) return;
      moved = true;
      start();
    }
    move(event, dx, dy);
  };
  function stop(event, cancelled) {
    if (motion?.id !== event.pointerId) return;
    motion = null;
    if (moved) end(event, cancelled);
  }
  button.onpointerup = event => stop(event, false);
  button.onpointercancel = button.onlostpointercapture = event => stop(event, true);
  return callback => { if (!moved) callback(); moved = false; };
}

function connectPorts(output, input) {
  if (!output || !input || !canWire(output, input)) return;
  commitFlow({action: "wire", source: output.dataset.id, output: output.dataset.port,
    target: input.dataset.id, input: input.dataset.port});
}

function flowPort(id, direction, spec) {
  const editable = !flow.model.managed;
  const port = element(editable ? "button" : "span", "flow-port", spec.label);
  if (editable) port.type = "button";
  Object.assign(port.dataset, {id, direction, port: spec.id, multiple: Boolean(spec.multiple)});
  const name = args.program.blocks[id].label;
  port.setAttribute("aria-label", spec.label + (direction === "out" ? " output from " : " input to ") + name);
  flow.ports.set(portKey(id, direction, spec.id), port);
  if (!editable) return port;
  if (direction === "in") {
    port.onclick = () => connectPorts(pendingOutput, port);
    return port;
  }
  function activate(output) {
    picked = "";
    pendingOutput = output;
    updateSockets();
  }
  let hovered = null;
  const click = pointerGesture(port, () => activate(port), event => {
    const candidate = document.elementFromPoint(event.clientX, event.clientY)?.closest('.flow-port[data-direction="in"]');
    hovered?.classList.remove("over");
    hovered = candidate && canWire(port, candidate) ? candidate : null;
    hovered?.classList.add("over");
    drawWires(event);
  }, (_, cancelled) => {
    hovered?.classList.remove("over");
    if (!cancelled && hovered) connectPorts(port, hovered);
    activate(null);
    drawWires();
  });
  port.onclick = () => click(() => {
    if (picked && canvasMatches(port, picked).length === 1) {
      commitFlow({action: "place", card: picked, source: id, output: spec.id,
        position: snapPlacement(port)});
    } else activate(pendingOutput === port ? null : port);
  });
  return port;
}

function flowNode(id, view) {
  const block = args.program.blocks[id];
  const node = element("article", "flow-node " + block.role);
  node.dataset.id = id;
  const pos = displayPosition(view.position, flow.model.managed);
  node.style.left = pos.left + flow.origin.left + "px";
  node.style.top = pos.top + flow.origin.top + "px";
  const input = element("div", "flow-ports inputs");
  const output = element("div", "flow-ports outputs");
  view.inputs.forEach(spec => input.append(flowPort(id, "in", spec)));
  view.outputs.forEach(spec => output.append(flowPort(id, "out", spec)));
  const title = element("button", "flow-step");
  title.type = "button";
  title.dataset.id = id;
  title.append(element("span", "verb", block.action), element("span", "name", block.label));
  if (!block.enabled) title.append(element("span", "hint", "Disabled"));
  if (view.inputs.length) node.append(input);
  node.append(title);
  if (view.outputs.length) node.append(output);
  title.onclick = () => selectBlock(id);
  if (!flow.model.managed) {
    let start;
    const click = pointerGesture(title, () => {
      start = [node.offsetLeft, node.offsetTop];
    }, (_, dx, dy) => {
      node.style.left = Math.max(24, start[0] + dx) + "px";
      node.style.top = Math.max(24, start[1] + dy) + "px";
      drawWires();
    }, (_, cancelled) => {
      if (cancelled) {
        node.style.left = start[0] + "px";
        node.style.top = start[1] + "px";
        drawWires();
      } else {
        commitFlow({action: "move", id, position: coordinates(node.offsetLeft, node.offsetTop)});
      }
    });
    title.onclick = () => click(() => selectBlock(id));
  }
  flow.nodes.set(id, node);
  return node;
}

function renderFlow(program) {
  const model = program.canvas;
  const section = element("section", model.managed ? "flow managed" : "flow");
  const head = element("div", "flow-head");
  const placeButton = element("button", "flow-place");
  placeButton.type = "button";
  placeButton.hidden = true;
  const hint = model.managed ? "Dashed wire returns the next batch to Docking." : "Connect an output to an input.";
  head.append(element("h2", "label", "Canvas"), element("span", "hint", hint), placeButton);
  const viewport = element("div", "flow-viewport");
  viewport.setAttribute("aria-label", "Pipeline canvas");
  viewport.setAttribute("role", "region");
  viewport.tabIndex = 0;
  const plane = element("div", "flow-plane");
  const canvas = element("canvas", "flow-wires");
  canvas.setAttribute("aria-hidden", "true");
  const nodes = Object.entries(model.nodes);
  const positions = nodes.map(([, node]) => displayPosition(node.position, model.managed));
  const origin = flowOrigin ??= {
    left: 24 - (positions.length ? Math.min(...positions.map(pos => pos.left)) : 0),
    top: 24 - (positions.length ? Math.min(...positions.map(pos => pos.top)) : 0),
  };
  const ancestors = {};
  for (const stage of program.stages) {
    for (const id of stage) {
      ancestors[id] = new Set();
      for (const edge of model.edges) {
        if (edge.target !== id || edge.loop_back) continue;
        ancestors[id].add(edge.source);
        for (const parent of ancestors[edge.source] || []) ancestors[id].add(parent);
      }
    }
  }
  flow = {viewport, plane, canvas, placeButton, origin, model, ancestors,
    edges: model.edges, nodes: new Map(), ports: new Map()};
  plane.style.width = Math.max(320, ...positions.map(pos => pos.left + origin.left + 242)) + "px";
  plane.style.height = Math.max(410, ...positions.map(pos => pos.top + origin.top + 200)) + "px";
  plane.append(canvas);
  for (const [id, descriptor] of nodes) plane.append(flowNode(id, descriptor));
  for (const edge of model.edges) flow.ports.get(portKey(edge.target, "in", edge.targetHandle || "in"))?.classList.add("wired");
  if (!nodes.length) plane.append(element("p", "flow-empty", "Drop your first block here"));
  viewport.append(plane);
  section.append(head, viewport);
  viewport.onscroll = () => { if (flow?.viewport === viewport) drawWires(); };
  placeButton.onclick = () => {
    if (!picked) return;
    const top = Math.max(24, ...[...flow.nodes.values()].map(node => node.offsetTop + node.offsetHeight + 42));
    commitFlow({action: "place", card: picked, position: coordinates(24, top)});
  };
  section.onkeydown = event => {
    if (event.key !== "Escape") return;
    picked = "";
    pendingOutput = null;
    updateSockets();
  };
  return section;
}

window.addEventListener("resize", () => { if (flow) drawWires(); });
