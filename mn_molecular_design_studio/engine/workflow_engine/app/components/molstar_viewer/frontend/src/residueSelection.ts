export type ParsedResidueSelection = {
  chainId: string;
  residues: number[];
  insertionCode: string | null;
};

const range = (start: number, end: number) => {
  const step = start <= end ? 1 : -1;
  return Array.from({ length: Math.abs(end - start) + 1 }, (_, i) => start + i * step);
};

/** Parse the viewer's existing chain:number syntax plus exact author insertion codes.
 * `^_` selects the blank insertion code; a selector without `^` retains the
 * legacy behavior of matching every residue carrying that author number.
 */
export function parseResidueSelection(value: string): ParsedResidueSelection {
  const item = String(value || "").trim();
  const auth = item.match(/^([^:\s^]+):(-?\d+)(?:\^([A-Za-z0-9_]))?(?:-(-?\d+))?$/);
  const legacy = auth ? null : item.match(/^([A-Za-z]+)(-?\d+)(?:-(-?\d+))?$/);
  const match = auth || legacy;
  if (!match) throw new Error(`Invalid format for selection: ${item}`);

  const insertionCode = auth && match[3] !== undefined ? (match[3] === "_" ? "" : match[3]) : null;
  const end = auth ? match[4] : match[3];
  if (insertionCode !== null && end !== undefined) {
    throw new Error("Exact insertion-code selections must identify one author residue.");
  }
  const startNumber = parseInt(match[2], 10);
  const residues = end ? range(startNumber, parseInt(end, 10)) : [startNumber];
  return {
    chainId: match[1] === "_" ? "" : match[1],
    residues,
    insertionCode,
  };
}

export function parseResidueSelections(selections: string[]): ParsedResidueSelection[] {
  return selections.map(parseResidueSelection);
}
