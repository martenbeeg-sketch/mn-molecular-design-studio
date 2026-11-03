const assert = require("node:assert/strict");
const fs = require("node:fs");
const path = require("node:path");
const test = require("node:test");
const ts = require("typescript");

const source = fs.readFileSync(path.join(__dirname, "../src/residueSelection.ts"), "utf8");
const javascript = ts.transpileModule(source, {
  compilerOptions: { module: ts.ModuleKind.CommonJS, target: ts.ScriptTarget.ES2020 },
}).outputText;
const loaded = { exports: {} };
new Function("exports", "require", "module", javascript)(loaded.exports, require, loaded);
const { parseResidueSelection, parseResidueSelections } = loaded.exports;

test("preserves existing author-number range selectors", () => {
  assert.deepEqual(parseResidueSelection("A:10-12"), {
    chainId: "A", residues: [10, 11, 12], insertionCode: null,
  });
});

test("selects an exact author insertion code or the blank insertion code", () => {
  assert.deepEqual(parseResidueSelection("A:10^B"), {
    chainId: "A", residues: [10], insertionCode: "B",
  });
  assert.deepEqual(parseResidueSelection("A:10^_"), {
    chainId: "A", residues: [10], insertionCode: "",
  });
});

test("preserves old selectors and rejects insertion-coded ranges", () => {
  assert.deepEqual(parseResidueSelection("A10"), {
    chainId: "A", residues: [10], insertionCode: null,
  });
  assert.throws(() => parseResidueSelection("A:10^A-12"));
  assert.equal(parseResidueSelections(["A:10^A", "B:7^_"]).length, 2);
});
