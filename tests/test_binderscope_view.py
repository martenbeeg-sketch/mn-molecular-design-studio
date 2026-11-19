import copy
import hashlib
import json
import re
from pathlib import Path

import gemmi
import pytest

from mn_molecular_design_studio.binderscope_view import build_view_payload, inspector_document


@pytest.fixture
def subject(tmp_path):
    from test_binderscope import write_complex
    from mn_molecular_design_studio.engine.workflow_engine.tools.binderscope import normalize_complex

    path = tmp_path / "design.pdb"
    write_complex(path)
    structure = gemmi.read_structure(str(path))
    structure[0]['X'][0].seqid = gemmi.SeqId(306, 'A')
    structure.write_pdb(str(path))
    mapping, sequences = normalize_complex(path, tmp_path / 'canonical.pdb', 'Y', 'X')
    text = path.read_text()
    count = len(mapping['B'])
    result = {'candidate_id': 'one', 'source': {'source': {'record_name': 'Saved design'},
        'residue_map': mapping, 'input_structure_digest': hashlib.sha256(text.encode()).hexdigest(),
        'assessment': {'b_factor_interpretation': 'plddt_0_100', 'component_format': 'linear protein'}},
        'analysis': {'plddt_res': [80.] * count, 'dsasa_res': [0.] * count, 't_dsasa': {'1': 30.},
            'ss': 'H' * count, 'tc': [1], 'bc': [2], 'contacts': [{'t': 1, 'b': 2, 'd': 4., 'n': 2}],
            'hbonds': [{'t': 1, 'b': 2, 'ta': 'O', 'ba': 'N', 'd': 3.2}], 'salt': [], 'buns': [],
            'liab': [{'kind': 'Example', 'pos': 2, 'reason': 'Example test'}]},
        'hotspot_compliance': {'status': 'assessed', 'residues': [{'canonical_number': 1, 'requested': 'X306A', 'status': 'missed'}]}}
    return result, text


def test_payload_preserves_exact_residue_identity_and_only_renumbers_a_display_copy(subject):
    result, text = subject
    original = copy.deepcopy(result)
    payload = build_view_payload(result, text)
    target = payload['residues']['A'][0]
    assert (target['aa'], target['label'], target['author'], target['number']) == ('A', 'A306A', 'X:306A', 1)
    assert target['confidence'] is None and payload['residues']['B'][0]['confidence'] == 80.
    assert target['burial'] == 30. and payload['residues']['B'][0]['burial'] == 0.
    rendered = gemmi.read_pdb_string(payload['pdb'])
    assert [chain.name for chain in rendered[0]] == ['A', 'B']
    assert rendered[0]['A'][0].seqid.num == 1
    assert '306A' in text and result == original
    assert len(payload['contacts']) == 1 and len(payload['interactions']) == 1
    assert payload['interactions'][0]['atoms'] == ['O', 'N']
    assert payload['hotspots'][0]['id'] == 'A:1'


@pytest.mark.parametrize('problem', ['checksum', 'missing_residue', 'duplicate_number', 'duplicate_author', 'noncontiguous_number', 'wrong_chain', 'extra_residue'])
def test_payload_never_guesses_changed_or_ambiguous_coordinates(subject, problem):
    result, text = subject
    if problem == 'checksum':
        text += '\n'
    elif problem == 'missing_residue':
        result['source']['residue_map']['A'][0]['number'] = 9999
    elif problem == 'duplicate_number':
        result['source']['residue_map']['A'][1]['canonical_number'] = 1
    elif problem == 'duplicate_author':
        result['source']['residue_map']['A'][1]['number'] = result['source']['residue_map']['A'][0]['number']
        result['source']['residue_map']['A'][1]['insertion_code'] = result['source']['residue_map']['A'][0]['insertion_code']
    elif problem == 'noncontiguous_number':
        result['source']['residue_map']['A'][1]['canonical_number'] = 99
    elif problem == 'wrong_chain':
        for row in result['source']['residue_map']['B']:
            row['chain'] = 'absent'
    else:
        result['source']['residue_map']['A'].pop()
    with pytest.raises(ValueError):
        build_view_payload(result, text)


@pytest.mark.parametrize('problem', ['unknown', 'missing', 'incomplete', 'outside_range'])
def test_unknown_confidence_remains_unknown_in_tracks_and_coloring(subject, problem):
    result, text = subject
    if problem == 'unknown':
        result['source']['assessment']['b_factor_interpretation'] = 'uninterpreted'
    elif problem == 'missing':
        result['analysis'].pop('plddt_res')
    elif problem == 'incomplete':
        result['analysis']['plddt_res'].pop()
    else:
        result['analysis']['plddt_res'][0] = 1000.
    payload = build_view_payload(result, text)
    assert not payload['confidence_confirmed']
    assert all(row['confidence'] is None for row in payload['residues']['B'])


def test_native_cif_without_occupancy_is_read_on_a_copy(subject):
    result, text = subject
    structure = gemmi.read_pdb_string(text)
    block = structure.make_mmcif_block()
    category = block.get_mmcif_category('_atom_site.')
    category.pop('occupancy')
    block.set_mmcif_category('_atom_site.', category)
    cif = block.as_string()
    result['source']['input_structure_digest'] = hashlib.sha256(cif.encode()).hexdigest()
    payload = build_view_payload(result, cif)
    assert len(payload['residues']['A']) == len(result['source']['residue_map']['A'])
    assert '_atom_site.occupancy' not in cif


def test_missing_atom_overlay_is_marked_unavailable_not_drawn_between_other_atoms(subject):
    result, text = subject
    result['analysis']['hbonds'][0]['ta'] = 'ABSENT'
    payload = build_view_payload(result, text)
    assert not payload['interactions'] and payload['unmapped_interactions'] == 1


def test_salt_overlay_uses_recorded_charged_group_centroids(subject):
    result, text = subject
    structure = gemmi.read_pdb_string(text)
    target = structure[0]['X'][0]
    target.name = 'ARG'
    for name, coordinates in [('NE', (1., 2., 3.)), ('NH1', (2., 3., 4.)), ('NH2', (3., 4., 5.)), ('CZ', (4., 5., 6.))]:
        atom = gemmi.Atom()
        atom.name = name
        atom.element = gemmi.Element(name[0])
        atom.pos = gemmi.Position(*coordinates)
        target.add_atom(atom)
    binder = structure[0]['Y'][1]
    binder.name = 'GLU'
    for name, coordinates in [('OE1', (5., 6., 7.)), ('OE2', (3., 4., 5.))]:
        atom = gemmi.Atom()
        atom.name = name
        atom.element = gemmi.Element('O')
        atom.pos = gemmi.Position(*coordinates)
        binder.add_atom(atom)
    text = structure.make_pdb_string()
    result['source']['input_structure_digest'] = hashlib.sha256(text.encode()).hexdigest()
    result['analysis']['salt'] = [{'t': 1, 'b': 2, 'tn': 'ARG', 'bn': 'GLU', 'd': 2.6}]
    payload = build_view_payload(result, text)
    interaction = next(row for row in payload['interactions'] if row['kind'] == 'Salt bridge')
    assert interaction['start'] == [2.5, 3.5, 4.5]
    assert interaction['end'] == [4., 5., 6.]
    assert 'centroids' in interaction['basis']


def test_invalid_coordinate_format_has_a_readable_fallback_error(subject):
    result, text = subject
    text = 'data_empty\n'
    result['source']['input_structure_digest'] = hashlib.sha256(text.encode()).hexdigest()
    with pytest.raises(ValueError, match='could not be parsed'):
        build_view_payload(result, text)


def test_document_is_self_contained_escapes_data_and_does_not_rescan_template_tokens(subject):
    result, text = subject
    payload = build_view_payload(result, text, title='</script><script>bad()</script> __STUDIO_SCRIPT__ & hello')
    document = inspector_document(payload)
    data = re.search(r'<script type="application/json" id="assessment-data">(.*?)</script>', document, re.S).group(1)
    assert json.loads(data)['title'] == payload['title']
    assert '<script>bad()' not in document
    assert '__STUDIO_SCRIPT__' in json.loads(data)['title']
    assert not re.search(r'<script\b[^>]*\bsrc=', document)
    assert "connect-src 'none'" in document
    assert document.count('id="contact-diagram"') == 1
    assert 'setHoverable' in document and 'setClickable' in document
    assert 'inspector.js' not in document


def test_native_inspector_serializes_the_verified_model_and_retains_complete_exports(subject, tmp_path):
    from streamlit.testing.v1 import AppTest

    result, text = subject
    path = tmp_path / 'saved-complex.pdb'
    path.write_text(text)
    result['source']['source']['structure_id'] = 'saved-sql-complex'
    app = AppTest.from_string(f'''
from pathlib import Path
from mn_molecular_design_studio.binderscope_inspector import render_inspector
class Database:
    def structure_file(self, user, identifier):
        assert user['id'] == 'authorized' and identifier == 'saved-sql-complex'
        return Path({str(path)!r})
render_inspector({result!r}, database=Database(), user={{'id':'authorized'}})
''').run()
    assert not app.exception
    frames = app.get('iframe')
    assert len(frames) == 1
    document = frames[0].proto.srcdoc
    data = json.loads(re.search(r'id="assessment-data">(.*?)</script>', document, re.S).group(1))
    assert data['source_sha256'] == hashlib.sha256(text.encode()).hexdigest()
    assert len(data['contacts']) == 1
    assert data['residues']['A'][0]['author'] == 'X:306A'
    downloads = {button.proto.label for button in app.get('download_button')}
    assert {'Download assessed complex', 'Download all contacts (CSV)', 'Download this assessment (JSON)'} <= downloads
