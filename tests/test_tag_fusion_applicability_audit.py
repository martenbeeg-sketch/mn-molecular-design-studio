"""AUD-012: one component may have a multi-chain interface partner."""
import pytest

from mn_molecular_design_studio.complex_analysis_batch import assessment_applicability


@pytest.mark.parametrize("binder,target,valid", [
    (["B"], ["A"], True), (["B"], ["A", "C"], True),
    (["B"], ["B"], False), (["B"], [], False), (["B", "C"], ["A"], False),
])
def test_tag_fusion_matches_component_and_partner_assembly_contract(binder, target, valid):
    assert assessment_applicability({"kind": "protein-protein"}, "analysis.tag_fusion_suitability", binder, target)[0] is valid
    if len(target) > 1:
        for method in ("analysis.rosetta", "analysis.ipsae"):
            assert assessment_applicability({"kind": "protein-protein"}, method, binder, target)[0] is False
