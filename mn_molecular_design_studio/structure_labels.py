"""Display labels supported by the structure's producing record."""
from __future__ import annotations

import json
from pathlib import PurePosixPath


BINDCRAFT2_RANKED_LABEL = "BindCraft2 ranked refolded complex"


def is_bindcraft2_ranked_path(path):
    if not isinstance(path, str):
        return False
    source = PurePosixPath(path.replace("\\", "/"))
    return (source.parts[-4:-1] == ("bindcraft2", "output", "3_Ranked")
            and source.suffix.lower() in {".cif", ".mmcif", ".pdb"}
            and not source.stem.endswith("_monomer"))


def ranked_bindcraft2_labels(connection, structure_ids, *, visible_groups=None):
    """Resolve legacy labels within the caller's authorized structure scope."""
    if not structure_ids:
        return {}
    groups = json.dumps(visible_groups) if visible_groups is not None else None
    rows = connection.execute(
        """SELECT s.id, json_extract(o.context,'$.candidate.complex_pdb') AS source
           FROM molecule_structures s
           JOIN molecule_observations o ON o.molecule_id=s.molecule_id AND o.group_id=s.group_id
           WHERE s.id IN (SELECT value FROM json_each(?))
             AND (? IS NULL OR s.group_id IN (SELECT value FROM json_each(?)))
             AND s.label='Complex Pdb' AND json_extract(s.context,'$.tool')='bindcraft2'
             AND o.namespace='bindcraft2' AND json_extract(o.context,'$.representation_id')=s.id
             AND o.run_id=json_extract(s.context,'$.run_id')
             AND json_extract(o.context,'$.candidate_id')=json_extract(s.context,'$.candidate_id')
             AND json_extract(o.context,'$.candidate.stage')='complex_refolding'
             AND json_extract(o.context,'$.structure_digest')=COALESCE(
                 NULLIF(json_extract(s.context,'$.portable_sha256'),''),
                 NULLIF(json_extract(s.context,'$.content_sha256'),''),
                 NULLIF(json_extract(s.context,'$.sha256'),''))""",
        (json.dumps(structure_ids), groups, groups),
    )
    return {row["id"]: BINDCRAFT2_RANKED_LABEL for row in rows if is_bindcraft2_ranked_path(row["source"])}
