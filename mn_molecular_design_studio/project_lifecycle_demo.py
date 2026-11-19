"""Explicitly disposable data for testing project actions in a private preview."""

import json
from pathlib import Path
import uuid

from mn_molecular_design_studio import molecule_imports
from mn_molecular_design_studio.cluster import config
from mn_molecular_design_studio.molecular_database import MolecularDatabase
from mn_molecular_design_studio.project_lifecycle import ProjectCoordinator, ProjectNode, stamp


def demo_workspace(root, user, *, add=False):
    root = Path(root).resolve()
    root.mkdir(parents=True, exist_ok=True, mode=0o700)
    organization_path = root / 'organization/organization.json'
    if not organization_path.exists() or add:
        organization_path.parent.mkdir(parents=True, exist_ok=True, mode=0o700)
        organization = json.loads(organization_path.read_text()) if organization_path.exists() else {'name': 'Disposable cleanup sandbox', 'projects': []}
        suffix = uuid.uuid4().hex[:8]
        project = {'id': uuid.uuid4().hex, 'name': 'Disposable cleanup demo ' + suffix,
                   'description': 'SIMULATED DATA ONLY. Two finished tool executions and an explicitly saved collection.',
                   'created_by': user['id'], 'created_at': stamp()}
        organization['projects'].append(project)
        shared = {'id': uuid.uuid4().hex, 'name': 'Protected shared-data demo ' + suffix,
                  'description': 'SIMULATED DATA ONLY. Shares a chemical identity, not a run folder.',
                  'created_by': user['id'], 'created_at': stamp()}
        organization['projects'].append(shared)
        config.atomic_json(organization_path, organization)
        database = MolecularDatabase(root / 'database')
        source_reference = 'design/design-' + suffix
        for reference, tool in ((source_reference, 'demo_design'), ('analysis/analysis-' + suffix, 'demo_analysis')):
            directory = root / 'runs' / reference
            directory.mkdir(parents=True)
            metadata = {'project_id': project['id'], 'owner_id': user['id'], 'status': 'completed', 'tool': tool,
                        'campaign_name': 'SIMULATED ' + tool, 'created_at': stamp(), 'updated_at': stamp(),
                        'studio_run_id': 'simulated-run-' + suffix, 'studio_run_name': 'Simulated design + post-analysis'}
            config.atomic_json(directory / 'metadata.json', metadata)
            config.atomic_json(directory / 'input.json', {'source_job': source_reference} if tool == 'demo_analysis' else {})
            (directory / 'result.txt').write_text('SIMULATED DATA: not an experimental or design result.\n')
            imported = database.import_records(user, molecule_imports.parse_sequences('>Simulated molecule\nACDEFG'),
                name='SIMULATED ' + tool, project_id=project['id'])
            task_group, run_id = reference.split('/')
            with database.connection() as connection:
                connection.execute('UPDATE molecule_groups SET origin=? WHERE id=?',
                    (json.dumps({'type': 'job', 'task_group': task_group, 'run_id': run_id}), imported['group_id']))
            database.add_observation(imported['group_id'], imported['molecule_ids'][0], tool, run_id, {'simulated_score': .5},
                                     context={'task_group': task_group, 'run_id': run_id})
        database.import_records(user, molecule_imports.parse_sequences('ACDEFG'), name='SIMULATED saved selection', project_id=project['id'])
        database.import_records(user, molecule_imports.parse_sequences('ACDEFG'), name='SIMULATED shared molecule', project_id=shared['id'])
    database = MolecularDatabase(root / 'database')

    def jobs():
        from mn_molecular_design_studio.job_dependencies import local_row
        rows = []
        for path in (root / 'runs').glob('*/*/metadata.json'):
            row = local_row(path.parent, node_id='sandbox', root=root / 'runs')
            rows.append({**row, 'node_name': 'Disposable sandbox'})
        return rows

    node = ProjectNode(home=root, runs=root / 'runs', database=database, node_id='sandbox', job_reader=jobs)
    return node, ProjectCoordinator(node, peers=[]), jobs


def add_cancellation_demo(root, user):
    root = Path(root).resolve()
    path = root / 'organization/organization.json'
    organization = json.loads(path.read_text())
    suffix = uuid.uuid4().hex[:8]
    project = {'id': uuid.uuid4().hex, 'name': 'Disposable cancellation demo ' + suffix,
               'description': 'SIMULATED ONLY. A design, downstream refolding and analysis, a completed result and an independent branch.',
               'created_by': user['id'], 'created_at': stamp()}
    organization['projects'].append(project)
    config.atomic_json(path, organization)
    refs = {'design': 'design/design-' + suffix, 'refolding': 'refolding/refolding-' + suffix,
            'analysis': 'analysis/analysis-' + suffix, 'completed': 'prediction/completed-' + suffix,
            'independent': 'design/independent-' + suffix}
    for stage, status, source in (('design', 'running', None), ('refolding', 'queued', 'design'),
                                   ('analysis', 'queued', 'refolding'), ('completed', 'completed', 'design'),
                                   ('independent', 'queued', None)):
        run = root / 'runs' / refs[stage]
        run.mkdir(parents=True)
        config.atomic_json(run / 'metadata.json', {'project_id': project['id'], 'owner_id': user['id'],
            'status': status, 'tool': 'demo_' + stage, 'campaign_name': 'SIMULATED ' + stage,
            'created_at': stamp(), 'updated_at': stamp(), 'simulated': True})
        config.atomic_json(run / 'input.json', {'source_job': refs[source]} if source else {})
        (run / 'result.txt').write_text('SIMULATED — no scientific computation or worker was started.\n')
    return project['id']


def cancel_simulated_job(node, user, row):
    from mn_molecular_design_studio.job_dependencies import descendants, identity, scope, TERMINAL
    if node.node_id != 'sandbox' or node.runs != node.home / 'runs':
        raise ValueError('Cancellation simulation requires the disposable sandbox.')
    current = next((value for value in node.job_reader() if identity(value) == identity(row)), None)
    if not current or current.get('owner_id') != user['id'] or scope(current) != scope(row):
        raise ValueError('Simulated job is unavailable.')
    selected = descendants(node.job_reader(), [current])
    changed = []
    for value in selected:
        if str(value.get('status') or '').lower() in TERMINAL:
            continue
        path = node.run_path(value['task_group'] + '/' + value['run_id']) / 'metadata.json'
        metadata = json.loads(path.read_text())
        if not metadata.get('simulated'):
            raise ValueError('This job is not part of a simulated cancellation fixture.')
        metadata.update(status='cancelled', cancelled_at=stamp(), updated_at=stamp(), cancel_reason='SIMULATED upstream cancellation')
        config.atomic_json(path, metadata)
        changed.append(identity(value))
    return {'cancelled_jobs': changed, 'downstream_count': sum(value != identity(current) for value in changed)}
