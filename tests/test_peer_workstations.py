"""Independent installations, simulated network partitions, and rejoining peers.

These tests use real local account databases, synchronization stores and HTTP
routes. Scientific execution is replaced with a recorded submission.
"""
from __future__ import annotations

import copy
import hashlib
import json
import uuid
from contextlib import contextmanager
from pathlib import Path

import pytest
from fastapi.testclient import TestClient

from mn_molecular_design_studio import auth, workspace
from mn_molecular_design_studio.cluster import config, service, transport
from mn_molecular_design_studio.cluster.agent import Agent
from mn_molecular_design_studio.cluster.replication import Replica, capture_identity, sync_peer
from mn_molecular_design_studio.cluster.store import Store

PASSWORD = 'a-long-local-password'
NEW_PASSWORD = 'another-local-password'
SNAPSHOT = {'accept_jobs':True, 'scheduler_online':True, 'storage_available':True,
            'images':{},'references':{},'code_revision':'test',
            'resources':{'free_gpus':[],'cpu_free':16,'ram_available_gb':64,'scratch_free_gb':100}}
NEEDS = {'gpu_count':0,'cpu_threads':1,'code_revision':'test','images':{},'references':[]}


@pytest.fixture
def peers(monkeypatch, tmp_path):
    nodes = [{'id':str(uuid.uuid4()),'name':name,'token':uuid.uuid4().hex+uuid.uuid4().hex,'url':f'http://node-{i}.test:8766'} for i,name in enumerate(('pc002975','4032','4033'))]

    @contextmanager
    def on(index):
        root = tmp_path/str(index)
        with monkeypatch.context() as mp:
            for key,value in {
                'MDS_CLUSTER_HOME':root/'cluster', 'MN_MOLECULAR_DESIGN_STUDIO_HOME':root/'studio',
                'MN_MOLECULAR_DESIGN_STUDIO_APP_HOME':root/'studio','MN_MOLECULAR_DESIGN_STUDIO_RUN_DIR':root/'runs',
                'MN_MOLECULAR_DESIGN_STUDIO_REFERENCE_DIR':root/'references',
                'MN_MOLECULAR_DESIGN_STUDIO_RESULTS_PARENT':root/'results',
                'MN_MOLECULAR_DESIGN_STUDIO_AUTH_MODE':'production',
            }.items():
                mp.setenv(key,str(value))
            yield
    for i,node in enumerate(nodes):
        with on(i):
            workspace.ensure_studio_workspace()
            (config.home()/'node-id').write_text(node['id'])
            config.save({'mode':'peer','name':node['name'],'node_token':node['token'],'coordinator_url':node['url'],
                         'default_placement':'auto','peers':[{'id':p['id'],'name':p['name'],'url':p['url']} for p in nodes if p!=node]})
            with Store().connect() as db:
                for p in nodes:
                    db.execute('INSERT INTO nodes(id,name,token_hash) VALUES (?,?,?)',(p['id'],p['name'],hashlib.sha256(p['token'].encode()).hexdigest()))
            Store().heartbeat(node['id'],copy.deepcopy(SNAPSHOT),[])
    return nodes,on


def transfer(on, source, target):
    with on(source):
        replica=Replica()
        pages=[]; since=0
        while True:
            page=replica.export(since)
            pages.extend(page['records']); since=page['cursor']
            if not page['more']:
                break
    with on(target):
        for index in range(0,len(pages),500):
            Replica().receive(pages[index:index+500])
        Replica().retry_conflicts()
    return pages


def account(on):
    with on(0):
        user=auth.create_initial_admin('scientist','Scientist',PASSWORD)
        project=workspace.create_project('Experiment',owner_id=user['id'])
    return user,project


def network(peers, monkeypatch, offline):
    nodes,on=peers
    from mn_molecular_design_studio import api
    def request(path,payload=None,*,base_url=None,token='',method='GET',raw=False,**kwargs):
        address=base_url or transport._command_server.get() or config.load()['coordinator_url']
        index=next(i for i,p in enumerate(nodes) if p['url']==address)
        if index in offline:
            raise transport.ClusterUnavailable('Simulated disconnected workstation')
        with on(index):
            with TestClient(api.app) as client:
                response=client.request(method,path,json=payload,headers={'Authorization':'Bearer '+token})
        if response.status_code>=400:
            if response.status_code in {401,403,404,409,422}:
                raise ValueError(response.text)
            raise transport.ClusterUnavailable(response.text)
        return response.content if raw else response.json()
    monkeypatch.setattr(transport,'request',request)


def test_login_project_creation_and_local_submission_work_without_any_peer(peers,monkeypatch):
    nodes,on=peers
    user,project=account(on)
    transfer(on,0,1)
    def no_network(*a,**k):
        raise AssertionError('Local operation attempted to use a remote service')
    monkeypatch.setattr(transport,'request',no_network)
    monkeypatch.setattr(service,'schema_check',lambda *a,**k:None)
    from mn_molecular_design_studio.cluster import capabilities
    monkeypatch.setattr(capabilities,'requirements',lambda *a,**k:NEEDS)
    with on(1):
        logged_in,token,_=auth.authenticate('scientist',PASSWORD)
        assert logged_in['id']==user['id']
        assert auth.resolve_session(token)['id']==user['id']
        local_project=workspace.create_project('Created while disconnected',owner_id=user['id'])
        payload={'schema_version':1,'workflow':'analysis.ipsae','parameters':{},'project_id':local_project['id'],
                 'execution':{'placement':'local','idempotency_key':'offline-request'}}
        result=service.submit_from_api(payload,logged_in,request_dir=config.home(),token=token,authorize_path=lambda p:None)
        repeated=service.submit_from_api(payload,logged_in,request_dir=config.home(),token=token,authorize_path=lambda p:None)
        assert result['job']['cluster_job_id']==repeated['job']['cluster_job_id']
        assert result['job']['node_id']==nodes[1]['id']
        assert workspace.validate_project_id(project['id'])==project['id']


def test_two_peers_exchange_projects_and_accounts_while_primary_is_offline(peers,monkeypatch):
    nodes,on=peers
    user,project=account(on)
    transfer(on,0,1); transfer(on,0,2)
    with on(2):
        _,old_session,_=auth.authenticate('scientist',PASSWORD)
    with on(1):
        project2=workspace.create_project('Created on 4032',owner_id=user['id'])
        auth.change_password(user['id'],NEW_PASSWORD,current_password=PASSWORD)
    network(peers,monkeypatch,{0})
    with on(1):
        sync_peer(config.load()['peers'][1],copy.deepcopy(SNAPSHOT))
    with on(2):
        assert auth.authenticate('scientist',NEW_PASSWORD)[0]['id']==user['id']
        assert auth.resolve_session(old_session) is None
        assert workspace.validate_project_id(project2['id'])==project2['id']
        assert {n['name']:n['availability'] for n in Store().nodes()}['4032']=='online'
    # The original workstation catches up from either surviving peer.
    transfer(on,2,0)
    with on(0):
        assert auth.authenticate('scientist',NEW_PASSWORD)[0]['id']==user['id']
        assert workspace.validate_project_id(project2['id'])==project2['id']


def test_concurrent_project_edits_converge_without_using_wall_clock(peers):
    _,on=peers
    user,project=account(on)
    transfer(on,0,1); transfer(on,0,2)
    for i in (1,2):
        with on(i):
            org=workspace.load_organization()
            org['projects'][0]['description']=f'edited on {i}'
            workspace.save_organization(org)
    transfer(on,1,2); transfer(on,2,1); transfer(on,2,0)
    descriptions=[]
    for i in range(3):
        with on(i):
            descriptions.append(workspace.load_organization()['projects'][0]['description'])
    assert len(set(descriptions))==1


def test_membership_revocations_sync_as_tombstones_and_are_not_resurrected(peers):
    _,on=peers
    admin,project=account(on)
    with on(0):
        member=auth.create_user('member','Member',PASSWORD,actor=admin)
        auth.grant_project_access(project['id'],member['id'],actor=admin)
    stale=transfer(on,0,1); transfer(on,0,2)
    with on(1):
        auth.revoke_project_access(project['id'],member['id'],actor=admin)
    transfer(on,1,2); transfer(on,2,0)
    for i in range(3):
        with on(i):
            Replica().receive(stale)
            assert auth.accessible_project_ids(member)==[]


def test_owner_execution_state_outranks_old_origin_placeholders(peers):
    nodes,on=peers
    admin,project=account(on)
    with on(0):
        store=Store(); store.heartbeat(nodes[1]['id'],copy.deepcopy(SNAPSHOT),[])
        request={'kind':'workflow','payload':{'schema_version':1,'workflow':'analysis.ipsae','parameters':{}},'project_id':project['id'],'bundle_id':None}
        job=store.submit(request,NEEDS,admin,origin_node=nodes[0]['id'],placement=nodes[1]['id'])['job']
        Replica().capture_jobs()
    stale=transfer(on,0,1)
    with on(1):
        Store().heartbeat(nodes[1]['id'],copy.deepcopy(SNAPSHOT),[{'cluster_job_id':job['cluster_job_id'],'task_group':'design','run_id':'success',
            'owner_id':admin['id'],'project_id':project['id'],'status':'completed','metrics':{'test_score':1},'run_dir':'/node1/results/design/success'}])
        Replica().capture_jobs()
    transfer(on,1,2); transfer(on,2,0)
    with on(0):
        result=Store().job(job['cluster_job_id'],admin)
        assert result['status']=='completed'
        assert result['node_id']==nodes[1]['id']
        assert result['run_dir'].startswith('cluster://')
        assert result['location']=='/node1/results/design/success'
    with on(2):
        Replica().receive(stale)
        assert Store().job(job['cluster_job_id'],admin)['status']=='completed'
        exported=next(r for r in Replica().export()['records'] if r['kind']=='job')
        spoof=copy.deepcopy(exported); spoof['author']=nodes[2]['id']; spoof['counter']+=100
        Replica().receive([spoof])
        assert Store().job(job['cluster_job_id'],admin)['node_id']==nodes[1]['id']
        with Replica().db() as db:
            assert db.execute('SELECT count(*) FROM conflicts').fetchone()[0]==1


def test_agent_starts_when_every_other_machine_and_local_api_are_unavailable(peers,monkeypatch):
    _,on=peers
    def offline(*a,**k):
        raise transport.ClusterUnavailable('offline')
    monkeypatch.setattr(transport,'request',offline)
    monkeypatch.setattr(Agent,'snapshot',lambda self:copy.deepcopy(SNAPSHOT))
    monkeypatch.setattr(Agent,'catalog',lambda self:[])
    monkeypatch.setattr(Agent,'storage_online',lambda self:True)
    with on(1):
        Agent().run(once=True)
        assert Store().nodes()[1:]
        assert json.loads((config.home()/'agent-status.json').read_text())['online']


def test_peer_agent_processes_local_queue_and_replays_receipt_without_primary(peers,monkeypatch):
    nodes,on=peers
    admin,project=account(on)
    transfer(on,0,1); transfer(on,0,2)
    network(peers,monkeypatch,{0})
    calls=[]
    def submit(self,command):
        calls.append(command['id'])
        return {'job':{'task_group':'design','run_id':'recorded-submission','run_dir':'/local/results/design/recorded-submission','status':'queued'}}
    monkeypatch.setattr(Agent,'submit',submit)
    with on(1):
        request={'kind':'workflow','payload':{'workflow':'analysis.ipsae'},'project_id':project['id'],'bundle_id':None}
        store=Store()
        result=store.submit(request,NEEDS,admin,origin_node=nodes[1]['id'],placement='local')
        agent=Agent()
        assert agent.poll()
        assert len(calls)==1
        command_id=calls[0]
        # Even a replay after losing a response uses the persisted local receipt.
        response,error=Agent().execute({'id':command_id,'kind':'workflow','job_id':result['job']['cluster_job_id']})
        assert error is None and response['job']['run_id']=='recorded-submission'
        assert len(calls)==1


def test_peer_metadata_routes_require_machine_credentials_and_keep_user_endpoints_private(peers):
    nodes,on=peers
    from mn_molecular_design_studio import api
    with on(1),TestClient(api.app) as client:
        path='/api/v1/cluster/node/peer/records'
        assert client.get(path).status_code==401
        assert client.get(path,headers={'Authorization':'Bearer invented'}).status_code==403
        headers={'Authorization':'Bearer '+nodes[2]['token']}
        assert client.get(path,headers=headers).status_code==200
        assert client.get('/api/v1/cluster/jobs',headers=headers).status_code==401
        response=client.post('/api/v1/cluster/node/peer/hello',headers=headers,json={'node_id':nodes[0]['id'],'snapshot':SNAPSHOT})
        assert response.status_code==403


def test_conflicting_new_usernames_are_preserved_for_admin_review(peers):
    _,on=peers
    account(on)
    transfer(on,0,1); transfer(on,0,2)
    users=[]
    for i in (1,2):
        with on(i):
            admin=auth.list_users()[0]
            users.append(auth.create_user('same-new-name',f'Scientist {i}',PASSWORD,actor=admin))
    transfer(on,1,2)
    with on(2):
        assert next(u for u in auth.list_users() if u['username']=='same-new-name')['id']==users[1]['id']
        with Replica().db() as db:
            assert db.execute('SELECT count(*) FROM conflicts').fetchone()[0]==1


def test_missing_results_storage_never_dispatches_a_local_peer_job(peers):
    nodes,on=peers
    user,project=account(on)
    with on(1):
        Store().heartbeat(nodes[1]['id'],dict(SNAPSHOT,accept_jobs=False,storage_available=False),[])
        request={'kind':'workflow','payload':{'workflow':'analysis.ipsae'},'project_id':project['id'],'bundle_id':None}
        job=Store().submit(request,NEEDS,user,origin_node=nodes[1]['id'],placement='local')['job']
        assert job['status']=='waiting_for_workstation'
        assert Store().claim(nodes[1]['id']) is None


def test_unclaimed_automatic_work_can_move_but_claimed_work_remains_pinned(peers):
    nodes,on=peers
    user,project=account(on)
    with on(1):
        store=Store()
        store.heartbeat(nodes[0]['id'],copy.deepcopy(SNAPSHOT),[])
        busy=copy.deepcopy(SNAPSHOT); busy['resources']['cpu_free']=0
        store.heartbeat(nodes[1]['id'],busy,[])
        request={'kind':'workflow','payload':{'workflow':'analysis.ipsae'},'project_id':project['id'],'bundle_id':None}
        first=store.submit(request,NEEDS,user,origin_node=nodes[1]['id'])['job']
        assert first['node_id']==nodes[0]['id']
        with store.connect() as db:
            db.execute('UPDATE nodes SET last_seen=0 WHERE id=?',(nodes[0]['id'],))
        store.heartbeat(nodes[2]['id'],copy.deepcopy(SNAPSHOT),[])
        assert store.job(first['cluster_job_id'],user)['node_id']==nodes[2]['id']
        assert store.claim(nodes[0]['id']) is None
        claimed=store.claim(nodes[2]['id'])
        assert claimed['job_id']==first['cluster_job_id']
        with store.connect() as db:
            db.execute('UPDATE nodes SET last_seen=0 WHERE id=?',(nodes[2]['id'],))
        store.heartbeat(nodes[0]['id'],copy.deepcopy(SNAPSHOT),[])
        assert store.job(first['cluster_job_id'],user)['node_id']==nodes[2]['id']
        assert store.claim(nodes[0]['id']) is None


def test_incremental_sync_catches_up_after_disconnection_without_duplicate_records(peers,monkeypatch):
    nodes,on=peers
    user,project=account(on)
    network(peers,monkeypatch,set())
    with on(0):
        peer=config.load()['peers'][0]
        sync_peer(peer,copy.deepcopy(SNAPSHOT))
        first_cursor=Replica().state('peer:'+nodes[1]['id'])['push']
        new=workspace.create_project('Later project',owner_id=user['id'])
        assert Replica().export(first_cursor)['records']
        sync_peer(peer,copy.deepcopy(SNAPSHOT))
        cursor=Replica().state('peer:'+nodes[1]['id'])['push']
        sync_peer(peer,copy.deepcopy(SNAPSHOT))
        assert Replica().state('peer:'+nodes[1]['id'])['push']==cursor
    with on(1):
        assert workspace.validate_project_id(new['id'])==new['id']
        with Replica().db() as db:
            assert db.execute('SELECT count(*) FROM records WHERE key=?',('project:'+new['id'],)).fetchone()[0]==1


def test_migration_preserves_identity_and_backs_up_accounts_and_registry(peers,tmp_path):
    nodes,on=peers
    from typer.testing import CliRunner
    from mn_molecular_design_studio.cluster.cli import app
    with on(1):
        config.save({'mode':'worker','coordinator_url':'http://127.0.0.1:18766'})
        own=config.identity(); token=config.load()['node_token']
        with Store().connect() as db:
            members=[dict(r) for r in db.execute('SELECT id,name,token_hash,disabled FROM nodes')]
        manifest=tmp_path/'private-manifest.json'
        manifest.write_text(json.dumps({'schema_version':1,'node_id':own,'members':members,'peers':[{'id':p['id'],'name':p['name'],'url':p['url']} for p in nodes if p['id']!=own]}))
        result=CliRunner().invoke(app,['enable-peers',str(manifest)])
        assert result.exit_code==0,result.output
        assert config.identity()==own and config.load()['node_token']==token
        assert config.is_peer() and not config.is_worker()
        assert list(config.home().glob('before-peers-*/coordinator.sqlite3'))
        assert token not in result.output
