import os, sys, gzip, json
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
import importlib
import utils.brenda as brenda

def test_query_brenda(tmp_path, monkeypatch):
    sample = tmp_path/'brenda.json.gz'
    with gzip.open(sample, 'wt') as fh:
        json.dump([
            {"ec_number": "1.1.1.1", "substrate": ["ethanol"], "organism": ["Escherichia coli"], "pathway": ["ethanol degradation"]}
        ], fh)
    monkeypatch.setenv('BRENDA_JSON_PATH', str(sample))
    importlib.reload(brenda)
    res = brenda.query_brenda('ethanol', host='Escherichia coli')
    assert res['found']
    assert '1.1.1.1' in res['ec_numbers']
    assert res['host_matches'] == ['1.1.1.1']