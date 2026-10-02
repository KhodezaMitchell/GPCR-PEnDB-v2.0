email = ''
api_key = ''

def esearch(*args, **kwargs):
    raise SystemExit('OFFLINE CHECK: PubMed ESearch attempted; historical cache miss')

def efetch(*args, **kwargs):
    raise SystemExit('OFFLINE CHECK: PubMed EFetch attempted; historical cache miss')

def read(*args, **kwargs):
    raise SystemExit('OFFLINE CHECK: Entrez parser called without cached record')
