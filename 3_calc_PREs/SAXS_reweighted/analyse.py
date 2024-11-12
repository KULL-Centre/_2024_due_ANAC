import pandas as pd
import numpy as np
import mdtraj as md
import itertools

def initProteins():
    proteins = pd.DataFrame(columns=['labels','wh','temp','obs','pH','expPREs','ionic','fasta'])
    #proteins = pd.DataFrame(columns=['labels','wh','temp','obs','pH','ionic','fasta','path'])
    fasta_ANAC046 = """NAPSTTITTTKQLSRIDSLDNIDHLLDFSSLPPLIDPGFLGQPGPSFSGARQQHDLKPVLHHPTTAPVDNTYLPTQALNFPYHSVHNSGSDFGYGAGSGNNNKGMIKLEHSLVSVSQETGLSSDVNTTATPEISSYPMMMNPAMMDGSKSACDGLDDLIFWEDLYTS""".replace('\n', '')
    proteins.loc['ANAC046'] = dict(labels=[45,87,151],
                               wh=600,temp=283,obs='ratio',pH=7.0,fasta=list(fasta_ANAC046),ionic=0.1,expPREs=None)
    return proteins

