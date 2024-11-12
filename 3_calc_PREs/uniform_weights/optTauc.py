from analyse import *
import MDAnalysis
import time
import os
import glob
import sys
from DEERPREdict.PRE import PREpredict
from DEERPREdict.utils import Operations
from argparse import ArgumentParser
from scipy.stats import spearmanr
import logging
import psutil
import shutil

parser = ArgumentParser()
parser.add_argument('--name',dest='name',type=str,required=True)
args = parser.parse_args()

proteins = pd.read_pickle('proteins.pkl').astype(object)
delay=5.4e-3

def loadExpPREs(name,prot):
    value = {}
    error = {}
    resnums = np.arange(1,len(prot.fasta)+1)
    for label in prot.labels:
        value[label], error[label] = np.loadtxt('../expPREs/exp-{:d}.dat'.format(label),unpack=True, usecols=(0,1))
    v = pd.DataFrame(value,index=resnums)
    v.rename_axis('residue', inplace=True)
    v.rename_axis('label', axis='columns',inplace=True)
    e = pd.DataFrame(error,index=resnums)
    e.rename_axis('residue', inplace=True)
    e.rename_axis('label', axis='columns',inplace=True)
    return pd.concat(dict(value=v,error=e),axis=1)

def calcChi2(name,prot):
    obs = 1 if prot.obs=='ratio' else 2
    chi2 = 0
    for label in prot.labels:
        y = np.loadtxt('calcPREs/res-{:d}.dat'.format(label))[:,obs]
        chi = (prot.expPREs.value[label].values - y) / prot.expPREs.error[label].values
        chi = chi[~np.isnan(chi)]
        chi2 += np.nansum( np.power( chi, 2) ) / chi.size
    return chi2 / len(prot.labels)

def loadInitPREs(name,prot):
    obs = 1 if prot.obs=='ratio' else 2
    value = {}
    resnums = np.arange(1,len(prot.fasta)+1)
    for label in prot.labels:
        value[label] = np.loadtxt('calcPREs/res-{:d}.dat'.format(label))[:,obs]
    v = pd.DataFrame(value,index=resnums)
    v.rename_axis('residue', inplace=True)
    v.rename_axis('label', axis='columns',inplace=True)
    return v

def optTauC(name,prot,chi2_results,rs_results,tau_c_range,delay):
    chi2list = []
    for tc in tau_c_range:
        chi2 = 0
        rs = 0
        for label in prot.labels:
            x,y = np.loadtxt('calcPREs/res-{:d}.dat'.format(label),usecols=(0,1),unpack=True)
            measured_resnums = np.where(~np.isnan(y))[0]
            data = pd.read_pickle('calcPREs/res-{:d}.pkl'.format(label), compression='gzip')
            gamma_2_av = np.full(y.size, fill_value=np.NaN)
            s_pre = np.power(data['r3'], 2)/data['r6']*data['angular']
            gamma_2 = Operations.calc_gamma_2(data['r6'], s_pre, tau_c = tc * 1e-9, tau_t = 1e-10, wh = prot.wh, k = 1.23e16)
            gamma_2 = np.ma.MaskedArray(gamma_2, mask = np.isnan(gamma_2))
            gamma_2_av[measured_resnums] = np.ma.average(gamma_2, axis=0).data
            if prot.obs == 'ratio':
                y = 10 * np.exp(-gamma_2_av * delay) / ( 10 + gamma_2_av )
            else:
                y = gamma_2_av
            mask = np.isfinite(y)&np.isfinite(prot.expPREs.value[label].values)
            chi = (prot.expPREs.value[label].values[mask] - y[mask]) / prot.expPREs.error[label].values[mask]
            chi2 += np.nansum( np.power( chi, 2) ) / chi.size
            rs += spearmanr(y[mask],prot.expPREs.value[label].values[mask])[0]
        chi2list.append(chi2 / len(prot.labels))
        chi2_results.loc[tc] = chi2 / len(prot.labels)
        rs_results.loc[tc] = rs / len(prot.labels)
    tc_min = tau_c_range[np.argmin(chi2list)]
    chi2_results.to_pickle('chi2_'+name+'_tc.pkl')
    rs_results.to_pickle('rs_'+name+'_tc.pkl')

    for label in prot.labels:
        x,y = np.loadtxt('calcPREs/res-{:d}.dat'.format(label),usecols=(0,1),unpack=True)
        measured_resnums = np.where(~np.isnan(y))[0]
        data = pd.read_pickle('calcPREs/res-{:d}.pkl'.format(label),compression='gzip')
        gamma_2_av = np.full(y.size, fill_value=np.NaN)
        s_pre = np.power(data['r3'], 2)/data['r6']*data['angular']
        gamma_2 = Operations.calc_gamma_2(data['r6'], s_pre, tau_c = tc_min * 1e-9, tau_t = 1e-10, wh = prot.wh, k = 1.23e16)
        gamma_2 = np.ma.MaskedArray(gamma_2, mask = np.isnan(gamma_2))
        gamma_2_av[measured_resnums] = np.ma.average(gamma_2, axis=0).data
        i_ratio = 10 * np.exp(-gamma_2_av * delay) / ( 10 + gamma_2_av )
        np.savetxt('calcPREs/res-{:d}.dat'.format(label),np.c_[x,i_ratio,gamma_2_av])
    return tc_min, calcChi2(name,prot)

for name in proteins.index:
    proteins.at[name,'expPREs'] = loadExpPREs(name,proteins.loc[name])

tau_c_range = np.arange(1.0, 20.1, 1.0)
chi2_results = pd.DataFrame(index=tau_c_range)
rs_results = pd.DataFrame(index=tau_c_range)

# run PRE calculations
time0 = time.time()
tau_c, chi2_pre = optTauC(args.name,proteins.loc[args.name],chi2_results,rs_results,tau_c_range,delay)
PREs = {}
PREs[args.name] = {'chi2':chi2_pre, 'tau_c':tau_c, 'expPREs':proteins.at[args.name,'expPREs'], 'calcPREs':loadInitPREs(args.name,proteins.loc[args.name])}
print(args.name,tau_c,chi2_pre)
pd.DataFrame(PREs).to_pickle('{:s}_PREs.pkl'.format(args.name))
print('Timing {:.3f}'.format(time.time()-time0))
