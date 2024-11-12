import mdtraj as md
import numpy as np
import os
import pickle as pkl
recalculate = True
ensemble_subsampled_size = 50
def load_pickle(filename):
    with open(filename, 'rb') as f:
        loaded_obj = pkl.load(f)
        
    return loaded_obj

def save_pickle(filename, pickle_obj):
    with open(filename, 'wb') as f:
        pkl.dump(pickle_obj, f)
weights = load_pickle('BME_rew_SAXSonly/weights.pkl')
if recalculate == True:
    frames_sel = np.random.choice(np.arange(0,len(weights)), size=ensemble_subsampled_size, replace=True, p=weights)
os.system('rm -r subsampled_frames_ANAC046 subsampled_frames_ANAC046.tar.gz')
os.system('mkdir subsampled_frames_ANAC046')
save_pickle('subsampled_frames_ANAC046/idx_frames_sel.pkl', frames_sel)
traj = md.load('../Simulations/traj_AA.xtc', top='../Simulations/top_AA.pdb')
superpose_atoms = traj.top.select('resi 73 to 93')
traj= traj.superpose(traj, frame=0, atom_indices=superpose_atoms, ref_atom_indices=superpose_atoms, parallel=True)

for i,frame in enumerate(frames_sel):
    traj[frame].save_pdb(f'subsampled_frames_ANAC046/frame_{i}.pdb')
os.system('tar -zcvf subsampled_frames_ANAC046.tar.gz subsampled_frames_ANAC046')
recalculate = False

