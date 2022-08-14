from __future__ import print_function
from dis import dis
import mdtraj as md
import numpy as np
import sys
import subprocess
import math
import MDAnalysis as mda
import pandas as pd

NATIVE_CUTOFF = 0.5 
traj_name='md_red.xtc'
traj=md.load(traj_name,top='start.pdb')
arr = np.empty((1))
#ref='start.pdb'
#refpdb=md.load(ref)
top_pdb=traj.topology
frames=traj.n_frames
#print(frames)
id_prot=top_pdb.select("protein and (element O or element P or element C or element N or element S)")
rna=[]
id_rna=top_pdb.select("not protein and (element O or element P or element C or element N or element S)")
rna.append(id_rna)
#print(rna)
dfrna=pd.DataFrame(rna).T
#print(dfrna)
#print(id_rna)
    #all pairs of protein-RNA index

frame_list=list(range(0,frames))
res=[]
for i in id_rna:
    atom = top_pdb.atom(i)
    res.append(atom.residue.name)
    #print(atom.residue.name)
dfres=pd.DataFrame(res)
#print(dfres)
df_resrna=[dfrna,dfres]
rna_res=pd.concat(df_resrna, axis=1, ignore_index=True)
#print(rna_res)

for i in frame_list:
    pairs = np.array(
        [(i,j) for  i in id_prot  for j in id_rna])
    #print(pairs)
    
    heavy_pairs_distances = md.compute_distances(traj, pairs)[[frame_list]]
    
    df4=pd.DataFrame(res)
    df1 = pd.DataFrame(pairs)
    distance=np.array(heavy_pairs_distances)
    df2 = pd.DataFrame(distance.T)
    dfs=[df1, df2]
    all = pd.concat(dfs, axis=1, ignore_index=True)




frame1=[]
ntry=np.arange(2,(frames+2))
for i in ntry:
    idxmin = (all.groupby([1])[i]).idxmin()
    #print(idxmin)
    new_df=all.loc[idxmin]
    
    frame1.append(new_df[i].values)

#np.savetxt('check_res.dat', all, fmt='%10.2f')
#print(new_df.shape)
#np.savetxt('min_dist.dat', frame1, fmt='%10.2f')

dist=pd.DataFrame(new_df[i].values)

ind_res_rna=np.arange(len(id_rna))
ind_res_prot=np.arange(len(id_prot))
for i in range(0,len(id_prot)):
            ind_res_prot[i]=top_pdb.atom(id_prot[i]).residue.index
            
for i in range(0,len(id_rna)):
            ind_res_rna[i]=top_pdb.atom(id_rna[i]).residue.index
            

for i in ind_res_rna:
    newval= pd.DataFrame(ind_res_rna - ind_res_prot[-1])

frame1= pd.DataFrame(frame1).T

for i in frame1:
    last_df2=[newval,frame1]

all3 = pd.concat(last_df2, axis=1, ignore_index=True)



for i in all3:
    idxmin = (all3.groupby([0])[i]).idxmin()
    #print(idxmin)
    new_df3=all3.loc[idxmin]
#print(new_df3.T)
#print(new_df3)
#last_df=new_df3.drop(columns=new_df3.columns[0], axis=1, inplace=True)
#print(last_df)
np.savetxt('min_dist_last.dat', new_df3.T, fmt='%10.2f')
#print(new_df3)



for i in all3:
    stack=new_df3[1].append(new_df3[i]).reset_index(drop=True)
#print(stack)

np.savetxt('dist_stack.dat', stack, fmt='%10.2f')
#file=open('amp_trj.dat')
#pucker= pd.DataFrame(file)
#print(pucker)