# -*- coding: utf-8 -*-
"""
Created on Tue Jun 21 14:38:37 2022

@author: ZR48SA
"""
#%% change directory to script directory (should work on windows and mac)
import os
from pathlib import Path
from inspect import getsourcefile
os.chdir(str(Path(os.path.abspath(getsourcefile(lambda:0))).parents[0]))
script_dir=os.getcwd()
print(os.getcwd())

basedir=os.getcwd()


import re

from pathlib import Path
from inspect import getsourcefile
import pandas as pd
import numpy as np

#%% Inputs
    


idxfiles=[] #input idxml files

#%% Read modifications

mfile="path/to/unimod.txt"
with open(mfile,"r") as f:
    m=f.read()

ms=m.split("[Term]")[2:]
unimod_names=[i.split("name: ")[1].split("\n")[0] for i in ms]
unimod_masses=[i.split('xref: delta_mono_mass "')[1].split('"')[0] for i in ms]
unimod_df=pd.DataFrame(list(zip(unimod_names,unimod_masses)),columns=["unimod_name","unimod_mass"])
unimod_df["unimod_mass"]=unimod_df["unimod_mass"].astype(float)

#%%
for idxfile in idxfiles:
    with open (idxfile,"r") as f:
        xml_string=f.read()
        
            
    
    #get proteins
    ps=[i.split(">")[0] for i in xml_string.split("<ProteinHit ")][1:]
    parse_targets=["id","accession","score","sequence"]
    ps=[[i.split(p+'=')[1].split('"')[1] for i in ps[1:]] for p in parse_targets]
    protdf=pd.DataFrame(ps,index=parse_targets).T
    protdf=protdf[["id","accession"]].set_index("id")
    
    #get peptides
    ps=[i.split("</PeptideHit>")[0] for i in xml_string.split("<PeptideIdentification score_type=")]
    parse_targets=["sequence","spectrum_reference","charge","MZ","RT","score","protein_refs"]
    ps=[[i.split(p+'=')[1].split('"')[1] for i in ps[1:]] for p in parse_targets]
    pepdf=pd.DataFrame(ps,index=parse_targets).T
    pepdf["First Scan"]=pepdf["spectrum_reference"].apply(lambda x: x.split("=")[-1])
    
    #homogenize Nterm modifications
    pepdf.sequence="."+pepdf.sequence
    pepdf["sequence"]=pepdf["sequence"].str.replace("...",".",regex=False)
    
    #add modifications
    
    #determine bracket type
    if  pepdf.sequence.str.contains("[",regex=False).any(): #not sure if the need for an escape character is platform dependent
        bo="[" #bracket open
        bc="]" #bracket close
        
    if  pepdf.sequence.str.contains("(",regex=False).any(): #not sure if the need for an escape character is platform dependent
        bo="(" #bracket open
        bc=")" #bracket close
    
    
    mods=[]
    for ix_s,s in enumerate(pepdf["sequence"]):
        
        if bo in s: 
        
    
    
            counter=0
            for ix,i in enumerate(s):
                if i==bo:
                   b_open=ix
                if i==bc:
                   b_close=ix
                   counter+=b_close-b_open
                   mods.append([ix_s,s[(b_open+1):b_close],s[b_open-1],ix-1-counter]) #pep_ix mod preceding, ix in s
                   
                   
           
    mods=pd.DataFrame(mods,columns=["peptide_index","modification","preceding_AA","sequence_index"])
    mods["preceding_AA"]=mods["preceding_AA"].replace(".","N-Term(Prot)",regex=False) #if N-term, replace preceding_AA with N-Term(Prot)
    mods["str_ix"]=mods["sequence_index"].astype(str).replace("0","")
    
    
    try: #if listed as masses, try to match to the closest unimod mass
        mods["modification"]=mods["modification"].astype(float) #if float: convert using closeness of mass
        u=mods["modification"].drop_duplicates().reset_index(drop=True)
        um=unimod_df.iloc[np.argmin(abs(unimod_df["unimod_mass"].values-u.values.reshape(-1,1)),axis=1)].reset_index(drop=True)
        um["u"]=u
        mods["modification"]=mods.merge(um,left_on="modification",right_on="u",how="left")["unimod_name"]
        
    except:
        pass
    
    mods["mod_str"]=mods["preceding_AA"]+mods["str_ix"]+"("+mods["modification"]+")"
    pepdf["Modifications"]=mods.sort_values(by=["peptide_index","sequence_index"]).groupby("peptide_index",sort=False).apply(lambda x: "; ".join(x["mod_str"]))
    
    #make underscored peptide sequence
    pepdf["Sequence"]=pepdf.sequence.apply(lambda x: re.sub("[\[\[].*?[\]\]]", "", x).replace(",","")).str.replace(".","",regex=False) #remove ptms in peptides
    pepdf["Sequence"]=pepdf.Sequence.apply(lambda x: re.sub("[\(\[].*?[\)\]]", "", x).replace(",","")).str.replace(".","",regex=False) #remove ptms in peptides
    pepdf["Annotated_Sequence"]=pepdf["Sequence"]
    mods["sequence_index"]=mods["sequence_index"]-1
    mods.loc[mods["sequence_index"]==-1,"sequence_index"]=0
    for index,i in mods.iterrows():
        pepdf["Annotated_Sequence"].iloc[i["peptide_index"]][0:i["sequence_index"]]+pepdf["Annotated_Sequence"].iloc[i["peptide_index"]][i["sequence_index"]].lower()+pepdf["Annotated_Sequence"].iloc[i["peptide_index"]][(i["sequence_index"]+1):]
    
    
    
        
    pepdf["Protein Accessions"]=pepdf["protein_refs"].str.split()
    pepdf=pepdf.explode("Protein Accessions")
    pepdf["Protein Accessions"]=protdf.loc[pepdf["Protein Accessions"]].reset_index(drop=True)
    
    oi=pepdf.groupby(pepdf.columns[:-1].tolist()).size().reset_index(drop=True)
    pepdf=pepdf.groupby(pepdf.columns[:-1].tolist())["Protein Accessions"].apply(lambda x:" ".join(x)).reset_index()
    pepdf["# Proteins"]=oi
    
    
    #fill remainder
    
    pepdf["Master Protein Accessions"]=pepdf["protein_refs"]
    pepdf["PSM Ambiguity"]="Unambiguous"
    pepdf["Confidence"]="High"
    pepdf["Charge"]=pepdf["charge"].astype(int)
    pepdf["m/z [Da]"]=pepdf["MZ"].astype(float)
    pepdf["MH+ [Da]"]=pepdf["m/z [Da]"]*pepdf["Charge"]-pepdf["Charge"]*1.007825319+1.0078250319
    pepdf["RT [min]"]=pepdf["RT"].astype(float)/60
    pepdf["XCorr"]=pepdf["score"]
    pepdf["Spectrum File"]=Path(idxfile).stem+".raw"
    
    
    
    #for this specific diamond pipeine
    
    
    #pepdf["Peptide"]=pepdf["Annotated Sequence"]
    pepdf.to_csv("target_"+Path(idxfile).stem+".csv")
