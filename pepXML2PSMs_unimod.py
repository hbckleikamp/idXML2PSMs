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
import csv
from pathlib import Path
from inspect import getsourcefile
import pandas as pd
import numpy as np

#%% Inputs
    
std_aa_mass = {'G': 57.02146, 'A': 71.03711, 'S': 87.03203, 'P': 97.05276, 'V': 99.06841,
               'T': 101.04768,'C': 103.00919,'L': 113.08406,'I': 113.08406,'J': 113.08406,
               'N': 114.04293,'D': 115.02694,'Q': 128.05858,'K': 128.09496,'E': 129.04259,
               'M': 131.04049,'H': 137.05891,'F': 147.06841,'U': 150.95364,'R': 156.10111,
               'Y': 163.06333,'W': 186.07931,'O': 237.14773}  

idxfiles=[] #input pepxml files

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
    
        
    
    #%%
    
    
    ps=[i.split("</spectrum_query>")[0] for i in xml_string.split("<spectrum_query ")]
    parse_targets=["start_scan","assumed_charge","retention_time_sec","spectrum","precursor_neutral_mass"]
    s=[[i.split(p+'=')[1].split('"')[1] for i in ps[1:]] for p in parse_targets]
    specdf=pd.DataFrame(s,index=parse_targets).T
    
    pepdfs=[]
    for ip,ppp in enumerate(ps[1:]):
  
        pp=[i.split("</search_hit>")[0] for i in ppp.split("<search_hit ")[1:]]
        parse_targets=["peptide","calc_neutral_pep_mass","hit_rank","tot_num_ions","massdiff",'"hyperscore" value']
        s=[[i.split(p+'=')[1].split('"')[1] for i in pp] for p in parse_targets]
        pepdf=pd.DataFrame(s,index=parse_targets).T
        pepdf["Members"]=[i.split('Members=')[1].split('"')[0] for i in pp]
        pepdf["# Proteins"]=[i.split('protein_descr="n=')[1].split()[0] for i in pp]
        
        #nested list comprehension is too hard so make an actual loop
        mods=[]
        for i in pp:
            r=[]
            if "modification_info" in i:
                m=i.split('<mod_aminoacid_mass mass="')[1:]
                p=i.split('position="')[1:]
                for ij,_ in enumerate(m):
                    r.append([m[ij].split('"')[0],p[ij].split('"')[0]])
                    
                    
            mods.append(r)
    
        pepdf["Modification"]=mods
        pepdf["ix"]=ip
        pepdfs.append(pepdf)
    pepdfs=pd.concat(pepdfs).set_index("ix")
    pepdf=specdf.merge(pepdfs,how="right",left_index=True,right_index=True)
    
    #fill remainder
    pepdf["First Scan"]=pepdf["start_scan"]
    pepdf["Protein Accessions"]=pepdf["Members"]
    pepdf["Master Protein Accessions"]=pepdf["Members"]
    pepdf["PSM Ambiguity"]="Unambiguous"
    pepdf["Confidence"]="High"

    pepdf["Charge"]=pepdf["assumed_charge"].astype(int)
    pepdf["m/z [Da]"]=(pepdf["precursor_neutral_mass"].astype(float)+pepdf["Charge"]*1.007825319)/pepdf["Charge"]
    pepdf["MH+ [Da]"]=pepdf["precursor_neutral_mass"].astype(float)+1.0078250319
    
    pepdf["RT [min]"]=pepdf['retention_time_sec'].astype(float)/60
    pepdf["XCorr"]=pepdf['"hyperscore" value']
    pepdf["Spectrum File"]=Path(idxfile).stem+".raw"
    
    
    #%% hardest part: modifications
    mods=[]
    for index,i in pepdf.iterrows():
        if i['Modification']!=[]:
 
            peptide=i["peptide"]
            a=np.array(i['Modification']).astype(float)
            for ai in a:
                mods.append([i.start_scan,
                             i.hit_rank,
                             ai[0]-std_aa_mass.get(peptide[int(ai[1])-1]),
                             peptide[int(ai[1])-1],
                             ai[1]-1]   )
    
    mods=pd.DataFrame(mods,columns=["start_scan","hit_rank","modification","preceding_AA","sequence_index"])

    mods["modification"]=mods["modification"].astype(float) #if float: convert using closeness of mass
    u=mods["modification"].drop_duplicates().reset_index(drop=True)
    um=unimod_df.iloc[np.argmin(abs(unimod_df["unimod_mass"].values-u.values.reshape(-1,1)),axis=1)].reset_index(drop=True)
    um["u"]=u
    mods["modification"]=mods.merge(um,left_on="modification",right_on="u",how="left")["unimod_name"]
    mods["str_ix"]=mods["sequence_index"]+1


    
    mods["mod_str"]=mods["preceding_AA"]+mods["str_ix"].astype(int).astype(str)+"("+mods["modification"]+")"
    mmods=mods.sort_values(by=["start_scan","hit_rank","sequence_index"]).groupby(["start_scan","hit_rank"],sort=False).apply(lambda x: "; ".join(x["mod_str"]))
    mmods.name="Modifications"
    pepdf=pepdf.merge(mmods,on=["start_scan","hit_rank"],how="left")

    #make underscored peptide sequence
    pepdf["Annotated Sequence"]=pepdf["peptide"]

    for index,i in pepdf.iterrows():
        if i['Modification']!=[]:
            #break
            peptide=i["peptide"]
            a=np.array(i['Modification']).astype(float)
            for ai in a:
                pepdf.loc[index,"Annotated Sequence"]= pepdf.loc[index,"Annotated Sequence"][:int(ai[1])-1]+pepdf.loc[index,"Annotated Sequence"][int(ai[1]-1)].lower() +pepdf.loc[index,"Annotated Sequence"][int(ai[1]):]

    
    #%%
        

    
    pepdf=pepdf[[
          "Annotated Sequence" ,
          "Confidence"    ,
          "PSM Ambiguity" ,
          "Modifications" ,
          "First Scan",
          "Spectrum File"  ,
          "Charge" ,
          "m/z [Da]" ,
          "MH+ [Da]" ,
          "RT [min]" ,
          "XCorr" ,
          "Protein Accessions" ,
          "Master Protein Accessions" ,
          "# Proteins" ]].fillna("").drop_duplicates()
    
    pepdf=pepdf[pepdf["Annotated Sequence"]!=""]
    
    
    pepdf.columns=['"'+c+'"' for c in pepdf.columns]
    pepdf='"'+pepdf.astype(str)+'"'
    pepdf.to_csv(Path(idxfile).stem+"_CALISP_PSMS.txt",sep="\t",index=False,quoting=csv.QUOTE_NONE)
