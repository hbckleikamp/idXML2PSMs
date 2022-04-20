# -*- coding: utf-8 -*-
"""
Created on Tue Apr 19 17:07:50 2022

@author: ZR48SA
"""


import numpy as np
import pandas as pd
import re


files=["*.idXML"]

for file in files:
    with open (file,"r") as f:
        xml_string=f.read()
    
    
    
    #%% get modifications
    ms=[i.split(">")[0] for i in xml_string.split('_modifications" value="[')]
    ms=[i.split("]")[0] for i in ms[1:]]
    ms=[i.replace(" (",")").split(")")[:-1] for i in ms if i !=""] #remove empy entries, replace " (" with ")", split at ")"
    ms={i[0]:i[1] for i in ms} #make into name list and s
    
    #%% get proteins
    ps=[i.split(">")[0] for i in xml_string.split("<ProteinHit id=")]
    ps=[i.replace('=', ' ').split()[::2] for i in ps[1:]]
    pdf=pd.DataFrame(ps,columns=["protein","accession","score","sequence"])
    pdf=pdf.applymap(lambda x: x[1:-1]).set_index("protein")["accession"]
    
    #%%
    ps=[i.split("</PeptideHit>")[0] for i in xml_string.split("<PeptideIdentification score_type=")]
    possible_modifications={"Carbamidomethyl":"c","Oxidation":"o"} #this would need to be manually altered for each file if others are used
    
    parse_targets=["sequence","spectrum_reference","charge","MZ","RT","score","protein_refs"]
    ps=[[i.split(p+'=')[1].split('"')[1] for i in ps[1:]] for p in parse_targets]
    pepdf=pd.DataFrame(ps)
    pepdf=pepdf.T
    pepdf.columns=parse_targets
    
    pepdf["First Scan"]=pepdf["spectrum_reference"].apply(lambda x: x.split("=")[-1])
    
    
    #%% sequence and modifications
    
    #"standard" modifications
    pepdf["Annotated Sequence"]=pepdf.sequence
    for k,v in ms.items():
        l=v.lower()
        pepdf["Annotated Sequence"]=pepdf["Annotated Sequence"].str.replace('('+k+')',l,regex=False)
    
    #custom n and cterm modifications
    pepdf["Sequence"]=pepdf.sequence.apply(lambda x: re.sub("[\(\[].*?[\)\]]", "", x).replace(",","")) #remove ptms in peptides
    
    qc=pepdf["Annotated Sequence"].str.contains("neutr_cterm")
    qn=pepdf["Annotated Sequence"].str.contains("neutr_nterm")
    pepdf["Annotated Sequence"]=pepdf["Annotated Sequence"].apply(lambda x: re.sub("[\(\[].*?[\)\]]", "", x).replace(",","")) #remove ptms in peptides
    pepdf.loc[qn,"Annotated Sequence"]=pepdf.loc[qn,"Annotated Sequence"].apply(lambda x: x[0].lower()+x[1:])
    pepdf.loc[qc,"Annotated Sequence"]=pepdf.loc[qn,"Annotated Sequence"].apply(lambda x: x[:-1]+x[-1].lower())
    
    
    jas=[]
    for k in ms.keys():
        c=pepdf[pepdf.sequence.str.contains(k)].sequence
        ja=c.apply(lambda x: np.cumsum([len(re.sub("[\(\[].*?[\)\]]", "", i)) for i in x.split("("+k+")")])).reset_index()
        ja["prefix"]=ms.get(k)
        ja["suffix"]="("+k+")"
        ja=ja.explode("sequence")
        ja["order"]=ja.sequence
        jas.append(ja)
        
    #find c-term n-term
    c=pepdf[pepdf.sequence.str.contains("_neutr_nterm")].sequence
    ja=c.apply(lambda x: x.split("_neutr_nterm")[0].str.split("(")[1]+"_neutr_nterm").reset_index()
    ja["prefix"]="N-Term("
    ja["suffix"]="_neutr_nterm)"
    ja["order"]=0
    jas.append(ja)
    
    c=pepdf[pepdf.sequence.str.contains("_neutr_cterm")].sequence
    ja=c.apply(lambda x: x.split("_neutr_cterm")[0].str.split("(")[1]+"_neutr_cterm").reset_index()
    ja["prefix"]="C-Term("
    ja["suffix"]="_neutr_cterm)"
    ja["order"]=100000 #arbitraty big nr
    jas.append(ja)
    
    jas=pd.concat(jas)
    
    jas["sequence"]=jas.sequence.astype(str)
    pepdf["Modifications"]=jas.sort_values(by=["index","order"]).groupby("index",sort=False).apply(lambda x: "; ".join(x["prefix"]+x["sequence"]+x["suffix"].tolist()))
    #%% add proteins
    
    pepdf["Protein Accessions"]=pepdf["protein_refs"].str.split()
    pepdf=pepdf.explode("Protein Accessions")
    pepdf["Protein Accessions"]=pdf.loc[pepdf["Protein Accessions"]].reset_index(drop=True)
    
    oi=pepdf.groupby(pepdf.columns[:-1].tolist()).size().reset_index(drop=True)
    pepdf=pepdf.groupby(pepdf.columns[:-1].tolist())["Protein Accessions"].apply(lambda x:" ".join(x)).reset_index()
    pepdf["# Proteins"]=oi
    
    
    #%% fill remainder
    
    pepdf["Master Protein Accessions"]=pepdf["protein_refs"]
    pepdf["PSM Ambiguity"]="Unambiguous"
    pepdf["Confidence"]="High"
    pepdf["Charge"]=pepdf["charge"].astype(int)
    pepdf["m/z [Da]"]=pepdf["MZ"].astype(float)
    pepdf["MH+ [Da]"]=pepdf["m/z [Da]"]*pepdf["Charge"]-pepdf["Charge"]*1.007825319+1.0078250319
    pepdf["RT [min]"]=pepdf["RT"]/60
    pepdf["XCorr"]=pepdf["score"]
    pepdf["Spectrum File"]=file
    
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
          "# Proteins" ]]
    
    pepdf.to_csv(file.replace(".idXML","_PSMS.txt"),sep="\t")
