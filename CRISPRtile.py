#!/usr/bin/env python
# coding: utf-8

import glob
import pandas as pd
import numpy as np
import plotly.graph_objects as go
import matplotlib
import kaleido
matplotlib.rcParams['pdf.fonttype'], matplotlib.rcParams['ps.fonttype'] = 42, 42
def findsgrnacolumn(dataframe):
    findsgrnacolumn=list(dataframe.T[0])
    for k in range(len(findsgrnacolumn)):
        if type(findsgrnacolumn[k])==str and len(findsgrnacolumn[k])==20:
            counter=0
            for m in range(len(findsgrnacolumn[k])):
                if findsgrnacolumn[k][m] not in ['A','T','C','G']:
                    counter=1
                if m==19 and counter==0:
                    return(k)
if glob.glob('input/input.csv')!=['input/input.csv']:
    print('Generating input.csv file')
    specificitygroups=input("Enter the column header name of specificity scores in guide specificity file with example in quotation \"SPs1,SPs2,SPs3,SPs5,SPs6,SPs7,SPs8,SPs9\":").split(',')
    funcgroups=input("Enter the column header name of Day6;Day 14;top group;bottom group where groups are seperated by ; and replicates are seperated by , with example in quotation \"D6R1,D6R2,D6R3;D14R1,D14R2,D14R3;AR1,AR2,AR3;UAR1,UAR2,UAR3\". Enter 0 if you have custom scores.:").split(';')
    genecol=input("Enter the column header name of Genes in score file with example in quotation \"Gene\":")
    custom=0
    if funcgroups==['0']:
        funcgroups=input("Enter the column header name of scores where groups are seperated by ; and replicates are seperated by , with example in quotation \"D6R1,D6R2,D6R3;D14R1,D14R2,D14R3;AR1,AR2,AR3;UAR1,UAR2,UAR3\":").split(';')
        custom=1
    csvlist, funccol =glob.glob('input/*.csv'), sum([funcgroups[m].split(',') for m in range(len(funcgroups))],[])
    for k in range(len(csvlist)):
        checkdf=pd.read_csv(csvlist[k])
        if(all(x in checkdf.columns for x in specificitygroups)):
            specificityfile=checkdf.dropna()
        if(all(x in checkdf.columns for x in funccol)):
            functionalfile=checkdf
    specificityfilecol, functionalfilecol=list(specificityfile.columns), list(functionalfile.columns)
    specind=[specificityfilecol[k] for k in range(len(specificityfilecol)) if len(set(specificityfile[specificityfilecol[k]]))==len(specificityfile[specificityfilecol[k]])]
    funcind=[functionalfilecol[m] for m in range(len(functionalfilecol)) if len(set(functionalfile[functionalfilecol[m]]))==len(functionalfile[functionalfilecol[m]])]
    specificitycombine, functionalcombine, maxlen = 0, 0, 0
    for k in range(len(specind)):
        for m in range(len(funcind)):
            intlen=len(set(specificityfile[specind[k]]).intersection(set(functionalfile[funcind[m]])))
            if intlen>maxlen:
                specificitycombine, functionalcombine, maxlen = specind[k], funcind[m], intlen
    controlid=list(set(functionalfile[functionalcombine])-set(specificityfile[specificitycombine]))
    specdf, funcdf=specificityfile.set_index(specificitycombine), functionalfile.set_index(functionalcombine)
    finalfunc=pd.DataFrame(funcdf[funcgroups[0].split(',')].mean(axis=1),columns=[funcgroups[0]])
    for k in range(1,len(funcgroups)):
        finalfunc=finalfunc.join(pd.DataFrame(funcdf[funcgroups[k].split(',')].mean(axis=1),columns=[funcgroups[k]]))
    if custom==0:
        def controlcorrect(a,b):
            matrix=(finalfunc[funcgroups[a]]/finalfunc[funcgroups[b]])
            control=matrix[controlid]
            control.replace([np.inf, -np.inf, 0], np.nan, inplace=True)
            control=control.dropna()
            return(np.log2(matrix/np.median(control)))
        finalfunc=pd.DataFrame(controlcorrect(1,0), columns=['Dropout']).join(pd.DataFrame(controlcorrect(2,3), columns=['Function']))
    graphcol, funcdf =list(finalfunc.columns), funcdf.join(finalfunc)
    def findcategory(df):
        dfcol, colgroup =list(df.columns), []
        for k in range(len(dfcol)):
            dfsubset=df[dfcol[k]]
            if type(dfsubset[1])==str and len(set(dfsubset))<len(dfsubset):
                colgroup.append(dfcol[k]) 
        return(colgroup)
    colgroup=findcategory(funcdf)
    for k in range(len(colgroup)):
        colgrouplist, sgrnaindex, colgrouplistnew=list(funcdf[colgroup[k]]), list(funcdf.index), []
        for m in range(len(colgrouplist)):
            if sgrnaindex[m] in controlid:
                colgrouplistnew.append('Control')
            else:
                colgrouplistnew.append(colgrouplist[m])
        funcdf[colgroup[k]]=colgrouplistnew
    funcgraph, filterlist=funcdf.reset_index(), []
    pregfilter=funcgraph[funccol[6:]].values.tolist()
    for k in range(len(pregfilter)):
        if 0 in pregfilter[k]:
            filterlist+=pregfilter[k]
    percentilecut, filteredlist=(sorted(sum(pregfilter, [])))[int(len((sorted(sum(pregfilter, []))))*.05)], []
    for k in range(len(pregfilter)):
        for m in range(len(pregfilter[k])):
            if percentilecut>pregfilter[k][m]:
                filteredlist.append(funcgraph[functionalcombine][k])
    funcgraph=funcgraph[~funcgraph[functionalcombine].isin(filteredlist)]
    for m in range(len(graphcol)):
        fig = go.Figure()
        for k in range(len(colgroup)):
            genes = set(funcgraph[genecol])
            for gene in genes:
                fig.add_trace(go.Violin(x=funcgraph[genecol][funcgraph[genecol] == gene],y=funcgraph[graphcol[m]][funcgraph[genecol] == gene],
                                        name=gene,box_visible=True,meanline_visible=True, line_color='black',showlegend=False))
                fig.update_layout(xaxis= {'anchor': 'y', 'domain': [0.0, 1.0], 'title': {'text': colgroup[k]}},
                                 yaxis = {'anchor': 'x', 'domain': [0.0, 1.0], 'title': {'text': graphcol[m]}})
            fig.write_image('output/'+colgroup[k]+' '+graphcol[m]+' violin.pdf')
            fig.write_html('output/'+colgroup[k]+' '+graphcol[m]+' violin interactive.html')
    preinputdf=specdf.join(funcdf).reset_index()
    inputdf=preinputdf[specificitygroups+[preinputdf.columns[findsgrnacolumn(preinputdf)]]+list(finalfunc.columns)].set_index(specificitygroups[0])
    inputdf.to_csv('input/input.csv')
    print('input.csv file has been generated.')
print('input.csv file has been found in the input folder. Resuming progress from this point. If you would like to start over, delete all files in the output folder.')
inputdf=pd.read_csv('input/input.csv')
sgrnacolumn=findsgrnacolumn(inputdf)
scoreheader, guideheader=list(inputdf.columns)[sgrnacolumn+1:], list(inputdf.columns)[:sgrnacolumn]
if glob.glob('output/annotation.csv')!=['output/annotation.csv']:
    print('Annotation file has not been generated yet in output folder. Generating annotation file. ')
    import mdtraj as md
    import gzip
    import csv
    def analyzepdb(file):
        pdb = md.load(pdbfiles[file])
        dssps, dssp = list(md.compute_dssp(pdb, simplified=False)[0]), list(md.compute_dssp(pdb)[0])
        ssd=pd.DataFrame(zip(dssp,dssps),index=[(k+1) for k in range(len(dssp))], columns=['Secondary Structure','Secondary Structure Detailed'])
        sasa = md.shrake_rupley(pdb, mode = 'residue')[0]
        sasadf=pd.DataFrame(sasa,index=[(k+1) for k in range(len(sasa))],columns=['Solvent-Accessible Surface Area'])
        d = {'CYS': 'C', 'ASP': 'D', 'SER': 'S', 'GLN': 'Q', 'LYS': 'K',
         'ILE': 'I', 'PRO': 'P', 'THR': 'T', 'PHE': 'F', 'ASN': 'N', 
         'GLY': 'G', 'HIS': 'H', 'LEU': 'L', 'ARG': 'R', 'TRP': 'W', 
         'ALA': 'A', 'VAL':'V', 'GLU': 'E', 'TYR': 'Y', 'MET': 'M'}
        pdbdf=pdb.topology.to_dataframe()[0].replace({"resName": d}).join(pd.DataFrame(pdb.xyz[0], columns=['X_coord','Y_coord','Z_coord']))
        dfxyz=pdbdf[pdbdf['name'].isin(['N','CA','C','O'])]
        atoms=list(set(dfxyz['name']))
        header=['resSeq']+['resName']+[atoms[k]+'_x' for k in range(len(atoms))]+[atoms[k]+'_y' for k in range(len(atoms))]+[atoms[k]+'_z' for k in range(len(atoms))]
        dfxyzvalues, finalxyz, tempxyz=dfxyz.values.tolist(), [], [0 for k in range(len(header))]
        previd=dfxyzvalues[0][3]
        for k in range(len(dfxyzvalues)):
            if dfxyzvalues[k][3]!=previd:
                finalxyz.append(tempxyz)
                tempxyz=[0 for k in range(len(header))]
            previd, tempxyz[0], tempxyz[1]=dfxyzvalues[k][3], dfxyzvalues[k][3], dfxyzvalues[k][4]
            tempxyz[header.index(dfxyzvalues[k][1]+'_x')], tempxyz[header.index(dfxyzvalues[k][1]+'_y')], tempxyz[header.index(dfxyzvalues[k][1]+'_z')]=dfxyzvalues[k][7], dfxyzvalues[k][8], dfxyzvalues[k][9]
            if k==len(dfxyzvalues)-1:
                finalxyz.append(tempxyz)
        pdbanalyzed=pd.DataFrame(finalxyz, columns=header).sort_values(by=['resSeq']).set_index('resSeq').join(ssd).join(sasadf)
        return(pdbanalyzed)
    pdbfiles, aminolist, pdbdfmatrix=glob.glob('input/*.pdb'), [], []
    for k in range(len(pdbfiles)):
        dfanalyzepdb=analyzepdb(k)
        pdbdfmatrix.append(dfanalyzepdb)
        aminolist.append(list(dfanalyzepdb['resName'])+['*'])
    geneid=[[] for k in range(len(aminolist))]
    def annotation():
        flist, prev, prevAA, aminodone=[], 0, 0, []
        for k in range(1,14):
            with gzip.open('annotation/annotation'+str(k)+'.csv.gz', mode='rt') as file:
                reader = csv.reader(file)
                if k==1:
                    for line in reader:
                        if 'ensembl_peptide_id' in line:
                            columnlabel=line
                            break
                for line in reader:
                    if line[2]!=prev:
                        templist, tempamino=[line], [line[7]]
                    if line[2]==prev:
                        templist.append(line)
                        if line[6]!=prevAA:
                            tempamino.append(line[7])
                    prev, prevAA=line[2], line[6]
                    if tempamino in aminolist:
                        geneid[aminolist.index(tempamino)].append(line[4])
                        flist+=templist
                        print(line[4]+' annotation is complete')
                        if tempamino not in aminodone:
                            aminodone.append(tempamino)
                            if len(aminodone)==len(aminolist):
                                return(pd.DataFrame(flist,columns=columnlabel))
                        templist, tempamino=[], []
            print(str(7.6*k)+'% complete')
    annotation1, masterlistmatrix = annotation(), []
    for k in range(len(geneid)):
        pdbannotation=annotation1[annotation1['transcript_name'].isin([geneid[k][0]])].copy().reset_index()
        pdbannotation['index']=pdbannotation['position'].astype(int)
        masterlistmatrix.append(pdbannotation.set_index('index').join(pdbdfmatrix[k]))
    premlist=pd.concat([masterlistmatrix[k] for k in range(len(masterlistmatrix))]).reset_index()
    premlist['index']=premlist['guide']
    mlist1=premlist.set_index('index')
    mlist2=mlist1.join(inputdf.set_index(list(inputdf.columns)[findsgrnacolumn(inputdf)]))
    listsecstrucdetailed, secstrucnew=list(mlist2['Secondary Structure Detailed']), []
    for k in range(len(listsecstrucdetailed)):
        if listsecstrucdetailed[k]==' ':
            secstrucnew.append('N')
        else:
            secstrucnew.append(listsecstrucdetailed[k])
    mlist2['Secondary Structure Detailed'], mlist2['SecStruct']=secstrucnew, mlist2['Secondary Structure'].replace(['E'], 'B')
    mlist2=mlist2.sort_values(by=['transcript_name','gene_fraction']).reset_index().drop(columns=['index','resName', 'Secondary Structure']).dropna(subset=['Solvent-Accessible Surface Area'])
    MLprep=mlist2.filter(items=['transcript_name','AA', 'offtarget_score', 'doench_score', 'oof_score', 'position', 'provean_score', 'disorder_score']
                     +((list(mlist2.columns)[list(mlist2.columns).index('control')+1:])))
    onehotindex, transcript, AA, secstruc=list(MLprep.index), MLprep['transcript_name'], MLprep['AA'], MLprep['Secondary Structure Detailed']
    transcriptset, AAset, secstrucset=list(set(transcript)), list(set(AA)), list(set(secstruc))
    onehotset=transcriptset+['AA-'+str(AAset[k]) for k in range(len(AAset))]+['SS-'+str(secstrucset[k]) for k in range(len(secstrucset))]
    zeromatrix=[0 for k in range(len(onehotset))]
    tempzero, onehotlist=zeromatrix.copy(), []
    for k in range(len(AA)):
        if list(transcript)[k] in onehotset:
            tempzero[onehotset.index(list(transcript)[k])]=1
        if 'AA-'+list(AA)[k] in onehotset:
            tempzero[onehotset.index('AA-'+list(AA)[k])]=1
        if 'SS-'+list(secstruc)[k] in onehotset:
            tempzero[onehotset.index('SS-'+list(secstruc)[k])]=1
        onehotlist.append([onehotindex[k]]+tempzero)
        tempzero=zeromatrix.copy()
    onehotdf=pd.DataFrame(onehotlist, columns=['numberindex']+onehotset).set_index('numberindex')
    MLfinal=(onehotdf.join(MLprep)).drop(columns=["transcript_name","AA","Secondary Structure Detailed"])
    MLfinal.to_csv('output/train.csv')
    mlist2.to_csv('output/annotation.csv')
    print('Annotation file has been saved. You can resume progress from this point by running the program again.')
print('Generated annotation file has been found in the output folder. Resuming progress from this point. If you would like to start over, delete all files in the output folder.')
from autogluon.tabular import TabularDataset, TabularPredictor
if glob.glob('output/predict.csv')!=['output/predict.csv']:
    import torch
    MLfinal=pd.read_csv('output/train.csv').drop(columns=['numberindex'])
    MLfinal.replace([np.inf, -np.inf, 'NA'], np.nan, inplace=True)
    for k in range(len(scoreheader)):
        if glob.glob('output/Models/'+scoreheader[k])!=['output/Models/'+scoreheader[k]]:
            MLtrain=MLfinal.dropna(subset=[scoreheader[k]]).drop(columns=(scoreheader[:k]+scoreheader[k+1:]))
            predictor=TabularPredictor(label=scoreheader[k], path="output/Models/"+scoreheader[k]).fit(MLtrain, presets="best_quality", num_gpus=torch.cuda.device_count())
    correctlist=guideheader+['offtarget_score', 'doench_score', 'oof_score']
    for k in range(len(correctlist)):
        if type(MLfinal[correctlist[k]].mode()[0])==str:
            MLfinal[correctlist[k]]=[MLfinal[correctlist[k]].mode()[0] for m in range(len(MLfinal))]
        if type(MLfinal[correctlist[k]].mode()[0])!=str:
            MLfinal[correctlist[k]]=[MLfinal[correctlist[k]].median() for m in range(len(MLfinal))]
    mlist2=pd.read_csv('output/annotation.csv').drop(columns=['Unnamed: 0'])
    for k in range(len(scoreheader)):
        predictor = TabularPredictor.load("output/Models/"+scoreheader[k])
        predictdf=pd.DataFrame(predictor.predict(MLfinal))
        predictdf.rename(columns = {scoreheader[k]:scoreheader[k]+' predict'}, inplace = True)
        mlist2=mlist2.join(predictdf)
    mlist2.to_csv('output/predict.csv')
import plotly.express as px
from plotly import subplots
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA
mlist2=pd.read_csv('output/predict.csv').drop(columns=['Unnamed: 0'])
mlist2.replace(['#NAME?'], np.nan, inplace=True)
mlist2['Interpro_Description'].replace([np.nan], 'NA', inplace=True)
mlist2['Secondary Structure Detailed'].replace(['N'], 'Loops and irregular elements', inplace=True)
mlist2['Secondary Structure Detailed'].replace(['H'], 'Alpha helix', inplace=True)
mlist2['Secondary Structure Detailed'].replace(['B'], 'Residue in isolated beta-bridge', inplace=True)
mlist2['Secondary Structure Detailed'].replace(['E'], 'Extended strand, participates in beta ladder', inplace=True)
mlist2['Secondary Structure Detailed'].replace(['G'], '3-helix (3/10 helix)', inplace=True)
mlist2['Secondary Structure Detailed'].replace(['I'], '5 helix (pi helix)', inplace=True)
mlist2['Secondary Structure Detailed'].replace(['T'], 'hydrogen bonded turn', inplace=True)
mlist2['Secondary Structure Detailed'].replace(['S'], 'bend', inplace=True)
for k in range(len(scoreheader)):
    mlist2[scoreheader[k]]=mlist2[scoreheader[k]].astype(float)
transcriptset, predictcol=list(set(mlist2['transcript_name'])), [scoreheader[m]+' predict' for m in range(len(scoreheader))]
parallelcollist=['position','provean_score','disorder_score','Solvent-Accessible Surface Area']+predictcol
for k in range(len(transcriptset)):
    plotdf=mlist2[mlist2['transcript_name'].isin([transcriptset[k]])]
    fig = go.Figure(data=go.Parcoords(line = dict(color = plotdf['position'],colorscale = px.colors.diverging.Tealrose),
        dimensions = [dict(label = parallelcollist[k], values = plotdf[parallelcollist[k]]) for k in range(len(parallelcollist))]),
                   layout = go.Layout(autosize=False,width=800,height=500, title=transcriptset[k]))
    fig.write_image('output/'+transcriptset[k]+'_parallelplot.pdf')
    fig.write_html('output/'+transcriptset[k]+' parallelplot interactive.html')
    def pca(columns, whichpca):
        X = plotdf[columns]
        pca = PCA(n_components=2)
        components = pca.fit_transform(X)
        loadings = pca.components_.T * np.sqrt(pca.explained_variance_)
        fig = px.scatter(components, x=0, y=1, color=plotdf[predictcol[m]], size=plotdf['Solvent-Accessible Surface Area'])
        for i, feature in enumerate(columns):
            fig.add_annotation(ax=0, ay=0,axref="x", ayref="y",x=loadings[i, 0],y=loadings[i, 1],showarrow=True,arrowsize=2,arrowhead=2,xanchor="right",yanchor="top")
            fig.add_annotation(x=loadings[i, 0],y=loadings[i, 1],ax=0, ay=0,xanchor="center", yanchor="bottom",text=feature,yshift=5)
        fig.update_layout(title={'text': transcriptset[k],'y':0.95,'x':0.5,'xanchor': 'center','yanchor': 'top'},
                          xaxis= {'anchor': 'y', 'domain': [0.0, 1.0], 'title': {'text': whichpca+'PC1 ('+str(round(pca.explained_variance_ratio_[0], 3)*100)+'% Variance)'}},
                         yaxis = {'anchor': 'x', 'domain': [0.0, 1.0], 'title': {'text': whichpca+'PC2 ('+str(round(pca.explained_variance_ratio_[1], 3)*100)+'% Variance)'}},
                         coloraxis = {'colorbar': {'title': {'text': scoreheader[m]}}})
        fig.write_image('output/'+transcriptset[k]+' '+scoreheader[m]+' '+whichpca+'PCA.pdf')
        fig.write_html('output/'+transcriptset[k]+' '+scoreheader[m]+' '+whichpca+'PCA interactive.html')
    for m in range(len(predictcol)):
        pca(['CA_x','CA_y','CA_z'], '3D ')
        pca(['provean_score','disorder_score','Solvent-Accessible Surface Area'], 'PDS ')
        violingroup=['AA','Interpro_Description','Secondary Structure Detailed']
        for n in range(len(violingroup)):
            fig = go.Figure()
            groups = set(plotdf[violingroup[n]])
            if n==0:
                for group in groups:
                    fig.add_trace(go.Violin(x=plotdf[violingroup[n]][plotdf[violingroup[n]] == group],y=plotdf[predictcol[m]][plotdf[violingroup[n]] == group],
                                    name=group,box_visible=True,meanline_visible=True, line_color='black',showlegend=False))
            else:
                for group in groups:
                    fig.add_trace(go.Violin(x=plotdf[violingroup[n]][plotdf[violingroup[n]] == group],y=plotdf[predictcol[m]][plotdf[violingroup[n]] == group],
                                    name=group,box_visible=True,meanline_visible=True))
            fig.update_layout(title={'text': transcriptset[k],'y':0.95,'x':0.5,'xanchor': 'center','yanchor': 'top'},
                                  xaxis= {'anchor': 'y', 'domain': [0.0, 1.0], 'title': {'text': violingroup[n]}},
                                 yaxis = {'anchor': 'x', 'domain': [0.0, 1.0], 'title': {'text': scoreheader[m]}})
            fig.write_image('output/'+transcriptset[k]+' '+scoreheader[m]+' '+violingroup[n]+' violin.pdf')
            fig.write_html('output/'+transcriptset[k]+' '+scoreheader[m]+' '+violingroup[n]+' violin interactive.html')
        fig2 = px.scatter(plotdf, x="position", y=predictcol[m], color="Interpro_Description", marginal_y="box", size="Solvent-Accessible Surface Area")
        fig2.update_yaxes(title_text=scoreheader[m], col=1)
        fig2.update_xaxes(title_text='Position', col=1)
        fig2.update_layout(title_text=transcriptset[k])
        fig2.write_image('output/'+transcriptset[k]+' '+scoreheader[m]+' SASAsize.pdf')
        fig2.write_html('output/'+transcriptset[k]+' '+scoreheader[m]+' SASAsize interactive.html')
    position, exon, lineargroups=list(plotdf['position']), list(plotdf['Exon']), scoreheader+['provean_score','disorder_score']
    fig=subplots.make_subplots(rows=(len(lineargroups)),cols=1, vertical_spacing = 0.01)
    for t in range(len(lineargroups)):
        if lineargroups[t] in scoreheader:
            fig.add_trace(go.Scatter(x=position, y=plotdf[lineargroups[t]], name = lineargroups[t]+' raw', mode='markers', legendgroup = str(t+1)), row=t+1, col=1)
            predict=list(plotdf[lineargroups[t]+' predict'])
            predictor = TabularPredictor.load("output/Models/"+lineargroups[t])
            error=abs(predictor.leaderboard()['score_val'][0])
            errorupper, errorlower =[predict[k]+error for k in range(len(predict))], [predict[k]-error for k in range(len(predict))]
            fig.add_trace(go.Scatter(x=position+position[::-1], y=errorupper+errorlower[::-1], name = 'k-fold cross validation RMSE', fill='toself', line=dict(width=0), legendgroup = str(t+1)), row=t+1, col=1)
            fig.add_trace(go.Scatter(x=position, y=predict, name = lineargroups[t]+' guide efficiency corrected', legendgroup = str(t+1)), row=t+1, col=1)
        if lineargroups[t] == 'provean_score':
            setdomain=list(set(plotdf["Interpro_Description"]))
            for m in range(len(setdomain)):
                domainsubset=plotdf[plotdf['Interpro_Description'].isin([setdomain[m]])]
                fig.add_trace(go.Scatter(x=domainsubset["position"], y=domainsubset[lineargroups[t]], name = setdomain[m], mode='markers', legendgroup = str(t+1)), row=t+1, col=1)
        if lineargroups[t] == 'disorder_score':
            setss=list(set(plotdf["Secondary Structure Detailed"]))
            for m in range(len(setss)):
                ssubset=plotdf[plotdf["Secondary Structure Detailed"].isin([setss[m]])]
                fig.add_trace(go.Scatter(x=ssubset["position"], y=ssubset[lineargroups[t]], name = setss[m], mode='markers', legendgroup = str(t+1)), row=t+1, col=1)
        fig.update_yaxes(title_text=lineargroups[t], row=t+1, col=1)
    fig.update_layout(height=275*len(lineargroups), width=1100, legend_tracegroupgap = 180, title_text=transcriptset[k], colorway=px.colors.qualitative.Dark24)
    fig.update_xaxes(range=[0,max(list(position))+1],showgrid=False,showticklabels=False)
    fig.update_xaxes(title_text='Position', showticklabels=True, row=len(lineargroups), col=1)
    eold, pold=exon[0], position[0]
    for e in range(len(exon)):
        if exon[e]!=eold:
            fig.add_vline(x=(position[e]+pold)/2, line_dash="dot", row='all', col=1, line_color="#000000", line_width=2)
            pold, eold=position[e], exon[e]
    fig.write_image('output/'+transcriptset[k]+' linear tracks.pdf')
    fig.write_html('output/'+transcriptset[k]+' linear tracks interactive.html')

