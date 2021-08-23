import pandas as pd
import numpy as np

data = pd.read_table('/Users/catherinehogg/Documents/Semester3/Project/Scripts/isoformproject/local/localdata/PvR_isoformfpkm_all.txt',sep='\t',header=0,index_col=0)
patients= pd.read_table('/Users/catherinehogg/Documents/Semester3/Project/Scripts/isoformproject/local/localdata/test_patients.txt',sep='\t',header=None,index_col=0).index
del data['GeneName']
del data['GeneType']


for c in data.columns:
    id=c[:c.find('_')]
    if id not in patients:
        del data[c]

d={}
d['all']=data

for di in d:
    all_values=[]
    for col in d[di]:
        all_values.extend(list(d[di][col]))
    all_values=pd.DataFrame(all_values)
    lq=all_values[all_values >0].quantile(0.25)
    print(str(di)+"_"+str(lq[0]))
    ps=[]
    for i in range(0,len(d[di])):
        abovep=d[di].iloc[i][d[di].iloc[i] >= lq[0]][d[di].iloc[i][d[di].iloc[i] >= lq[0]].index.str.contains('_P')].count()
        totalp=d[di].iloc[i][d[di].iloc[i].index.str.contains('_P')].count()
        proportionp=abovep/totalp
        abover=d[di].iloc[i][d[di].iloc[i] >= lq[0]][d[di].iloc[i][d[di].iloc[i] >= lq[0]].index.str.contains('_R')].count()
        totalr=d[di].iloc[i][d[di].iloc[i].index.str.contains('_R')].count()
        proportionr=abover/totalr

        #Mark rows based on primary proportion or recurrent proportion being >=0.5
        if proportionp >=0.2 or proportionr >=0.2:
            ps.append('True_'+str(proportionp)+'_'+str(abovep)+'_'+str(totalp)+'_'+str(proportionr)+'_'+str(abover)+'_'+str(totalr))
        else:
            ps.append('False_'+str(proportionp)+'_'+str(abovep)+'_'+str(totalp)+'_'+str(proportionr)+'_'+str(abover)+'_'+str(totalr))

    #Filter data based on marked rows
    data['ProportionAbove0_'+di]=ps
    d[di]=d[di][data['ProportionAbove0_'+di].str.contains('True')]

    #Create outputs
    log2fc=pd.DataFrame()
    log2fc['genes']=d[di].index
    primout=pd.DataFrame()
    primout['genes']=d[di].index
    recuout=pd.DataFrame()
    recuout['genes']=d[di].index

    for patient in patients:
        prim=''
        recu=''
        if patient +'_Primary_FPKM' in d[di].columns:
            prim=d[di][patient+'_Primary_FPKM']
            recu=d[di][patient+'_Recurrent_FPKM']
        elif patient +'_P_FPKM' in d[di].columns:
            prim=d[di][patient+'_P_FPKM']
            recu=d[di][patient+'_R_FPKM']
        else:
            print(patient)
            continue
        log2fc[patient]=list(np.log2((recu+0.01)/(prim+0.01)))
        primout[patient]=list(prim)
        recuout[patient]=list(recu)
    log2fc.to_csv('/Users/catherinehogg/Documents/Semester3/Project/Results/localresults/test_log2fc_'+di+'.txt',sep='\t',header=True,index=False)
    primout.to_csv('/Users/catherinehogg/Documents/Semester3/Project/Results/localresults/test_primary_'+di+'.txt',sep='\t',header=True,index=False)
    recuout.to_csv('/Users/catherinehogg/Documents/Semester3/Project/Results/localresults/test_recurrent_'+di+'.txt',sep='\t',header=True,index=False)
