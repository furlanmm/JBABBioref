import matplotlib.pyplot as plt
import matplotlib as mpl
import seaborn as sns
import numpy as np
import pandas as pd
import os
import shutil
import sys
import math

#plt.rcParams['figure.constrained_layout.use'] = True

def listar_subpastas(pasta):
    caminhos = [os.path.join(pasta, nome) for nome in os.listdir(pasta)]
    subpastas = [subpasta for subpasta in caminhos if os.path.isdir(subpasta)]
    return subpastas

    #divValue = 1000
    #columnName = "Produced"
def createHeatChart(meansResDF,divValue,columnName, nameFile,nameFolder):
    lista_dict = {0:'VL',1:'L',2:'N',3:'H',4:'VH'}
    lista_dictValue = {0:'L/H',1:'L/N',2:'N/N',3:'H/N',4:'H/L'}
    meansResDF=meansResDF.rename(columns={"CakeInv":"Cake Inventory","Biogas":"Biogas Production","TurboGen":"Turbo Generator Capacity","byProd":"Vinasse and Cake Production","Biomet":"Biometane Production", "Value":"Biometane/Energy Value"})
    #meansResDF.to_excel("saida.xlsx")
    meansResDF[columnName]=meansResDF[columnName].div(divValue)
    Table1 = meansResDF.pivot_table(index="Biogas Production", columns="Cake Inventory", values=columnName,aggfunc="mean")[::-1]
    Table1.rename(index=lista_dict, columns=lista_dict, inplace=True)
    Table2 = meansResDF.pivot_table(index="Cake Inventory", columns="Biometane Production", values=columnName,aggfunc="mean")[::-1]
    Table2.rename(index=lista_dict, columns=lista_dict, inplace=True)
    Table3 = meansResDF.pivot_table(index="Cake Inventory", columns="Vinasse and Cake Production", values=columnName,aggfunc="mean")[::-1]
    Table3.rename(index=lista_dict, columns=lista_dict, inplace=True)
    Table4 = meansResDF.pivot_table(index="Cake Inventory", columns="Turbo Generator Capacity", values=columnName,aggfunc="mean")[::-1]
    Table4.rename(index=lista_dict, columns=lista_dict, inplace=True)
    Table5 = meansResDF.pivot_table(index="Biogas Production", columns="Turbo Generator Capacity", values=columnName,aggfunc="mean")[::-1]
    Table5.rename(index=lista_dict, columns=lista_dict, inplace=True)
    Table6 = meansResDF.pivot_table(index="Biometane Production", columns="Biogas Production", values=columnName,aggfunc="mean")[::-1]
    Table6.rename(index=lista_dict, columns=lista_dict, inplace=True)
    Table7 = meansResDF.pivot_table(index="Biometane Production", columns="Vinasse and Cake Production", values=columnName,aggfunc="mean")[::-1]
    Table7.rename(index=lista_dict, columns=lista_dict, inplace=True)
    Table8 = meansResDF.pivot_table(index="Biogas Production", columns="Vinasse and Cake Production", values=columnName,aggfunc="mean")[::-1]
    Table8.rename(index=lista_dict, columns=lista_dict, inplace=True)
    Table9 = meansResDF.pivot_table(index="Biometane Production", columns="Turbo Generator Capacity", values=columnName,aggfunc="mean")[::-1]
    Table9.rename(index=lista_dict, columns=lista_dict, inplace=True)
    Table10 = meansResDF.pivot_table(index="Vinasse and Cake Production", columns="Turbo Generator Capacity", values=columnName,aggfunc="mean")[::-1]
    Table10.rename(index=lista_dict, columns=lista_dict, inplace=True)
    Table11 = meansResDF.pivot_table(index="Biometane/Energy Value", columns="Cake Inventory", values=columnName,aggfunc="mean")#[::-1]
    Table11.rename(index=lista_dictValue, columns=lista_dict, inplace=True)
    Table12 = meansResDF.pivot_table(index="Biometane/Energy Value", columns="Biogas Production", values=columnName,aggfunc="mean")#[::-1]
    Table12.rename(index=lista_dictValue, columns=lista_dict, inplace=True)
    Table13 = meansResDF.pivot_table(index="Biometane/Energy Value", columns="Biometane Production", values=columnName,aggfunc="mean")#[::-1]
    Table13.rename(index=lista_dictValue, columns=lista_dict, inplace=True)
    Table14 = meansResDF.pivot_table(index="Biometane/Energy Value", columns="Vinasse and Cake Production", values=columnName,aggfunc="mean")#[::-1]
    Table14.rename(index=lista_dictValue, columns=lista_dict, inplace=True)
    Table15 = meansResDF.pivot_table(index="Biometane/Energy Value", columns="Turbo Generator Capacity", values=columnName,aggfunc="mean")#[::-1]
    Table15.rename(index=lista_dictValue, columns=lista_dict, inplace=True)


    fig, axes = plt.subplots(5, 3, figsize=(15, 25), gridspec_kw={'width_ratios': [5, 5, 3]})
    minValue = max(0,min(Table1.min().min(),Table2.min().min(),Table3.min().min(),Table4.min().min(),Table5.min().min(),Table6.min().min(),Table7.min().min(),Table8.min().min(),Table9.min().min(),Table10.min().min()))
    maxValue = max(Table1.max().max(),Table2.max().max(),Table3.max().max(),Table4.max().max(),Table5.max().max(),Table6.max().max(),Table7.max().max(),Table8.max().max(),Table9.max().max(),Table10.max().max())
    minValue = math.floor(minValue*10)/10
    maxValue = math.ceil(maxValue*10)/10
    spaces = (maxValue - minValue) / 10
    cbar_ax = fig.add_axes([0.95, 0.2, 0.01, 0.6]) 
    norm = mpl.colors.Normalize(vmin=minValue, vmax=maxValue)
    cbarTodas = plt.colorbar(mpl.cm.ScalarMappable(norm=norm, cmap="coolwarm"), cax=cbar_ax, orientation="vertical", location="right",format="%.1f")
    if spaces > 0.1:
        cbarTodas.set_ticks([minValue+i*spaces for i in range(11)])
    else:
        spaces = (maxValue - minValue) / 2
        cbarTodas.set_ticks([minValue+i*spaces for i in range(3)])
    #print(minValue,maxValue,spaces)
    #print([minValue+i*spaces for i in range(11)])
    #print(np.arange(minValue, maxValue + spaces + 1e-4, spaces))
    sns.heatmap(Table1, linewidths=.5, annot=True, fmt=".1f", ax=axes[1,0] , cmap="coolwarm",cbar=False, vmin=minValue, vmax=maxValue)
    sns.heatmap(Table2, linewidths=.5, annot=True, fmt=".1f", ax=axes[4,0], cmap="coolwarm",cbar=False, vmin=minValue, vmax=maxValue)
    sns.heatmap(Table3, linewidths=.5, annot=True, fmt=".1f", ax=axes[4,1], cmap="coolwarm",cbar=False, vmin=minValue, vmax=maxValue)
    sns.heatmap(Table4, linewidths=.5, annot=True, fmt=".1f", ax=axes[4,2], cmap="coolwarm",cbar=False, vmin=minValue, vmax=maxValue)
    sns.heatmap(Table5, linewidths=.5, annot=True, fmt=".1f", ax=axes[1,2], cmap="coolwarm",cbar=False, vmin=minValue, vmax=maxValue)
    sns.heatmap(Table6, linewidths=.5, annot=True, fmt=".1f", ax=axes[2,1], cmap="coolwarm",cbar=False, vmin=minValue, vmax=maxValue)
    sns.heatmap(Table7, linewidths=.5, annot=True, fmt=".1f", ax=axes[2,0], cmap="coolwarm",cbar=False, vmin=minValue, vmax=maxValue)
    sns.heatmap(Table8, linewidths=.5, annot=True, fmt=".1f", ax=axes[1,1], cmap="coolwarm",cbar=False, vmin=minValue, vmax=maxValue)
    sns.heatmap(Table9, linewidths=.5, annot=True, fmt=".1f", ax=axes[2,2], cmap="coolwarm",cbar=False, vmin=minValue, vmax=maxValue)
    sns.heatmap(Table10, linewidths=.5, annot=True, fmt=".1f", ax=axes[3,2], cmap="coolwarm",cbar=False, vmin=minValue, vmax=maxValue)

    sns.heatmap(Table11, linewidths=.5, annot=True, fmt=".1f", ax=axes[3,0], cmap="coolwarm",cbar=False, vmin=minValue, vmax=maxValue)
    sns.heatmap(Table12, linewidths=.5, annot=True, fmt=".1f", ax=axes[0,1], cmap="coolwarm",cbar=False, vmin=minValue, vmax=maxValue)
    sns.heatmap(Table13, linewidths=.5, annot=True, fmt=".1f", ax=axes[0,0], cmap="coolwarm",cbar=False, vmin=minValue, vmax=maxValue)
    sns.heatmap(Table14, linewidths=.5, annot=True, fmt=".1f", ax=axes[3,1], cmap="coolwarm",cbar=False, vmin=minValue, vmax=maxValue)
    sns.heatmap(Table15, linewidths=.5, annot=True, fmt=".1f", ax=axes[0,2], cmap="coolwarm",cbar=False, vmin=minValue, vmax=maxValue)
    plt.tight_layout(rect=[0, 0, 0.95, 1],pad=2)
    #plt.show()
    plt.savefig(nameFolder+nameFile)
    plt.close()



if __name__ == '__main__':

    folderRes = sys.argv[1]
    numPeriods = int(sys.argv[2])
    periodsN2 = ['Week'+str(t) for t in range(1,numPeriods+1)]
    listaPastas = listar_subpastas(folderRes)
    meansRes = np.zeros((len(listaPastas), 26))
    indicedf = []
    for pasta in listaPastas:
        fileResDF = pd.read_excel(pasta+"/Vars.xlsx",sheet_name="Sheet1",index_col=0,header=0)
        indicePasta = listaPastas.index(pasta)
        meansRes[indicePasta][0] = pasta[-6]
        meansRes[indicePasta][1] = pasta[-5]
        meansRes[indicePasta][2] = pasta[-4]
        meansRes[indicePasta][3] = pasta[-3]
        meansRes[indicePasta][4] = pasta[-2]
        meansRes[indicePasta][5] = pasta[-1]
        #indicedf.append(list(pasta[-5:]))
        #indicedf.append(pasta[-5:])
        meansRes[indicePasta][6] = sum(fileResDF[w]["EG"] for w in periodsN2)
        meansRes[indicePasta][7] = sum(fileResDF[w]["EB"] for w in periodsN2)
        meansRes[indicePasta][8] = sum(fileResDF[w]["EM"]+fileResDF[w]["EP"] for w in periodsN2)
        meansRes[indicePasta][9] = sum(fileResDF[w]["EE"] for w in periodsN2)
        meansRes[indicePasta][10] = fileResDF["Week 0"]["FO"]
        meansRes[indicePasta][11] = fileResDF["Week 0"]["Time"]
        meansRes[indicePasta][12] = sum(fileResDF[w]["PU_Especial"]+fileResDF[w]["PU_Extra"]+fileResDF[w]["PU_VHP"] for w in periodsN2)
        meansRes[indicePasta][13] = sum(fileResDF[w]["PU_AEHC"]+fileResDF[w]["PU_AEAC"] for w in periodsN2)
        meansRes[indicePasta][14] = sum(fileResDF[w]["M"] for w in periodsN2)
        meansRes[indicePasta][15] = sum(fileResDF[w]["PUvin"] for w in periodsN2)
        meansRes[indicePasta][16] = sum(fileResDF[w]["Fvin"] for w in periodsN2)
        meansRes[indicePasta][17] = sum(fileResDF[w]["PUcake"] for w in periodsN2)
        meansRes[indicePasta][18] = sum(fileResDF[w]["Icake"] for w in ["Week 0"]+periodsN2)/(1+numPeriods)

        meansRes[indicePasta][19] = sum(fileResDF[w]["Q"] for w in periodsN2)
        meansRes[indicePasta][20] = sum(fileResDF[w]["Ib"] for w in ["Week 0"]+periodsN2)/(1+numPeriods)
        meansRes[indicePasta][21] = sum(fileResDF[w]["Mb"] for w in periodsN2)
        meansRes[indicePasta][22] = sum(fileResDF[w]["PUbg"] for w in periodsN2)
        meansRes[indicePasta][23] = sum(fileResDF[w]["Mbg"] for w in periodsN2)
        meansRes[indicePasta][24] = sum(fileResDF[w]["PUbm"] for w in periodsN2)
        fileSolDF = pd.read_excel(pasta+"/Sols.xlsx",sheet_name="Sheet1",index_col=0,header=0)
        lastLine = len(fileSolDF.index)-1
        meansRes[indicePasta][25] = (fileSolDF["UB"][lastLine]-fileSolDF["LB"][lastLine])/fileSolDF["UB"][lastLine]
        
    meansResDF = pd.DataFrame(meansRes,index=None,columns=["CakeInv","Biogas","TurboGen","byProd","Biomet", "Value","EnergyBag", "EnergyBiogas","EnergyConsumed","Exported","FO","Time","Sugar","Etanol","M","Vin","FertVin","Cake","InvCake","Flair","Bagasse","BagCaldeira","BiogasProd","TurboGenBiogas","BiomethaneProd","Gap"])
    meansResDF.to_excel("saida.xlsx")
    colunas = ["EnergyBag", "EnergyBiogas","EnergyConsumed","Exported","FO","Time","Sugar","Etanol","M","Vin","FertVin","Cake","InvCake","Flair","Bagasse","BagCaldeira","BiogasProd","TurboGenBiogas","BiomethaneProd","Gap"]
    divisores = [1000,1000,1000,1000,1000000,1,1000,1000,1000,1000,1000,1000,1000,1000,1000,1000,1000000,1000000,1000000,0.01]
    for i in range(len(colunas)):
        createHeatChart(meansResDF,divisores[i],colunas[i], colunas[i]+".png","./graficos/")
    #createHeatChart(meansResDF,1000,"Exported", "Exported.png","./graficos/")
    #createHeatChart(meansResDF,1000,"Produced", "Produced.png","./graficos/")
    #createHeatChart(meansResDF,1000000,"FO", "FO.png","./graficos/")
    #createHeatChart(meansResDF,1,"Time", "Time.png","./graficos/")
    print(max(meansRes[i][25] for i in range(0,len(listaPastas))))


