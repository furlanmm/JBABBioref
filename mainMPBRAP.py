import pandas as pd
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.ticker as tk
import gurobipy_pandas as gppd
import gurobipy as gp
import itertools as it
import classPASUC as ps
import classMPBRAP as mr
import math
import sys
import os
import shutil
import time
from functools import partial

matplotlib.use("pgf")
matplotlib.rcParams.update({
    "pgf.texsystem": "pdflatex",
    'font.family': 'serif',
    'text.usetex': True,
    'pgf.rcfonts': False,
#    'mathtext.fontset': 'custom'
})

pd.options.mode.chained_assignment = None  # default='warn'

if __name__ == '__main__':

    processTimes = [('Start',time.time())]#Add process time value
    originalPasuc1 = ps.PASUCL1()
    RP2 = mr.MPBRAPL2()
    

    fileData = sys.argv[1]
    fileDataBio = sys.argv[2]
    resFolder = sys.argv[3]
    idTest = sys.argv[4]
    nameFolder = resFolder+"/"+idTest+"/"
    os.makedirs(nameFolder,exist_ok=True)
    #print("Time 1: ",time.time()-inicio)
    originalPasuc1.readDataCoop(fileData)
    #print("Time data L1: ",time.time()-inicio)
    processTimes += [('L1 data file',time.time())]#Add process time value
    modelPasuc = originalPasuc1.createModelL1("PASUC/N1")
    #print("Time model L1: ",time.time()-inicio)
    processTimes += [('L1 mathematical model',time.time())]#Add process time value
    modelPasuc.optimize()
    #print("Time solve L1: ",time.time()-inicio)
    processTimes += [('L1 solving process',time.time())]#Add process time value

    RP2.readDataBiodigester(fileDataBio)
    #print("Time read L2: ",time.time()-inicio)
    processTimes += [('L2 data file',time.time())]#Add process time value
    RP2.calcParamsL2(originalPasuc1)
    #print("Time data L2: ",time.time()-inicio)
    processTimes += [('L2 param construction',time.time())]#Add process time value
    modelPasuc2 = RP2.createModelL2("MPBRAP/N2",originalPasuc1)
    #print("Time model L2: ",time.time()-inicio)
    processTimes += [('L2 mathematical model',time.time())]#Add process time value
    #modelPasuc2.setParam("SolutionLimit",1)
    modelPasuc2.setParam("DualReductions",0)
    modelPasuc2.setParam("TimeLimit",600)
    modelPasuc2.optimize(RP2.data_cb)
    #print("Time solve L2: ",time.time()-inicio)
    processTimes += [('L2 solving process',time.time())]#Add process time value

    print(modelPasuc2.status)
    if modelPasuc2.status == gp.GRB.INFEASIBLE:
        modelPasuc2.computeIIS()
        modelPasuc2.feasRelaxS(1, True, True, True)
        modelPasuc2.setParam("TimeLimit",300)
        modelPasuc2.optimize(callback=RP2.data_cb)
        # vars = modelPasuc2.getVars()
        # ubpen = [1.0]*modelPasuc2.numVars
        # modelPasuc2.feasRelax(1, True, vars, None, ubpen, None, None)
        # modelPasuc2.optimize()
        # for aux in range(0,modelPasuc2.numVars):
        #     print(vars[aux],vars[aux].X)
        processTimes = [('Feasibility Process',time.time())]#Add process time value

    processTimes += [('End',time.time())]#Add process time value


    periodosGrafico = ["Week 0"]+RP2.periodsN2
    dadosSaida = pd.DataFrame(0.0,columns=periodosGrafico, index = ["K","M","Ts"]+["PU_"+p for p in RP2.ProductsMill]+["Im_"+p for p in RP2.ProductsMill]+["PUvin","Fvin","OMvin","PUcake","Icake","Mcake","Q","Ib","Mb","PUbg","Mbg","PUbm","EG","EM","EP","EB","EE","S","FO","Description","Time", "Sugar Revenue Value","Ethanol Revenue Value","Molasse Revenue Value","Biomethane Revenue Value","Diesel Substitution Value","Energy from Bagasse","Energy from Biogas","Raw Material Cost","Transport Cost","Process Cost","Backlog Cost","Episilon Parcel"])
    dadosSaida.loc["Mcake","Week 0"]=RP2.M0cake.X
    dadosSaida.loc["Icake","Week 0"]=RP2.I0cake
    dadosSaida.loc["Ts","Week 0"] = np.nan
    dadosSaida.loc["M","Week 0"] = np.nan
    dadosSaida.loc["PUcake","Week 0"] = np.nan
    dadosSaida.loc["PUbg","Week 0"] = np.nan
    dadosSaida.loc["Mbg","Week 0"] = np.nan
    dadosSaida.loc["PUbm","Week 0"] = np.nan
    dadosSaida.loc["Q","Week 0"] = np.nan
    dadosSaida.loc["EG","Week 0"] = np.nan
    dadosSaida.loc["EM","Week 0"] = np.nan
    dadosSaida.loc["EP","Week 0"] = np.nan
    dadosSaida.loc["EB","Week 0"] = np.nan
    dadosSaida.loc["EE","Week 0"] = np.nan
    dadosSaida.loc["FO","Week 0"] = RP2.modelL2.getObjective().getValue()
    #dadosSaida["Week 0"]["Time"] = RP2.modelL2.Runtime
    totalTime = processTimes[-1][1]-processTimes[0][1]
    dadosSaida.loc["Description","Week 0"] = 'Solution Time'
    dadosSaida.loc["Time","Week 0"] = totalTime
    dadosSaida.loc["Description","Week1"] = processTimes[0][0]
    dadosSaida.loc["Time","Week1"] = 0#processTimes[0][1]
    for i in range(1,len(processTimes)):
        dadosSaida.loc["Description","Week"+str(i)] = processTimes[i][0]
        dadosSaida.loc["Time","Week"+str(i)] = (processTimes[i][1]-processTimes[i-1][1])/totalTime

    for w in RP2.periodsSeason:
        dadosSaida.loc["K",w] = sum((RP2.process.index(k)+1)*RP2.x[k,w].X for k in RP2.process)
        dadosSaida.loc["M",w] = sum(RP2.M3[k,w].X for k in RP2.process)
        dadosSaida.loc["Ts",w] = sum(RP2.dataProcess.at[k,"Ts"]*RP2.x[k,w].X for k in RP2.process)
        for p in RP2.ProductsMill:
            dadosSaida.loc["PU_"+p,w] = sum(RP2.A.at[(p,k),w]*RP2.M3[k,w].X for k in RP2.process)    
            dadosSaida.loc["Im_"+p,w] = RP2.ImMill[p,w].X    
        dadosSaida.loc["PUvin",w] =sum(RP2.A.at[("Vinasse",k),w]*RP2.M3[k,w].X for k in RP2.process)
        dadosSaida.loc["Fvin",w] = sum(RP2.Fvin[k,w].X for k in RP2.process)
        dadosSaida.loc["OMvin",w] = sum(RP2.dataProcess.at[k,'OMvin']*(RP2.A.at[("Vinasse",k),w]*RP2.M3[k,w].X-RP2.Fvin[k,w].X) for k in RP2.process)
        dadosSaida.loc["PUcake",w] = sum(RP2.A.at[("Filter Cake",k),w]*RP2.M3[k,w].X for k in RP2.process)
        dadosSaida.loc["Icake",w] = RP2.Icake[w].X
        dadosSaida.loc["Mcake",w] = RP2.Mcake[w].X
        dadosSaida.loc["Q",w] = RP2.Q[w].X
        dadosSaida.loc["Ib",w] = RP2.Ib[w].X
        dadosSaida.loc["Mb",w] = RP2.Mb[w].X
        dadosSaida.loc["PUbg",w] = RP2.PUbg[w].X
        dadosSaida.loc["Mbg",w] = RP2.Mbg[w].X
        dadosSaida.loc["PUbm",w] = RP2.PUbm[w].X
        dadosSaida.loc["EG",w] = RP2.EG[w].X
        dadosSaida.loc["EM",w] = sum(RP2.cfe*RP2.M3[k,w].X for k in RP2.process)
        dadosSaida.loc["EP",w] = sum(RP2.dataProdMill.at[p,'cve']*RP2.PU2[p,w].X for p in RP2.ProductsMill)
        dadosSaida.loc["EB",w] = RP2.EB[w].X
        dadosSaida.loc["EE",w] = RP2.EE[w].X
        dadosSaida.loc["S",w] = RP2.S[w].X

    for w in range(RP2.season,RP2.numWeeks):
        dadosSaida.loc["Ts",RP2.periodsN2[w]] = np.nan
        dadosSaida.loc["Icake",RP2.periodsN2[w]] = RP2.Icake[RP2.periodsN2[w]].X
        dadosSaida.loc["Mcake",RP2.periodsN2[w]] = RP2.Mcake[RP2.periodsN2[w]].X
        dadosSaida.loc["Q",RP2.periodsN2[w]] = RP2.Q[RP2.periodsN2[w]].X
        dadosSaida.loc["Ib",RP2.periodsN2[w]] = RP2.Ib[RP2.periodsN2[w]].X
        dadosSaida.loc["Mb",RP2.periodsN2[w]] = RP2.Mb[RP2.periodsN2[w]].X
        dadosSaida.loc["PUbg",RP2.periodsN2[w]] = RP2.PUbg[RP2.periodsN2[w]].X
        dadosSaida.loc["Mbg",RP2.periodsN2[w]] = RP2.Mbg[RP2.periodsN2[w]].X
        dadosSaida.loc["PUbm",RP2.periodsN2[w]] = RP2.PUbm[RP2.periodsN2[w]].X
        dadosSaida.loc["EG",RP2.periodsN2[w]] = RP2.EG[RP2.periodsN2[w]].X
        dadosSaida.loc["EB",RP2.periodsN2[w]] = RP2.EB[RP2.periodsN2[w]].X
        dadosSaida.loc["EE",RP2.periodsN2[w]] = RP2.EE[RP2.periodsN2[w]].X
        dadosSaida.loc["S",RP2.periodsN2[w]] = RP2.S[RP2.periodsN2[w]].X        
    
    # Objective Function Parcels
    dadosSaida.loc["Sugar Revenue Value","Week 0"] = sum(RP2.PU2[pa,w].X/RP2.P2*RP2.R2*originalPasuc1.dataProd.at['ConvCrpaaa',pa] for w in RP2.periodsSeason for pa in RP2.SugarMill)
    dadosSaida.loc["Ethanol Revenue Value","Week 0"] = sum(RP2.PU2[pe,w].X/RP2.P3*RP2.R3*originalPasuc1.dataProd.at['ConvCrpaaa',pe] for w in RP2.periodsSeason for pe in RP2.EthanolMill)                          
    dadosSaida.loc["Molasse Revenue Value","Week 0"] = sum(RP2.PU2[pm,w].X/RP2.P4*RP2.R4*originalPasuc1.dataMills.at['Matr',RP2.millName] for w in RP2.periodsSeason for pm in RP2.MolassesMill)
    dadosSaida.loc["Biomethane Revenue Value","Week 0"] = sum(RP2.dataHorizon.at[w,"VbmPipe"]*RP2.dataHorizon.at[w,"DbmPipe"] for w in RP2.periodsN2)
    dadosSaida.loc["Diesel Substitution Value","Week 0"] = sum(RP2.dataHorizon.at[w,"VbmDiesel"]*(RP2.PUbm[w].X-RP2.dataHorizon.at[w,"DbmPipe"]) for w in RP2.periodsN2)
    dadosSaida.loc["Energy from Bagasse","Week 0"] = sum(RP2.ve*RP2.EE[w].X for w in RP2.periodsN2)
    dadosSaida.loc["Energy from Biogas","Week 0"] = sum(RP2.dataHorizon.at[w,"Vebio"]*RP2.EB[w].X for w in RP2.periodsN2)
    dadosSaida.loc["Raw Material Cost","Week 0"] = sum(RP2.Cm.at[m,w]*RP2.M1[m,w].X for m in RP2.sourceRaw for w in RP2.periodsSeason)
    dadosSaida.loc["Transport Cost","Week 0"] = sum(RP2.Cf.at[f,w]*RP2.M2[f,w].X for f in RP2.transport for w in RP2.periodsSeason)
    dadosSaida.loc["Process Cost","Week 0"]= sum(RP2.Ck.at[k,w]*RP2.M3[k,w].X for k in RP2.process for w in RP2.periodsSeason)
    dadosSaida.loc["Backlog Cost","Week 0"] = sum(originalPasuc1.ca*RP2.ImMill[p,w].X for p in RP2.ProductsMill for w in RP2.periodsSeason)
    dadosSaida.loc["Episilon Parcel","Week 0"] = -1e-4*RP2.maxEE.X
    dadosSaida.loc["Biogas Plant Cost","Week 0"] = sum(RP2.biogasPlantCost*RP2.PUbg[w].X for w in RP2.periodsSeason)
    dadosSaida.to_excel(nameFolder+"/Vars.xlsx")
    pd.DataFrame(RP2._data,index=None,columns=["Time","LB","UB"]).to_excel(nameFolder+"/Sols.xlsx")
    #print(RP2._data)
    #print(RP2.M0cake.X)

    '''
    ######################################################
    #Comentado para gerar apenas as tabelas de resultados#
    ###################################################### 
    
    dadosColuna = dadosSaida.T
    labelsEixoX=[0,4,8,12,16,20,24,28,32,36,40,44,48,52]
    pointsEixoX = [periodosGrafico[w] for w in labelsEixoX]

    divy = 10
    #tamanhofig = (7,2.125)
    tamanhofig = (7,2.8)
    plt.figure(figsize=tamanhofig)
    plt.grid(True,linestyle=':')
    plt.xlim(0,52)
    plt.xlabel("Periods (in Weeks)")
    ax1 = plt.gca()
    ax2 = ax1.twinx()
    ax1.set_ylabel("Sugarcane (in ton)")
    ax2.set_ylabel("Juice to Sugarmill (%)")
    #maxy = 55000
    maxy = math.ceil(dadosColuna["M"].max()/1000)*1000
    ax1.set_ylim(0,maxy)
    ax1.set_yticks([i*(maxy/divy) for i in range(0,divy+1)])
    maxy = 1
    ax2.set_ylim(0,maxy)
    ax2.set_yticks([i*(maxy/divy) for i in range(0,divy+1)],labels=[100*i*(maxy/divy) for i in range(0,divy+1)])
    lns1 = ax1.plot(periodosGrafico, dadosColuna["M"],label="Crushed Sugarcane",color="blue")
    lns2 = ax2.plot(periodosGrafico, dadosColuna["Ts"], color="red", label="Juice to SugarMill")
    lns = lns1 + lns2
    labels = [l.get_label() for l in lns]
    plt.xticks(ticks=pointsEixoX,labels=labelsEixoX)
    plt.legend(lns, labels, loc='upper right')
    plt.savefig(nameFolder+"/Sugarcane.pgf",bbox_inches='tight')
    #plt.show()
    plt.close()
    plt.figure(figsize=tamanhofig)
    plt.grid(True,linestyle=':')
    plt.xlim(0,52)
    #maxy = 15000
    maxy = max(math.ceil(dadosColuna["PUvin"].max()/1000)*1000, math.ceil(dadosColuna["Fvin"].max()/1000)*1000)
    plt.ylim(0,maxy)
    plt.yticks([i*(maxy/divy) for i in range(0,divy+1)])
    plt.plot(periodosGrafico, dadosColuna["PUvin"],label="Produced Vinasse",color="blue")
    plt.plot(periodosGrafico, dadosColuna["Fvin"],label="Fertirrigation", color="green")
    plt.plot(periodosGrafico, dadosColuna["PUvin"]-dadosColuna["Fvin"],label="Used Vinasse",color="red")
    plt.legend(loc='upper right')
    plt.xticks(ticks=pointsEixoX,labels=labelsEixoX)
    plt.xlabel("Periods (in Weeks)")
    plt.ylabel(r"Vinasse (in $m^3$)")
    plt.savefig(nameFolder+"/Vinasse.pgf",bbox_inches='tight')
    #plt.show()
    plt.close()
    plt.figure(figsize=tamanhofig)
    plt.grid(True,linestyle=':')
    plt.xlim(0,52)
    plt.xlabel("Periods (in Weeks)")
    ax1 = plt.gca()
    ax2 = ax1.twinx()
    ax1.set_ylabel("Filter Cake (in tons)")
    ax2.set_ylabel("Inventory Level (in tons)")
    #maxy = 1500
    maxy = max(math.ceil(dadosColuna["PUcake"].max()/100)*100, math.ceil(dadosColuna["Mcake"].max()/100)*100)
    ax1.set_ylim(0,maxy)
    ax1.set_yticks([i*(maxy/divy) for i in range(0,divy+1)])
    lns1 = ax1.plot(periodosGrafico, dadosColuna["PUcake"], label="Produced Filter Cake",color="blue")
    lns2 = ax1.plot(periodosGrafico, dadosColuna["Mcake"], label="Used Filter Cake",color="red")
    maxy = 16000
    maxy = math.ceil(dadosColuna["Icake"].max()/1000)*1000
    ax2.set_ylim(0,maxy)
    ax2.set_yticks([i*(maxy/divy) for i in range(0,divy+1)])
    lns3 = ax2.plot(periodosGrafico,dadosColuna["Icake"],label="Inventory of Filter Cake",color="purple")
    lns = lns1 + lns2 + lns3
    labels = [l.get_label() for l in lns]
    plt.xticks(ticks=pointsEixoX,labels=labelsEixoX)
    plt.legend(lns, labels, loc='upper right')
    plt.savefig(nameFolder+"/Cake.pgf",bbox_inches='tight')
    #plt.show()
    plt.close()
    plt.figure(figsize=tamanhofig)
    plt.grid(True,linestyle=':')
    plt.xlim(0,52)
    plt.xlabel("Periods (in Weeks)")
    ax1 = plt.gca()
    ax2 = ax1.twinx()
    ax1.set_ylabel(r"Biogas (in $m^3$)")
    ax2.set_ylabel(r"Biomethane (in $m^3$)")
    #maxy = 350000
    maxy = max(math.ceil(dadosColuna["PUbg"].max()/1000)*1000, math.ceil(dadosColuna["Q"].max()/1000)*1000, math.ceil(dadosColuna["Mbg"].max()/1000)*1000)
    ax1.set_ylim(0,maxy)
    ax1.set_yticks([i*(maxy/divy) for i in range(0,divy+1)])
    lns1 = ax1.plot(periodosGrafico,dadosColuna["PUbg"],label="Produced Biogas",color="blue")
    lns2 = ax1.plot(periodosGrafico,dadosColuna["Q"],label="Burned Biogas (flair)",color="green")
    lns3 = ax1.plot(periodosGrafico,dadosColuna["Mbg"],label="Burned Biogas (turbo gen)",color="red")
    #maxy = 85000
    maxy = math.ceil(dadosColuna["PUbm"].max()/900)*1000
    ax2.set_ylim(0,maxy)
    ax2.set_yticks([i*(maxy/divy) for i in range(0,divy+1)])
    lns4 = ax2.plot(periodosGrafico,dadosColuna["PUbm"],label="Produced Biomethane",color="deepskyblue")
    lns = lns1 + lns2 + lns3 + lns4
    labels = [l.get_label() for l in lns]
    plt.xticks(ticks=pointsEixoX,labels=labelsEixoX)
    plt.legend(lns, labels, loc='upper right')
    plt.savefig(nameFolder+"/Biogas.pgf",bbox_inches='tight')
    #plt.show()
    plt.close()

    plt.figure(figsize=tamanhofig)
    plt.grid(True,linestyle=':')
    plt.xlim(0,52)
    #maxy = 1700
    maxy = max(math.ceil(dadosColuna["EG"].max()/100)*100, math.ceil(dadosColuna["EB"].max()/100)*100)
    plt.ylim(0,maxy)
    plt.yticks([i*(maxy/divy) for i in range(0,divy+1)])
    plt.plot(periodosGrafico, dadosColuna["EG"],label="Boiler",color="blue")
    plt.plot(periodosGrafico, dadosColuna["EB"],label="Turbo Generator", color="green")
    plt.legend(loc='upper right')
    plt.xticks(ticks=pointsEixoX,labels=labelsEixoX)
    plt.xlabel("Periods (in Weeks)")
    plt.ylabel("Energy (in KWh)")
    plt.savefig(nameFolder+"/EnergyProd.pgf",bbox_inches='tight')
    #plt.show()
    plt.close()

    plt.figure(figsize=tamanhofig)
    plt.grid(True,linestyle=':')
    plt.xlim(0,52)
    #maxy = 2000
    maxy = max(math.ceil(dadosColuna["EM"].max()/100)*100, math.ceil(dadosColuna["EP"].max()/100)*100, math.ceil(dadosColuna["EE"].max()/100)*100, math.ceil(dadosColuna["EB"].max()/100)*100)
    plt.ylim(0,maxy)
    plt.yticks([i*(maxy/divy) for i in range(0,divy+1)])
    plt.plot(periodosGrafico, dadosColuna["EM"],label="Crushing", color="blue")
    plt.plot(periodosGrafico, dadosColuna["EP"],label="Production", color="red")
    plt.plot(periodosGrafico, dadosColuna["EE"],label="Export", color="green")
    plt.plot(periodosGrafico, dadosColuna["EB"],label="BioExport", color="limegreen")
    plt.legend(loc='upper right')
    plt.xticks(ticks=pointsEixoX,labels=labelsEixoX)
    plt.xlabel("Periods (in Weeks)")
    plt.ylabel("Energy (in KWh)")
    plt.savefig(nameFolder+"/EnergyUse.pgf",bbox_inches='tight')
    #plt.show()
    plt.close()



    dadosBarras = dadosSaida[dadosSaida.index.isin(['PU_Especial','PU_Extra','PU_VHP','PU_AEAC','PU_AEHC','PU_Melaco','Im_Especial','Im_Extra','Im_VHP','Im_AEAC','Im_AEHC','Im_Melaco'])].filter(items=RP2.periodsSeason).T
    # Configuração da largura das barras
    largura = 0.47
    x = np.arange(len(dadosBarras.index))+1
    # Criando o gráfico de barras
    #fig, ax = plt.subplots()
    plt.figure(figsize=tamanhofig)
    plt.xlim(0.5,25.5)
    plt.xticks(ticks=x,labels=x)
    ax = plt.gca()
    ax.xaxis.set_minor_locator(tk.AutoMinorLocator(2))
    #ax.xaxis.set_major_formatter(tk.NullFormatter())
    ax.xaxis.set_minor_formatter(tk.NullFormatter())
    ax.tick_params(axis='x', which='major', length=0)
    ax.tick_params(axis='x', which='minor', length=3)    
    plt.grid(visible=True,axis='y',linestyle=':')
    plt.grid(visible=True,axis='x',linestyle=':',which='minor')
    plt.xlabel("Periods (in Weeks)")
    maxy = max(math.ceil(max(dadosColuna["PU_Especial"]+dadosColuna["Im_Especial"])/90)*100, math.ceil(max(dadosColuna["PU_Extra"]+dadosColuna["Im_Extra"])/90)*100, math.ceil(max(dadosColuna["PU_VHP"]+dadosColuna["Im_VHP"])/90)*100, math.ceil(max(dadosColuna["PU_AEHC"]+dadosColuna["Im_AEHC"])/90)*100, math.ceil(max(dadosColuna["PU_AEAC"]+dadosColuna["Im_AEAC"])/90)*100)
    plt.ylim(0,maxy)
    plt.yticks([i*(maxy/divy) for i in range(0,divy+1)])
    barras1 = ax.bar(x - largura/2, dadosBarras['PU_Especial'], largura, label='Especial',color='orangered')
    barras2 = ax.bar(x - largura/2, dadosBarras['PU_Extra'], largura, label='Extra',bottom=dadosBarras['PU_Especial'],color='darkorange')
    barras3 = ax.bar(x - largura/2, dadosBarras['PU_VHP'], largura, label='VHP',bottom=dadosBarras['PU_Especial']+dadosBarras['PU_Extra'],color='gold')
    barras4 = ax.bar(x + largura/2, dadosBarras['PU_AEHC'], largura, label='AEHC',color='springgreen')
    barras5 = ax.bar(x + largura/2, dadosBarras['PU_AEAC'], largura, label='AEAC',bottom=dadosBarras['PU_AEHC'],color='green')
    #plt.xticks(ticks=x,labels=x)
    plt.ylabel(r"Prodution (in tons or $m^3$)")
    ax2 = ax.twinx()
    plt.ylim(maxy,0)
    plt.yticks([i*(maxy/divy) for i in range(0,divy+1)])
    barras6 = ax2.bar(x - largura/2, dadosBarras['Im_Especial'], largura, color='orangered')
    barras7 = ax2.bar(x - largura/2, dadosBarras['Im_Extra'], largura, bottom=dadosBarras['Im_Especial'],color='darkorange')
    barras8 = ax2.bar(x - largura/2, dadosBarras['Im_VHP'], largura, bottom=dadosBarras['Im_Especial']+dadosBarras['Im_Extra'],color='gold')
    barras9 = ax2.bar(x + largura/2, dadosBarras['Im_AEHC'], largura, color='springgreen')
    barras10 = ax2.bar(x + largura/2, dadosBarras['Im_AEAC'], largura, bottom=dadosBarras['Im_AEHC'],color='green')
    #ax.set_xticklabels(x)
    # ax.legend()
    #lns = barras1 + barras2 + barras3 + barras4 + barras5
    #labels = [l.get_label() for l in lns]
    #plt.legend(lns, labels, loc='upper right')
    ax.legend(loc='upper center', bbox_to_anchor=(0.5, 1.21), ncol=5, labelspacing=1, columnspacing=1)
    plt.ylabel(r"Backlog (in tons or $m^3$)")
    # Exibindo o gráfico
    plt.savefig(nameFolder+"/Production.pgf",bbox_inches='tight')
    plt.close()

    shutil.copy("baseSaida.tex", nameFolder+"/saida.tex")
    #os.system("pdflatex "+nameFolder+"/saida.tex")
    os.chdir(nameFolder)
    os.system("pdflatex saida.tex")
 
    #print(math.ceil(dadosColuna["EM"].max()/100)*100)
    print(processTimes)
    '''
