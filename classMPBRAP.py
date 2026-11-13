import pandas as pd
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.ticker as tk
import gurobipy_pandas as gppd
import gurobipy as gp
import itertools as it
import classPASUC as ps
import time


class MPBRAPL2(ps.PASUCL2):
    #Internal Data between models
    Rb: pd.DataFrame = ...
    #Data from mill
    BiodigesterMill: list = ...
    periodsSeason: list = ...
    #txvin: pd.DataFrame = ...
    #txcake:pd.DataFrame = ...
    Epbcake: float = ...
    txbiogasvin: float = ...
    txbiogascake: float = ...
    #tconvvin: float = ...
    tconvcake: float = ...
    BiogasCake: float = ...
    #OMvin: pd.DataFrame = ... #Organic Matter from vinasse acording to process decisions
    digmaxvin: float = ...
    deltaOMvin: int = ...
    deltacake: int = ...
    digmaxcake: float = ...
    tgmax: float = ...
    txbiomet: float = ...
    rcftg: float = ...
    I0cake: float = ...
    biogasPlantCost: float = ...
    season: int = ...
    #MPBRAPN2 Variables
    Fvin: gp.Var = ...
    Icake: gp.Var = ...
    Mcake: gp.Var = ...
    Q: gp.Var = ...
    PUbg: gp.Var = ...
    PUbm: gp.Var = ...
    Mbg: gp.Var = ...
    M0cake: gp.Var = ...
    EB: gp.Var = ...

    maxEE: gp.Var = ...

    _obj = None
    _bd = None
    _data = []
    _start = time.time()


    #Callback function to save bounds during optimize
    def data_cb(self, model, where):
        #print(time.time() - self._start, model, where)
        if where == gp.GRB.Callback.MIP:
            cur_obj = model.cbGet(gp.GRB.Callback.MIP_OBJBST)
            cur_bd = model.cbGet(gp.GRB.Callback.MIP_OBJBND)
            # Did objective value or best bound change?
            #print(time.time() - self._start, cur_obj, cur_bd)
            if self._obj != cur_obj or self._bd != cur_bd:
                self._obj = cur_obj
                self._bd = cur_bd
                self._data.append([time.time() - self._start, cur_obj, cur_bd])


    def readDataBiodigester(self, fileName):
        self.readDataMill(fileName)
        self.BiodigesterMill = self.generalN2.filter(like="BiodigesterMill", axis = 0).dropna(axis=1).values.tolist()[0]
        self.bioDig = pd.read_excel(fileName, sheet_name = "Biodigester",index_col=0, header=None, comment="#")
        self.Epbcake = self.bioDig.filter(like="EpbCake", axis = 0).dropna(axis=1).values.tolist()[0][0] 
        self.txbiogasvin = self.bioDig.filter(like="txBiogasVin", axis = 0).dropna(axis=1).values.tolist()[0][0] 
        self.txbiogascake = self.bioDig.filter(like="txBiogasCake", axis = 0).dropna(axis=1).values.tolist()[0][0] 
        self.digmaxvin = self.bioDig.filter(like="digMaxVin", axis = 0).dropna(axis=1).values.tolist()[0][0] 
        self.digmaxcake = self.bioDig.filter(like="digMaxCake", axis = 0).dropna(axis=1).values.tolist()[0][0] 
        self.deltaOMvin = self.bioDig.filter(like="deltaOMvin", axis = 0).dropna(axis=1).values.tolist()[0][0] 
        self.deltacake = self.bioDig.filter(like="deltacake", axis = 0).dropna(axis=1).values.tolist()[0][0] 
        self.tgmax = self.bioDig.filter(like="tgmax", axis = 0).dropna(axis=1).values.tolist()[0][0]
        self.txbiomet = self.bioDig.filter(like="txbiomet", axis = 0).dropna(axis=1).values.tolist()[0][0]
        self.rcftg = self.bioDig.filter(like="rcftg", axis = 0).dropna(axis=1).values.tolist()[0][0]
        #self.tconvvin = self.bioDig.filter(like="txConvVin", axis = 0).dropna(axis=1).values.tolist()[0][0]
        self.tconvcake = self.bioDig.filter(like="txConvCake", axis = 0).dropna(axis=1).values.tolist()[0][0]
        self.I0cake = self.bioDig.filter(like="I0Cake", axis = 0).dropna(axis=1).values.tolist()[0][0]
        self.biogasPlantCost = self.bioDig.filter(like="biogasPlantCost", axis = 0).dropna(axis=1).values.tolist()[0][0]
        self.season = self.generalN2.filter(like="Season", axis = 0).dropna(axis=1).values.tolist()[0][0]
        self.periodsSeason = ['Week'+str(t) for t in range(1,self.season+1)]
        #self.txvin = self.dataProcess['txvin']
        #self.txcake = self.dataProcess['txcake']
        #self.OMvin = self.dataProcess['OMvin']

    def calcParamsL2(self,modelL1):
        #This function uses data from Mill to calculate part of the parameters from PASUC/L2 
        #Begin of productive matrix values calculation (A_pkw, p:products, k:processes, w:weeks)
        inicio = time.time()
        arraySjM = np.zeros((len(self.process),len(self.periodsSeason)),np.float64)
        for k in self.process:
            for w in self.periodsSeason:
                arraySjM[self.process.index(k)][self.periodsSeason.index(w)] = self.dataProcess.at[k,'PzaS']*(self.dataHorizon.at[w,'PzaJ']-self.PzaM)/ (self.dataHorizon.at[w,'PzaJ']*(self.dataProcess.at[k,'PzaS']-self.PzaM))
        self.SjM = pd.DataFrame(arraySjM,index=self.process,columns=self.periodsSeason)  
        self.SmF = (self.Bmf*self.PzaM)/100
        arrayRd = np.zeros((len(self.EthanolMill),len(self.periodsSeason)),np.float64)
        for pe in self.EthanolMill:
            for w in self.periodsSeason:
                arrayRd[self.EthanolMill.index(pe)][self.periodsSeason.index(w)] = (self.dataProdMill.at[pe,'Rend']*self.dataHorizon.at[w,'Efd'])/100  
        self.Rd = pd.DataFrame(arrayRd,index=self.EthanolMill,columns=self.periodsSeason)
        arrayRa = np.zeros((len(self.SugarMill),len(self.process),len(self.periodsSeason)),np.float64)
        arrayRe = np.zeros((len(self.EthanolMill),len(self.process),len(self.periodsSeason)),np.float64)
        arrayRm = np.zeros((len(self.MolassesMill),len(self.process),len(self.periodsSeason)),np.float64)
        arrayRb = np.zeros((len(self.BiodigesterMill),len(self.process),len(self.periodsSeason)),np.float64)
        #arrayA = np.zeros((len(self.ProductsMill+["Vinasse","Filter Cake","VinasseBiogas"])*len(self.process),len(self.periodsSeason)),np.float64)
        # for w in self.periodsSeason:
        #     for k in self.process:
        #         ind = 0
        #         for pa in self.SugarMill:
        #             arrayA[ind*self.numProc+self.process.index(k)][self.periodsSeason.index(w)] = self.dataHorizon.at[w,'pc']*10*(self.dataHorizon.at[w,'Elbti']/100)*self.SjM.at[k,w]*self.dataProcess.at[k,'Ts']*self.dataProcess.at[k,pa]
        #             ind+=1
        #         for pe in self.EthanolMill:
        #             arrayA[ind*self.numProc+self.process.index(k)][self.periodsSeason.index(w)] = (self.dataHorizon.at[w,'Elbti']/10)*(self.Rd.at[pe,w]/100)*(((self.dataHorizon.at[w,'pc']/0.95)*(1-self.SjM.at[k,w])+self.dataHorizon.at[w,'Ar'])*self.dataProcess.at[k,'Ts']*self.dataProcess.at[k,'Tm']+
        #                                                                                                                               (self.dataHorizon.at[w,'pc']/0.95+self.dataHorizon.at[w,'Ar'])*(1-self.dataProcess.at[k,'Ts']))*self.dataProcess.at[k,pe]
        #             ind+=1
        #         for pm in self.MolassesMill:
        #             arrayA[ind*self.numProc+self.process.index(k)][self.periodsSeason.index(w)] = self.dataHorizon.at[w,'pc']*self.dataHorizon.at[w,'Elbti']*(1-self.SjM.at[k,w])*self.dataProcess.at[k,'Ts']*100*(1-self.dataProcess.at[k,'Tm'])/(10*self.SmF)
        #             ind+=1
        #         arrayA[ind*self.numProc+self.process.index(k)][self.periodsSeason.index(w)] = sum(self.dataProcess.at[k,'txvin']*arrayA[(self.EthanolMill.index(pe)+len(self.SugarMill))*self.numProc+self.process.index(k)][self.periodsSeason.index(w)] for pe in self.EthanolMill)
        #         ind+=1
        #         arrayA[ind*self.numProc+self.process.index(k)][self.periodsSeason.index(w)] = self.dataProcess.at[k,'txcake']*self.dataProcess.at[k,'Tc']*self.dataHorizon.at[w,'txJuice']
        #         ind+=1
        #         arrayA[ind*self.numProc+self.process.index(k)][self.periodsSeason.index(w)] = 1000*self.txbiogasvin*self.dataProcess.at[k,'OMvin']*self.dataHorizon.at[w,'txConvVin']
        for k in self.process:
            for w in self.periodsSeason:
                for pa in self.SugarMill:
                    arrayRa[self.SugarMill.index(pa)][self.process.index(k)][self.periodsSeason.index(w)] = self.dataHorizon.at[w,'pc']*10*(self.dataHorizon.at[w,'Elbti']/100)*self.SjM.at[k,w]*self.dataProcess.at[k,'Ts']*self.dataProcess.at[k,pa]
                for pe in self.EthanolMill:
                    arrayRe[self.EthanolMill.index(pe)][self.process.index(k)][self.periodsSeason.index(w)] = (self.dataHorizon.at[w,'Elbti']/10)*(self.Rd.at[pe,w]/100)*(((self.dataHorizon.at[w,'pc']/0.95)*(1-self.SjM.at[k,w])+self.dataHorizon.at[w,'Ar'])*self.dataProcess.at[k,'Ts']*self.dataProcess.at[k,'Tm']+
                                                                                                                                            (self.dataHorizon.at[w,'pc']/0.95+self.dataHorizon.at[w,'Ar'])*(1-self.dataProcess.at[k,'Ts']))*self.dataProcess.at[k,pe]                                                                                                            
                for pm in self.MolassesMill:
                    arrayRm[self.MolassesMill.index(pm)][self.process.index(k)][self.periodsSeason.index(w)] = self.dataHorizon.at[w,'pc']*self.dataHorizon.at[w,'Elbti']*(1-self.SjM.at[k,w])*self.dataProcess.at[k,'Ts']*100*(1-self.dataProcess.at[k,'Tm'])/(10*self.SmF)
                arrayRb[self.BiodigesterMill.index("Vinasse")][self.process.index(k)][self.periodsSeason.index(w)] = sum(self.dataProcess.at[k,'txvin']*arrayRe[self.EthanolMill.index(pe)][self.process.index(k)][self.periodsSeason.index(w)] for pe in self.EthanolMill)
                arrayRb[self.BiodigesterMill.index("Filter Cake")][self.process.index(k)][self.periodsSeason.index(w)] = self.dataProcess.at[k,'txcake']*self.dataProcess.at[k,'Tc']*self.dataHorizon.at[w,'txJuice']
                arrayRb[self.BiodigesterMill.index("VinasseBiogas")][self.process.index(k)][self.periodsSeason.index(w)] = 1000*self.txbiogasvin*self.dataProcess.at[k,'OMvin']*self.dataHorizon.at[w,'txConvVin']
                #arrayRb[self.BiodigesterMill.index("CakeBiogas")][self.process.index(k)][self.periodsSeason.index(w)] = 1000*self.txbiogascake*self.tconvcake
        self.Ra = pd.DataFrame([[arrayRa[self.SugarMill.index(pa)][self.process.index(k)][self.periodsSeason.index(w)] for w in self.periodsSeason] for pa in self.SugarMill for k in self.process],index=it.product(self.SugarMill,self.process),columns=self.periodsSeason)
        self.Re = pd.DataFrame([[arrayRe[self.EthanolMill.index(pe)][self.process.index(k)][self.periodsSeason.index(w)] for w in self.periodsSeason] for pe in self.EthanolMill for k in self.process],index=it.product(self.EthanolMill,self.process),columns=self.periodsSeason)
        self.Rm = pd.DataFrame([[arrayRm[self.MolassesMill.index(pm)][self.process.index(k)][self.periodsSeason.index(w)] for w in self.periodsSeason] for pm in self.MolassesMill for k in self.process],index=it.product(self.MolassesMill,self.process),columns=self.periodsSeason)
        self.Rb = pd.DataFrame([[arrayRb[self.BiodigesterMill.index(pb)][self.process.index(k)][self.periodsSeason.index(w)] for w in self.periodsSeason] for pb in self.BiodigesterMill for k in self.process],index=it.product(self.BiodigesterMill,self.process),columns=self.periodsSeason)
        self.A = pd.concat([self.Ra, self.Re, self.Rm, self.Rb]).div(1000)
        #self.A.to_excel("A.xlsx")
        # self.A = pd.DataFrame(arrayA,index=it.product(self.ProductsMill+["Vinasse","Filter Cake","VinasseBiogas"],self.process),columns=self.periodsSeason).div(1000)
        #End of productive matrix values
        #Begin of industrial cost of TRS produced in process k in week w  calculation(Ck)
        arrayConvA = np.zeros((len(self.SugarMill),len(self.process),len(self.periodsSeason)),np.float64)
        arrayConvD = np.zeros((len(self.EthanolMill),len(self.process),len(self.periodsSeason)),np.float64)
        arrayConvM = np.zeros((len(self.MolassesMill),len(self.process),len(self.periodsSeason)),np.float64)
        for k in self.process:
            for w in self.periodsSeason:
                for pa in self.SugarMill:
                    arrayConvA[self.SugarMill.index(pa)][self.process.index(k)][self.periodsSeason.index(w)] = self.Ra.at[(pa,k),w]*(self.dataProcess.at[k,'Pol']/100)*(1-self.dataProcess.at[k,'Umid']/100)*self.dataProcess.at[k,pa]
                for pe in self.EthanolMill:
                    arrayConvD[self.EthanolMill.index(pe)][self.process.index(k)][self.periodsSeason.index(w)] = self.Re.at[(pe,k),w]*self.dataProdMill.at[pe,'Fator']*self.dataProcess.at[k,pe]
                for pm in self.MolassesMill:
                    arrayConvM[self.MolassesMill.index(pm)][self.process.index(k)][self.periodsSeason.index(w)] = (self.Rm.at[(pm,k),w]*self.SmF/100)/0.95+self.Rm.at[(pm,k),w]*self.arm/100
        self.ConvA = pd.DataFrame([[arrayConvA[self.SugarMill.index(pa)][self.process.index(k)][self.periodsSeason.index(w)] for w in self.periodsSeason] for pa in self.SugarMill for k in self.process],index=it.product(self.SugarMill,self.process),columns=self.periodsSeason)
        self.ConvD = pd.DataFrame([[arrayConvD[self.EthanolMill.index(pe)][self.process.index(k)][self.periodsSeason.index(w)] for w in self.periodsSeason] for pe in self.EthanolMill for k in self.process],index=it.product(self.EthanolMill,self.process),columns=self.periodsSeason)
        self.ConvM = pd.DataFrame([[arrayConvM[self.MolassesMill.index(pm)][self.process.index(k)][self.periodsSeason.index(w)] for w in self.periodsSeason] for pm in self.MolassesMill for k in self.process],index=it.product(self.MolassesMill,self.process),columns=self.periodsSeason)
        self.ConvAll = pd.concat([self.ConvA, self.ConvD, self.ConvM])
        print("Conv ",time.time()-inicio)
        arrayARTproc = np.zeros((len(self.process),len(self.periodsSeason)),np.float64)
        arrayCustoEnsaque = np.zeros((len(self.process),len(self.periodsSeason)),np.float64)
        for k in self.process:
            for w in self.periodsSeason:
                arrayARTproc[self.process.index(k)][self.periodsSeason.index(w)] = sum(self.ConvAll.at[(p,k),w]/self.dataProdMill.at[p,'divART'] for p in self.ProductsMill)
                arrayCustoEnsaque[self.process.index(k)][self.periodsSeason.index(w)] = sum((self.Ra.at[(pa,k),w]*self.dataProcess.at[k,'CustoEnsaque'])/50.0 for pa in self.SugarMill)
        self.ARTproc = pd.DataFrame(arrayARTproc,index=self.process,columns=self.periodsSeason)
        self.CustoEnsaque = pd.DataFrame(arrayCustoEnsaque,index=self.process,columns=self.periodsSeason)
        self.Ck = self.ARTproc.mul(self.cproc).add(self.CustoEnsaque)
        #Begin of cost of TRS from supplier m in week w  calculation(Cm)
        arrayCm = np.zeros((len(self.sourceRaw),len(self.periodsSeason)),np.float64)
        for m in self.sourceRaw:
            for w in self.periodsSeason:
                arrayCm[self.sourceRaw.index(m)][self.periodsSeason.index(w)] = self.atrMill.at[m,w]*self.dataHorizon.at[w,'vAtr']+self.dataSourceRaw.at[m,'delta']
        self.Cm = pd.DataFrame(arrayCm,index=self.sourceRaw,columns=self.periodsSeason)
        #Input calculation for PASUC/L2 from outputs of PASUC/L1
        self.P1 = sum(modelL1.PU[p,u,t].X*modelL1.dataProd.at["Patr",p] for p in modelL1.products for u in modelL1.mills for t in modelL1.periods)
        self.P2 = sum(modelL1.PU[pa,u,t].X*modelL1.dataProd.at["ConvCrpaaa",pa] for pa in modelL1.sugar for u in modelL1.mills for t in modelL1.periods)
        self.P3 = sum(modelL1.PU[pd,u,t].X*modelL1.dataProd.at["ConvCrpaaa",pd] for pd in modelL1.ethanol for u in modelL1.mills for t in modelL1.periods)
        self.P4 = sum(modelL1.PU[pm,u,t].X*modelL1.dataMills.at["Matr",u] for pm in modelL1.molasses for u in modelL1.mills for t in modelL1.periods)
        self.R1 = sum(modelL1.vp.at[p,t]*modelL1.D[p,t].X for p in modelL1.products for t in modelL1.periods)
        self.R2 = sum(modelL1.vp.at[pa,t]*modelL1.D[pa,t].X for pa in modelL1.sugar for t in modelL1.periods)
        self.R3 = sum(modelL1.vp.at[pd,t]*modelL1.D[pd,t].X for pd in modelL1.ethanol for t in modelL1.periods)
        self.R4 = 0 #sum(modelL1.vp.at[pm,t]*modelL1.D[pm,t].X for pm in modelL1.molasses for t in modelL1.periods)
        arrayGoals = np.zeros((len(self.ProductsMill),len(self.periodsSeason)),dtype=np.float64)
        #Goals for mill self.millName
        for i in self.periodsGoals.index:
            for p in self.ProductsMill:
                arrayGoals[self.ProductsMill.index(p)][self.periodsSeason.index(self.periodsGoals.at[i,"Week"])] = modelL1.PU[p,self.millName,self.periodsGoals.at[i,"Month"]].X
        self.Goals = pd.DataFrame(arrayGoals,index=self.ProductsMill,columns=self.periodsSeason)

    def createModelL2(self,name,modelL1):
        #This function create the PASUC/L2 model
        self.modelL2 = gp.Model(name)
        #Variables instances
        self.x = self.modelL2.addVars(self.process,self.periodsSeason,name="Xkt",vtype=gp.GRB.BINARY)
        self.MMill = self.modelL2.addVars(self.periodsSeason,name="Mt")
        self.M1 = self.modelL2.addVars(self.sourceRaw,self.periodsSeason,name="M1mt")
        self.M2 = self.modelL2.addVars(self.transport,self.periodsSeason,name="M2ft")
        self.M3 = self.modelL2.addVars(self.process,self.periodsSeason,name="M3kt")
        self.PU2 = self.modelL2.addVars(self.ProductsMill,self.periodsSeason,name="PU2")
        self.C = self.modelL2.addVars(self.sourceRaw,self.periodsSeason,name="Cmt")
        self.S = self.modelL2.addVars(self.periodsN2,name="St")
        self.ImMill = self.modelL2.addVars(self.ProductsMill,self.periodsSeason,name="I-pt")
        self.Ib = self.modelL2.addVars(self.periodsN2,name="Ibt")
        self.Mb = self.modelL2.addVars(self.periodsN2,name="Mbt")
        self.VAP = self.modelL2.addVars(self.periodsN2,name="VAPt")
        self.EG = self.modelL2.addVars(self.periodsN2,name="EGt")
        self.EE = self.modelL2.addVars(self.periodsN2,name="EEt")
        self.EB = self.modelL2.addVars(self.periodsN2,name="EB")        
        #Variable used to smooth energy exportation
        self.maxEE = self.modelL2.addVar(name="maxEE")
        #Biodigester variables
        self.Fvin = self.modelL2.addVars(self.process,self.periodsSeason,name="Fvin")
        self.Icake = self.modelL2.addVars(self.periodsN2,name="Icake")
        self.Mcake = self.modelL2.addVars(self.periodsN2,name="Mcake")
        self.Q = self.modelL2.addVars(self.periodsN2,name="Q")
        self.PUbg = self.modelL2.addVars(self.periodsN2,name="PUbg")
        self.PUbm = self.modelL2.addVars(self.periodsN2,name="PUbm")
        self.Mbg = self.modelL2.addVars(self.periodsN2,name="Mbg")
        self.M0cake = self.modelL2.addVar(name="M0cake")
        self.BiogasCake = self.txbiogascake*self.tconvcake/self.txbiomet
        #End variables
        self.modelL2.update()
        #self.calcParamsL2(modelL1)
        self.modelL2.addConstrs((gp.quicksum(self.A.at[(p,k),w]*self.M3[k,w] for k in self.process) == self.PU2[p,w]  for w in self.periodsSeason for p in self.ProductsMill),name="constPU")
        self.modelL2.addConstrs((gp.quicksum(self.PU2[p,self.periodsSeason[wl]] for wl in range(0,w+1))+self.ImMill[p,self.periodsSeason[w]]
                            -gp.quicksum(self.Goals.at[p,self.periodsSeason[wl]] for wl in range(0,w+1))>=0
                            for p in self.ProductsMill for w in range(0,self.season)),name="const23")
        self.modelL2.addConstrs((gp.quicksum(self.M1[m,w] for m in self.sourceRaw) == self.MMill[w] for w in self.periodsSeason),name="const24a")
        self.modelL2.addConstrs((gp.quicksum(self.M2[f,w] for f in self.transport) == self.MMill[w] for w in self.periodsSeason),name="const24b")
        self.modelL2.addConstrs((gp.quicksum(self.M3[k,w] for k in self.process) == self.MMill[w] for w in self.periodsSeason),name="const24c")
        self.modelL2.addConstrs((self.C[m,self.periodsSeason[w]]==self.C[m,self.periodsSeason[w-1]]-self.M1[m,self.periodsSeason[w-1]] for m in self.sourceRaw for w in range(1,self.season)),name="const25")
        self.modelL2.addConstrs((self.C[m,self.periodsSeason[0]]== self.dataSourceRaw.at[m,'C0'] for m in self.sourceRaw),name="const25_0")
        self.modelL2.addConstrs((self.C[m,w]>=self.M1[m,w] for m in self.sourceRaw for w in self.periodsSeason),name="const25b")
        self.modelL2.addConstr((gp.quicksum(self.dataSourceRaw.at[m,'C0'] for m in self.sourceRaw)==gp.quicksum(self.M1[m,w] for m in self.sourceRaw for w in self.periodsSeason)),name="const26")
        self.modelL2.addConstrs((gp.quicksum(self.C[mp,w] for mp in self.selfRaw)<=self.dataHorizon.at[w,'disp'] for w in self.periodsSeason),name="const27")
        self.modelL2.addConstrs((gp.quicksum(self.M1[mf,w] for mf in self.externalRaw)<= self.dataHorizon.at[w,'alfa']*self.MMill[w] for w in self.periodsSeason),name="const28")
        self.modelL2.addConstrs((self.MMill[w] >= modelL1.dataMills.at["Mmin",self.millName]*self.dataHorizon.at[w,'n']*self.dataHorizon.at[w,'phi'] for w in self.periodsSeason),name="constr29a")
        self.modelL2.addConstrs((self.MMill[w] <= modelL1.dataMills.at["Mmax",self.millName]*self.dataHorizon.at[w,'n']*self.dataHorizon.at[w,'phi'] for w in self.periodsSeason),name="constr29b")
        #print(self.dataHorizon)
        self.modelL2.addConstrs((self.S[self.periodsSeason[w]]==self.S[self.periodsSeason[w-1]]+(gp.quicksum((self.dataProdMill.at[p,'rca']+self.dataProdMill.at[p,'Patr']*self.rcb)*self.PU2[p,self.periodsSeason[w]] for p in self.ProductsMill)+self.dataHorizon.at[self.periodsSeason[w],'rcc'])
                        -(gp.quicksum(self.Gm.at[m,self.periodsSeason[w]]*self.M1[m,self.periodsSeason[w]] for m in self.sourceRaw)+gp.quicksum(self.Gf.at[f,self.periodsSeason[w]]*self.M2[f,self.periodsSeason[w]] for f in self.transport)
                        +gp.quicksum(self.Gk.at[k,self.periodsSeason[w]]*self.M3[k,self.periodsSeason[w]] for k in self.process)+self.dataHorizon.at[self.periodsSeason[w],'Gfix']) for w in range(1,self.season)),name="const30Season")
        self.modelL2.addConstr((self.S[self.periodsSeason[0]] == self.initBudget+(gp.quicksum((self.dataProdMill.at[p,'rca']+self.dataProdMill.at[p,'Patr']*self.rcb)*self.PU2[p,self.periodsSeason[0]] for p in self.ProductsMill)+self.dataHorizon.at[self.periodsSeason[0],'rcc'])
                        -(gp.quicksum(self.Gm.at[m,self.periodsSeason[0]]*self.M1[m,self.periodsSeason[0]] for m in self.sourceRaw)+gp.quicksum(self.Gf.at[f,self.periodsSeason[0]]*self.M2[f,self.periodsSeason[0]] for f in self.transport)
                        +gp.quicksum(self.Gk.at[k,self.periodsSeason[0]]*self.M3[k,self.periodsSeason[0]] for k in self.process)+self.dataHorizon.at[self.periodsSeason[0],'Gfix'])),name="const30_0Season")
        self.modelL2.addConstrs((self.S[self.periodsN2[w]]==self.S[self.periodsN2[w-1]]+self.dataHorizon.at[self.periodsN2[w],'rcc']
                        -self.dataHorizon.at[self.periodsN2[w],'Gfix'] for w in range(self.season,self.numWeeks)),name="const30OffSeason")
        self.modelL2.addConstrs((self.M2[f,w]<=self.dataTransport.at[f,'fmax']*self.dataHorizon.at[w,'n']*self.beta.at[f,w] for f in self.transport for w in self.periodsSeason),name="const31")
        self.modelL2.addConstrs((self.PU2[p,w]<=modelL1.cpmax.at[p,self.millName]*self.dataHorizon.at[w,'n'] for p in self.ProductsMill for w in self.periodsSeason),name="const32")
        self.modelL2.addConstrs((gp.quicksum(self.PU2[pa,w] for pa in self.SugarMill) <= modelL1.dataMills.at['Cfmax',self.millName]*self.dataHorizon.at[w,'n'] for w in self.periodsSeason),name="const33")
        self.modelL2.addConstrs((gp.quicksum(self.PU2[pe,w] for pe in self.EthanolMill) <= modelL1.dataMills.at['Cdmax',self.millName]*self.dataHorizon.at[w,'n'] for w in self.periodsSeason),name="const34")
        self.modelL2.addConstrs((gp.quicksum(self.x[k,w] for k in self.process)==1 for w in self.periodsSeason),name="const35")
        self.modelL2.addConstrs((self.M3[k,w]<=modelL1.dataMills.at["Mmax",self.millName]*self.dataHorizon.at[w,'n']*self.dataHorizon.at[w,'phi']*self.x[k,w] for k in self.process for w in self.periodsSeason),name="const36")
        self.modelL2.addConstrs((self.Ib[self.periodsSeason[w]]==self.Ib[self.periodsSeason[w-1]]+gp.quicksum(self.fiber.at[m,self.periodsSeason[w]]/(1-self.dataHorizon.at[self.periodsSeason[w],'ub'])*self.M1[m,self.periodsSeason[w]] for m in self.sourceRaw)-self.Mb[self.periodsSeason[w]] for w in range(1,self.season)),name="const37Season")
        self.modelL2.addConstrs((self.Ib[self.periodsN2[w]]==self.Ib[self.periodsN2[w-1]]-self.Mb[self.periodsN2[w]] for w in range(self.season,self.numWeeks)),name="const37Season")
        self.modelL2.addConstr((self.Ib[self.periodsSeason[0]]==gp.quicksum(self.fiber.at[m,self.periodsSeason[0]]/(1-self.dataHorizon.at[self.periodsSeason[0],'ub'])*self.M1[m,self.periodsSeason[0]] for m in self.sourceRaw)-self.Mb[self.periodsSeason[0]]),name="const37_0")
        self.modelL2.addConstrs((self.Ib[w]>=gp.quicksum(self.fiber.at[m,w]/(1-self.dataHorizon.at[w,'ub'])*self.M1[m,w] for m in self.sourceRaw)*self.eb for w in self.periodsSeason),name="const38")
        self.modelL2.addConstr((self.Ib[self.periodsN2[-1]]>=self.Epb),name="const39")
        self.modelL2.addConstrs((self.rc*self.Mb[w]==self.VAP[w] for w in self.periodsN2),name="const40")
        self.modelL2.addConstrs((self.VAP[w] >= gp.quicksum(self.M3[k,w] for k in self.process)*self.cfvap + (self.EG[w]/self.rcf) for w in self.periodsSeason),name="const41aSeason")
        self.modelL2.addConstrs((self.VAP[self.periodsN2[w]] >= self.EG[self.periodsN2[w]]/self.rcf for w in range(self.season,self.numWeeks)),name="const41aOffSeason")
        self.modelL2.addConstrs((gp.quicksum(self.M3[k,w] for k in self.process)*self.cfvap + (self.EG[w]/self.rcf) >= gp.quicksum(self.dataProdMill.at[p,'cvap']*self.PU2[p,w] for p in self.ProductsMill) for w in self.periodsSeason),name="const41b")
        #self.modelL2.addConstrs((self.EG[w]-(gp.quicksum(self.cfe*self.M3[k,w] for k in self.process)+gp.quicksum(self.dataProdMill.at[p,'cve']*self.PU2[p,w] for p in self.ProductsMill))==self.EE[w] for w in self.periodsN2),name="const42")
        self.modelL2.addConstrs((self.VAP[w]<=self.vapmax*self.dataHorizon.at[w,'n'] for w in self.periodsN2),name="const43")
        self.modelL2.addConstrs((self.EG[w]<=self.egmax*self.dataHorizon.at[w,'n'] for w in self.periodsN2),name="const44")
        self.modelL2.addConstrs((self.A.at[("Vinasse",k),w]*self.M3[k,w]<=self.digmaxvin*self.dataHorizon.at[w,'nbio']+self.Fvin[k,w] for k in self.process for w in self.periodsSeason),name="bd1")
        self.modelL2.addConstrs((self.A.at[("Vinasse",k),w]*self.M3[k,w]-self.Fvin[k,w]>=0 for k in self.process for w in self.periodsSeason),name="bdaux2")
        #self.modelL2.addConstrs((self.Fvin[k,w]<=self.A.at[("Vinasse",k),w]*modelL1.dataMills.at["Mmax",self.millName]*self.dataHorizon.at[w,'n']*self.dataHorizon.at[w,'phi']*self.x[k,w] for k in self.process for w in self.periodsSeason),name="BigMFvin")
        self.modelL2.addConstr((self.Icake[self.periodsSeason[0]]==self.I0cake+gp.quicksum(self.A.at[("Filter Cake",k),self.periodsSeason[0]]*self.M3[k,self.periodsSeason[0]] for k in self.process)-self.Mcake[self.periodsSeason[0]]-self.M0cake),name="bdx1_0")
        self.modelL2.addConstrs((self.Icake[self.periodsSeason[w]]==self.Icake[self.periodsSeason[w-1]]+gp.quicksum(self.A.at[("Filter Cake",k),self.periodsSeason[w]]*self.M3[k,self.periodsSeason[w]] for k in self.process)-self.Mcake[self.periodsSeason[w]]for w in range(1,len(self.periodsSeason))),name="bdx1")
        self.modelL2.addConstrs((self.Icake[self.periodsN2[w]]==self.Icake[self.periodsN2[w-1]]-self.Mcake[self.periodsN2[w]]for w in range(self.season,self.numWeeks)),name="bdx1")
        #self.modelL2.addConstrs((self.Icake[self.periodsN2[w]]==self.Icake[self.periodperiodsN2sYear[w-1]]-self.Mcake[self.periodsN2[w]]for w in range(self.numWeeks,len(self.periodsYear))),name="bdx1_OffSeason")
        self.modelL2.addConstr((self.M0cake<=self.I0cake),name="bdx1_0")
        self.modelL2.addConstrs((self.Mcake[w]<=self.digmaxcake*self.dataHorizon.at[w,'nbio'] for k in self.process for w in self.periodsN2),name="bd2")
        self.modelL2.addConstrs((self.Mcake[self.periodsN2[w]]<=self.Mcake[self.periodsN2[w-1]]+self.deltacake*self.dataHorizon.at[self.periodsN2[w],'nbio'] for w in range(1,len(self.periodsN2))),name="bdCakeStab")
        self.modelL2.addConstr((self.Mcake[self.periodsN2[0]]<=self.M0cake+self.deltacake*self.dataHorizon.at[self.periodsN2[0],'nbio']),name="bdCakeStab0")
        self.modelL2.addConstr((self.M0cake<=self.Mcake[self.periodsN2[-1]]+self.deltacake*self.dataHorizon.at[self.periodsN2[0],'nbio']),name="bdCakeStabFinal")
        #self.modelL2.addConstrs((gp.quicksum(self.dataProcess.at[k,'OMvin']*(self.A.at[("Vinasse",k),self.periodsN2[w]]*self.M3[k,self.periodsN2[w]]-self.Fvin[k,self.periodsN2[w]])for k in self.process)>=
        #                         gp.quicksum(self.dataProcess.at[k,'OMvin']*(self.A.at[("Vinasse",k),self.periodsN2[w-1]]*self.M3[k,self.periodsN2[w-1]]-self.Fvin[k,self.periodsN2[w-1]]) for k in self.process)-self.deltaOM*self.dataHorizon.at[self.periodsN2[w],'n'] for w in range(1,len(self.periodsN2))),name="bd3")
        #self.modelL2.addConstr((gp.quicksum(self.dataProcess.at[k,'OMvin']*(self.A.at[("Vinasse",k),self.periodsN2[0]]*self.M3[k,self.periodsN2[0]]-self.Fvin[k,self.periodsN2[0]])for k in self.process)>=0),name="bd3_0")
        self.modelL2.addConstrs((gp.quicksum(self.dataProcess.at[k,'OMvin']*(self.A.at[("Vinasse",k),self.periodsSeason[w]]*self.M3[k,self.periodsSeason[w]]-self.Fvin[k,self.periodsSeason[w]])for k in self.process)<=
                                 gp.quicksum(self.dataProcess.at[k,'OMvin']*(self.A.at[("Vinasse",k),self.periodsSeason[w-1]]*self.M3[k,self.periodsSeason[w-1]]-self.Fvin[k,self.periodsSeason[w-1]]) for k in self.process)+self.deltaOMvin*self.dataHorizon.at[self.periodsSeason[w],'nbio'] for w in range(1,self.season)),name="bd4")
        self.modelL2.addConstr((gp.quicksum(self.dataProcess.at[k,'OMvin']*(self.A.at[("Vinasse",k),self.periodsSeason[0]]*self.M3[k,self.periodsSeason[0]]-self.Fvin[k,self.periodsSeason[0]])for k in self.process)<=self.deltaOMvin*self.dataHorizon.at[self.periodsSeason[0],'nbio']),name="bd4_0")
        self.modelL2.addConstr((self.Icake[self.periodsN2[-1]]>=self.Epbcake),name="bd5")
        self.modelL2.addConstr((self.PUbg[self.periodsSeason[0]]==self.BiogasCake*self.M0cake),name="bd6_0")
        self.modelL2.addConstrs((self.PUbg[self.periodsN2[w]]==gp.quicksum(self.A.at[("VinasseBiogas",k),self.periodsSeason[w-1]]*(self.A.at[("Vinasse",k),self.periodsSeason[w-1]]*self.M3[k,self.periodsSeason[w-1]]-self.Fvin[k,self.periodsSeason[w-1]]) for k in self.process)
                                                            + self.BiogasCake*(self.Mcake[self.periodsSeason[w-1]]) for w in range(1,self.season+1)),name="bd6Season")
        self.modelL2.addConstrs((self.PUbg[self.periodsN2[w]]==self.BiogasCake*(self.Mcake[self.periodsN2[w-1]]) for w in range(self.season+1,self.numWeeks)),name="bd6OffSeason")
        #self.modelL2.addConstrs((self.PUbg[self.periodsN2[w]]==gp.quicksum(self.A.at[("CakeBiogas",k),self.periodsN2[w-1]]*(self.Mcake[self.periodsN2[w]]) for k in self.process) for w in range(self.numWeeks,len(self.periodsYear))),name="bd6")
        self.modelL2.addConstrs((self.Mbg[w]<=self.tgmax*self.dataHorizon.at[w,'nbio'] for w in self.periodsN2),name="bd7")
        self.modelL2.addConstrs((self.PUbm[w]==self.txbiomet*(self.PUbg[w]-self.Mbg[w]-self.Q[w]) for w in self.periodsN2),name="bd8")
        self.modelL2.addConstrs((self.EE[w]==self.EG[w]-(gp.quicksum(self.cfe*self.M3[k,w] for k in self.process)+gp.quicksum(self.dataProdMill.at[p,'cve']*self.PU2[p,w] for p in self.ProductsMill)) for w in self.periodsSeason),name="const42Season")
        self.modelL2.addConstrs((self.EE[self.periodsN2[w]]==self.EG[self.periodsN2[w]] for w in range(self.season,self.numWeeks)),name="const42OffSeason")
        self.modelL2.addConstrs((self.EE[w]<=self.maxEE for w in self.periodsN2),name="maxEE")
        self.modelL2.addConstrs((self.EB[w] == self.rcftg*self.Mbg[w] for w in self.periodsN2),name="bigasEnergy")
        self.modelL2.addConstrs((self.PUbm[w]>=self.dataHorizon.at[w,"Dbmmin"]+self.dataHorizon.at[w,"DbmPipe"] for w in self.periodsN2),name="bdmin")
        self.modelL2.addConstrs((self.PUbm[w]<=self.dataHorizon.at[w,"Dbmmax"]+self.dataHorizon.at[w,"DbmPipe"] for w in self.periodsN2),name="bdmax")
        #self.modelL2.addConstrs((self.Fcake[k,w]<=self.G*self.x[k,w] for k in self.process for w in self.periodsN2),name="BigMFcake")
        self.modelL2.setObjective((gp.quicksum(self.PU2[pa,w]/self.P2*self.R2*modelL1.dataProd.at['ConvCrpaaa',pa] for w in self.periodsSeason for pa in self.SugarMill)
                            +gp.quicksum(self.PU2[pe,w]/self.P3*self.R3*modelL1.dataProd.at['ConvCrpaaa',pe] for w in self.periodsSeason for pe in self.EthanolMill)
                            +gp.quicksum(self.PU2[pm,w]/self.P4*self.R4*modelL1.dataMills.at['Matr',self.millName] for w in self.periodsSeason for pm in self.MolassesMill)
                            +gp.quicksum(self.dataHorizon.at[w,"VbmPipe"]*self.dataHorizon.at[w,"DbmPipe"] for w in self.periodsN2)
                            +gp.quicksum(self.dataHorizon.at[w,"VbmDiesel"]*(self.PUbm[w]-self.dataHorizon.at[w,"DbmPipe"]) for w in self.periodsN2)
                            +gp.quicksum(self.ve*self.EE[w] for w in self.periodsN2)
                            +gp.quicksum(self.dataHorizon.at[w,"Vebio"]*self.EB[w] for w in self.periodsN2)
                            -gp.quicksum(self.Cm.at[m,w]*self.M1[m,w] for m in self.sourceRaw for w in self.periodsSeason)
                            -gp.quicksum(self.Cf.at[f,w]*self.M2[f,w] for f in self.transport for w in self.periodsSeason)
                            -gp.quicksum(self.Ck.at[k,w]*self.M3[k,w] for k in self.process for w in self.periodsSeason)
                            -gp.quicksum(modelL1.ca*self.ImMill[p,w] for p in self.ProductsMill for w in self.periodsSeason)
                            -gp.quicksum(self.biogasPlantCost*self.PUbg[w] for w in self.periodsSeason)
                            #-gp.quicksum(self.Fvin[k,w] for k in self.process for w in self.periodsSeason)
                            -1e-4*self.maxEE
                            ),gp.GRB.MAXIMIZE)
        self.modelL2.update()
        #self.modelL2.write("MPBRAPL2.lp")
        return self.modelL2
