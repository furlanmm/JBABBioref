import pandas as pd
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.ticker as tk
import gurobipy_pandas as gppd
import gurobipy as gp
import itertools as it

class PASUCL1:
    #Data from Cooperative
    general: pd.DataFrame = ...
    mills: list = ...
    periods: list = ...
    products: list = ...
    sugar: list = ...
    ethanol: list = ...
    molasses: list = ...
    inventory: list = ...  
    phi: float = ...
    ca: int = ...
    n: pd.DataFrame = ...
    eatr: pd.DataFrame = ...
    atr: pd.DataFrame = ...
    dmin: pd.DataFrame = ...
    vp: pd.DataFrame = ...
    dataMills: pd.DataFrame = ...
    theta: pd.DataFrame = ...
    cpmax: pd.DataFrame = ...
    cp: pd.DataFrame = ...
    inv: pd.DataFrame = ...
    emax: pd.DataFrame = ...
    ce: pd.DataFrame = ...
    i0: pd.DataFrame = ...
    dataProd: pd.DataFrame = ...
    #PASUCL1 Variables
    M: gp.Var = ...
    Ip: gp.Var = ...
    Im: gp.Var = ...
    PU: gp.Var = ...
    D: gp.Var = ...
    #Models
    modelL1: gp.Model = ...


    def readDataCoop(self,fileName):
        #This function read the data from a xlsx file which containts the information from the whole group of mills acording to model PASUC/L1
        #PAIVA, Rafael P.; ROCCO, Cleber D.; MORABITO, Reinaldo. A hierarchical supply chain model for the sugar–alcohol energy sector with robust optimization analysis. International Transactions in Operational Research, v. 30, n. 4, p. 1789-1818, 2023.
        self.general = pd.read_excel(fileName, sheet_name = "General",index_col=0, header=None)
        self.mills = self.general.filter(like="Mills", axis = 0).dropna(axis=1).values.tolist()[0]
        self.periods = self.general.filter(like="Periods", axis = 0).dropna(axis=1).values.tolist()[0]
        self.products = self.general.filter(like="Products", axis = 0).dropna(axis=1).values.tolist()[0]
        self.sugar = self.general.filter(like="Sugar", axis = 0).dropna(axis=1).values.tolist()[0]
        self.ethanol = self.general.filter(like="Ethanol", axis = 0).dropna(axis=1).values.tolist()[0]
        self.molasses = self.general.filter(like="Molasses", axis = 0).dropna(axis=1).values.tolist()[0]
        self.inventory = self.general.filter(like="Inventory", axis = 0).dropna(axis=1).values.tolist()[0]  
        self.phi = self.general.filter(like="phi", axis = 0).dropna(axis=1).values.tolist()[0][0]
        self.ca = self.general.filter(like="Ca", axis = 0).dropna(axis=1).values.tolist()[0][0]
        self.n = pd.read_excel(fileName, sheet_name = "n",index_col=0)#.replace(np.nan,0)
        self.eatr = pd.read_excel(fileName, sheet_name = "Eatr",index_col=0)#.replace(np.nan,0)
        self.atr = pd.read_excel(fileName, sheet_name = "Atr",index_col=0)#.replace(np.nan,0)
        self.dmin = pd.read_excel(fileName, sheet_name = "Dmin",index_col=0).replace(np.nan,0)
        self.vp = pd.read_excel(fileName, sheet_name = "Vp",index_col=0).replace(np.nan,0)
        self.dataMills = pd.read_excel(fileName, sheet_name = "Mills",index_col=0)
        self.theta = pd.read_excel(fileName, sheet_name = "Theta",index_col=0)#.replace(np.nan,0)
        self.cpmax = pd.read_excel(fileName, sheet_name = "Cpmax",index_col=0).replace(np.nan,-1)
        self.cp = pd.read_excel(fileName, sheet_name = "Cp",index_col=0).replace(np.nan,0)
        self.inv = pd.read_excel(fileName, sheet_name = "Inventory",index_col=0).replace(np.nan,0)
        self.emax = pd.read_excel(fileName, sheet_name = "Emax",index_col=0).replace(np.nan,-1)
        self.ce = pd.read_excel(fileName, sheet_name = "Ce",index_col=0).replace(np.nan,99999)
        self.i0 = pd.read_excel(fileName, sheet_name = "I0",index_col=0).replace(np.nan,0)
        self.dataProd = pd.read_excel(fileName, sheet_name = "Products",index_col=0)

    def createModelL1(self,name):
        #This function create the PASUC/L1 model
        self.modelL1 = gp.Model(name)
        #Variables Instances
        self.M = self.modelL1.addVars(self.mills,self.periods,name="Mut")
        self.Ip = self.modelL1.addVars(self.products,self.inventory,self.periods,name="I+")
        self.Im = self.modelL1.addVars(self.products,self.periods,name="I-")
        self.PU = self.modelL1.addVars(self.products,self.mills,self.periods,name="PUput")
        self.D = self.modelL1.addVars(self.products,self.periods,name="Dpt")
        self.modelL1.update()
        self.modelL1.addConstrs((gp.quicksum(self.Ip[p,e,t] for e in self.inventory if self.emax.at[p,e] != -1) == 
                       gp.quicksum(self.Ip[p,e,self.periods[self.periods.index(t)-1]] for e in self.inventory if self.emax.at[p,e] != -1)+
                       gp.quicksum(self.PU[p,u,t] for u in self.mills if self.cpmax.at[p,u] != -1)-self.D[p,t]+
                       self.Im[p,t]-self.Im[p,self.periods[self.periods.index(t)-1]] for p in self.products for t in self.periods
                       if t != self.periods[0]), name = "const2")
        self.modelL1.addConstrs((gp.quicksum(self.Ip[p,e,self.periods[0]] for e in self.inventory if self.emax.at[p,e] != -1) == 
                       self.i0.filter(like=p,axis=0).replace(np.nan,0).values.sum() + gp.quicksum(self.PU[p,u,self.periods[0]] for u in self.mills if self.cpmax.at[p,u] != -1)-
                       self.D[p,self.periods[0]]+ self.Im[p,self.periods[0]]-self.inv.at[p,"Ainic"] for p in self.products), name = "const2_0")
        self.modelL1.addConstrs((self.D[p,t] >= self.dmin.at[p,t] for p in self.products for t in self.periods), name="const3")
        self.modelL1.addConstrs((gp.quicksum(self.D[p,t] for t in self.periods) <= self.dataProd.at["Dmax",p] for p in self.products),name="const4")
        self.modelL1.addConstrs((gp.quicksum(self.M[u,t] for t in self.periods) <= self.dataMills.at["C0",u] for u in self.mills),name="const5")
        self.modelL1.addConstrs((self.M[u,t] >= self.dataMills.at["Mmin",u]*self.n.at[u,t]*self.phi  for u in self.mills for t in self.periods),name="constr6a")
        self.modelL1.addConstrs((self.M[u,t] <= self.dataMills.at["Mmax",u]*self.n.at[u,t]*self.phi for u in self.mills for t in self.periods),name="constr6b")
        self.modelL1.addConstrs((self.Ip[p,e,t] <= self.emax.at[p,e] for p in self.products for e in self.inventory for t in self.periods if self.emax.at[p,e] != -1),name="constr7")
        self.modelL1.addConstrs((gp.quicksum(self.Ip[p,e,self.periods[-1]] for e in self.inventory if self.emax.at[p,e] != -1)>=self.inv.at[p,"Ep"] for p in self.products),name="constr8a")
        self.modelL1.addConstrs((self.Im[p,self.periods[-1]] == 0 for p in self.products),name="constr8b")
        self.modelL1.addConstrs((self.PU[p,u,t]<=self.cpmax.at[p,u]*self.n.at[u,t] for p in self.products for u in self.mills for t in self.periods if self.cpmax.at[p,u] != -1),name="constr9")
        self.modelL1.addConstrs((gp.quicksum(self.PU[pa,u,t] for pa in self.sugar if self.cpmax.at[pa,u] != -1) <= self.dataMills.at["Cfmax",u]*self.n.at[u,t] for u in self.mills for t in self.periods),name="const10")
        self.modelL1.addConstrs((gp.quicksum(self.PU[pe,u,t] for pe in self.ethanol if self.cpmax.at[pe,u] != -1) <= self.dataMills.at["Cdmax",u]*self.n.at[u,t] for u in self.mills for t in self.periods),name="const11")
        self.modelL1.addConstrs(((self.M[u,t]*self.atr.at[u,t]*self.eatr.at[u,t])/1000-gp.quicksum(self.PU[p,u,t]*self.dataProd.at["Patr",p] for p in self.products if self.cpmax.at[p,u] != -1) == 0 for u in self.mills for t in self.periods),name="const12")
        self.modelL1.addConstrs((gp.quicksum(self.PU[pe,u,t]*self.dataProd.at["Patr",pe] for pe in self.ethanol if self.cpmax.at[pe,u] != -1)+
                        gp.quicksum(self.PU[pm,u,t]*self.dataMills.at["Matr",u] for pm in self.molasses if self.cpmax.at[pm,u] != -1)-
                        gp.quicksum(self.PU[pa,u,t]*self.theta.at[pa,u]*self.dataMills.at["Matr",u] for pa in self.sugar if self.cpmax.at[pa,u] != -1) >= 0
                        for u in self.mills for t in self.periods),name="const13")

        self.modelL1.setObjective((gp.quicksum(self.vp.at[p,t]*self.D[p,t] for p in self.products for t in self.periods)
                        -gp.quicksum(self.cp.at[p,u]*self.PU[p,u,t] for p in self.products for u in self.mills for t in self.periods)
                        -gp.quicksum(self.ce.at[p,e]*self.Ip[p,e,t] for p in self.products for e in self.inventory for t in self.periods)
                        -gp.quicksum(self.ca*self.Im[p,t] for p in self.products for t in self.periods)
                        ),gp.GRB.MAXIMIZE)
        self.modelL1.update()

        return self.modelL1
    
class PASUCL2:
    #Internal Data between models
    Goals: pd.DataFrame = ...
    SjM: pd.DataFrame = ...
    SmF:float = ...
    Rd: pd.DataFrame = ...
    Ra: pd.DataFrame = ...
    Re: pd.DataFrame = ...
    Rm: pd.DataFrame = ...
    A: pd.DataFrame = ...
    ConvA: pd.DataFrame = ...
    ConvD: pd.DataFrame = ...
    ConvM: pd.DataFrame = ...
    ConvAll: pd.DataFrame = ...
    ARTproc: pd.DataFrame = ...
    CustoEnsaque: pd.DataFrame = ...
    Ck: pd.DataFrame = ...
    Cm: pd.DataFrame = ...
    P1:float = ...
    P2:float = ...
    P3:float = ...
    P4:float = ...
    R1:float = ...
    R2:float = ...
    R3:float = ...
    R4:float = ...
    G: int = 100000
    #Data from mill
    generalN2: pd.DataFrame = ...
    millName: str = ... 
    numProc: int = ...
    numWeeks: int = ...
    process: list = ...
    periodsN2: list = ...
    ProductsMill: list = ...
    SugarMill: list = ...
    EthanolMill: list = ...
    MolassesMill: list = ...
    sourceRaw: list = ...
    selfRaw: list = ...
    externalRaw: list = ...
    transport: list = ...
    periodsGoals: pd.DataFrame = ...
    dataHorizon: pd.DataFrame = ...
    dataProcess: pd.DataFrame = ...
    dataSourceRaw: pd.DataFrame = ...
    Gf: pd.DataFrame = ...
    Gm: pd.DataFrame = ...
    Gk: pd.DataFrame = ...
    dataTransport: pd.DataFrame = ...
    dataProdMill: pd.DataFrame = ...
    atrMill: pd.DataFrame = ...
    fiber: pd.DataFrame = ...
    Cf: pd.DataFrame = ...
    beta: pd.DataFrame = ...
    PzaM:float = ...
    Bmf:float = ...
    rcb:float = ...
    Gfix:float = ...
    initBudget:float = ...
    eb:float = ...
    Epb:float = ...
    rc:float = ...
    cfvap:float = ...
    rcf:float = ...
    cfe:float = ...
    vapmax:float = ...
    egmax:float = ...
    arm:float = ...
    cproc:float = ...
    ve:float = ...
    #PASUC2 Variables
    x: gp.Var = ...
    MMill: gp.Var = ...
    M1: gp.Var = ...
    M2: gp.Var = ...
    M3: gp.Var = ...
    PU2: gp.Var = ...
    C: gp.Var = ...
    S: gp.Var = ...
    ImMill: gp.Var = ...
    Ib: gp.Var = ...
    Mb: gp.Var = ...
    VAP: gp.Var = ...
    EG: gp.Var = ...
    EE: gp.Var = ...
    #Models
    modelL2: gp.Model = ...

    def readDataMill(self, fileName):
        #This function read the data from a xlsx file which containts the information from the whole group of mills acording to model PASUC/L2
        #PAIVA, Rafael P.; ROCCO, Cleber D.; MORABITO, Reinaldo. A hierarchical supply chain model for the sugar–alcohol energy sector with robust optimization analysis. International Transactions in Operational Research, v. 30, n. 4, p. 1789-1818, 2023.
        self.generalN2 = pd.read_excel(fileName, sheet_name = "General",index_col=0, header=None)
        self.millName = self.generalN2.filter(like="Mill", axis = 0).dropna(axis=1).values.tolist()[0][0] 
        self.numProc = self.generalN2.filter(like="numProc", axis = 0).dropna(axis=1).values.tolist()[0][0]
        self.numWeeks = self.generalN2.filter(like="Periods", axis = 0).dropna(axis=1).values.tolist()[0][0]
        self.process = ['Proc'+str(i) for i in range(1,self.numProc+1)]
        self.periodsN2 = ['Week'+str(t) for t in range(1,self.numWeeks+1)]
        self.ProductsMill = self.generalN2.filter(like="ProductsMill", axis = 0).dropna(axis=1).values.tolist()[0]
        self.SugarMill = self.generalN2.filter(like="SugarMill", axis = 0).dropna(axis=1).values.tolist()[0]
        self.EthanolMill = self.generalN2.filter(like="EthanolMill", axis = 0).dropna(axis=1).values.tolist()[0]
        self.MolassesMill = self.generalN2.filter(like="MolassesMill", axis = 0).dropna(axis=1).values.tolist()[0]
        self.sourceRaw = self.generalN2.filter(like="sourceRaw", axis = 0).dropna(axis=1).values.tolist()[0]
        self.selfRaw = self.generalN2.filter(like="SelfRaw", axis = 0).dropna(axis=1).values.tolist()[0]
        self.externalRaw = self.generalN2.filter(like="ExternalRaw", axis = 0).dropna(axis=1).values.tolist()[0]
        self.transport = self.generalN2.filter(like="Transport", axis = 0).dropna(axis=1).values.tolist()[0]
        self.periodsGoals = pd.read_excel(fileName, sheet_name = "PeriodsGoals",index_col=0, header=0)
        self.dataHorizon = pd.read_excel(fileName, sheet_name = "Horizon",index_col=0, header=0)
        self.dataProcess = pd.read_excel(fileName, sheet_name = "Process",index_col=0, header=0)
        self.dataSourceRaw = pd.read_excel(fileName, sheet_name = "sourceRaw",index_col=0, header=0)
        self.Gf = pd.read_excel(fileName, sheet_name = "Gf",index_col=0, header=0)
        self.Gm = pd.read_excel(fileName, sheet_name = "Gm",index_col=0, header=0)
        self.Gk = pd.read_excel(fileName, sheet_name = "Gk",index_col=0, header=0)
        self.dataTransport = pd.read_excel(fileName, sheet_name = "Transport",index_col=0, header=0)
        self.dataProdMill = pd.read_excel(fileName, sheet_name = "Products",index_col=0, header=0)
        self.atrMill = pd.read_excel(fileName, sheet_name = "atr",index_col=0, header=0)
        self.fiber = pd.read_excel(fileName, sheet_name = "Fiber",index_col=0, header=0)
        self.Cf = pd.read_excel(fileName, sheet_name = "Cf",index_col=0, header=0)
        self.beta = pd.read_excel(fileName, sheet_name = "beta",index_col=0, header=0)
        self.dataHorizon['PzaJ'] = self.dataHorizon['PzaC'] - 1
        self.dataProcess['PzaS'] = self.dataProcess['Pol']/(1-(self.dataProcess['Umid']/100))
        self.PzaM = self.generalN2.filter(like="PzaM", axis = 0).dropna(axis=1).values.tolist()[0][0] 
        self.Bmf = self.generalN2.filter(like="Bmf", axis = 0).dropna(axis=1).values.tolist()[0][0] 
        self.rcb = self.generalN2.filter(like="rcb", axis = 0).dropna(axis=1).values.tolist()[0][0]
        #self.Gfix = self.generalN2.filter(like="Gfix", axis = 0).dropna(axis=1).values.tolist()[0][0]
        self.initBudget = self.generalN2.filter(like="Budget", axis = 0).dropna(axis=1).values.tolist()[0][0]
        self.eb = self.generalN2.filter(like="eb", axis = 0).dropna(axis=1).values.tolist()[0][0]
        self.Epb = self.generalN2.filter(like="Epb", axis = 0).dropna(axis=1).values.tolist()[0][0]
        self.rc = self.generalN2.filter(regex="rc$", axis = 0).dropna(axis=1).values.tolist()[0][0]
        self.cfvap = self.generalN2.filter(like="cfvap", axis = 0).dropna(axis=1).values.tolist()[0][0]
        self.rcf = self.generalN2.filter(like="rcf", axis = 0).dropna(axis=1).values.tolist()[0][0]
        self.cfe = self.generalN2.filter(like="cfe", axis = 0).dropna(axis=1).values.tolist()[0][0]
        self.vapmax = self.generalN2.filter(like="vapmax", axis = 0).dropna(axis=1).values.tolist()[0][0]
        self.egmax = self.generalN2.filter(like="egmax", axis = 0).dropna(axis=1).values.tolist()[0][0]
        self.arm = self.generalN2.filter(like="arm", axis = 0).dropna(axis=1).values.tolist()[0][0]
        self.cproc = self.generalN2.filter(like="cproc", axis = 0).dropna(axis=1).values.tolist()[0][0]
        self.ve = self.generalN2.filter(like="Ve", axis = 0).dropna(axis=1).values.tolist()[0][0]

    def calcParamsL2(self,modelL1):
        #This function uses data from Mill to calculate part of the parameters from PASUC/L2 
        #Begin of productive matrix values calculation (A_pkw, p:products, k:processes, w:weeks)
        arraySjM = np.zeros((len(self.process),len(self.periodsN2)),np.float64)
        for k in self.process:
            for w in self.periodsN2:
                arraySjM[self.process.index(k)][self.periodsN2.index(w)] = self.dataProcess.at[k,'PzaS']*(self.dataHorizon.at[w,'PzaJ']-self.PzaM)/ (self.dataHorizon.at[w,'PzaJ']*(self.dataProcess.at[k,'PzaS']-self.PzaM))
        self.SjM = pd.DataFrame(arraySjM,index=self.process,columns=self.periodsN2)  
        self.SmF = (self.Bmf*self.PzaM)/100
        arrayRd = np.zeros((len(self.EthanolMill),len(self.periodsN2)),np.float64)
        for pe in self.EthanolMill:
            for w in self.periodsN2:
                arrayRd[self.EthanolMill.index(pe)][self.periodsN2.index(w)] = (self.dataProdMill.at[pe,'Rend']*self.dataHorizon.at[w,'Efd'])/100  
        self.Rd = pd.DataFrame(arrayRd,index=self.EthanolMill,columns=self.periodsN2)
        arrayRa = np.zeros((len(self.SugarMill),len(self.process),len(self.periodsN2)),np.float64)
        arrayRe = np.zeros((len(self.EthanolMill),len(self.process),len(self.periodsN2)),np.float64)
        arrayRm = np.zeros((len(self.MolassesMill),len(self.process),len(self.periodsN2)),np.float64)
        #print(dataHorizon,dataProcess)
        for k in self.process:
            for w in self.periodsN2:
                for pa in self.SugarMill:
                    arrayRa[self.SugarMill.index(pa)][self.process.index(k)][self.periodsN2.index(w)] = self.dataHorizon.at[w,'pc']*10*(self.dataHorizon.at[w,'Elbti']/100)*self.SjM.at[k,w]*self.dataProcess.at[k,'Ts']*self.dataProcess.at[k,pa]
                for pe in self.EthanolMill:
                    arrayRe[self.EthanolMill.index(pe)][self.process.index(k)][self.periodsN2.index(w)] = (self.dataHorizon.at[w,'Elbti']/10)*(self.Rd.at[pe,w]/100)*(((self.dataHorizon.at[w,'pc']/0.95)*(1-self.SjM.at[k,w])+self.dataHorizon.at[w,'Ar'])*self.dataProcess.at[k,'Ts']*self.dataProcess.at[k,'Tm']+
                                                                                                                                            (self.dataHorizon.at[w,'pc']/0.95+self.dataHorizon.at[w,'Ar'])*(1-self.dataProcess.at[k,'Ts']))*self.dataProcess.at[k,pe]                                                                                                            
                for pm in self.MolassesMill:
                    arrayRm[self.MolassesMill.index(pm)][self.process.index(k)][self.periodsN2.index(w)] = self.dataHorizon.at[w,'pc']*self.dataHorizon.at[w,'Elbti']*(1-self.SjM.at[k,w])*self.dataProcess.at[k,'Ts']*100*(1-self.dataProcess.at[k,'Tm'])/(10*self.SmF)
            self.Ra = pd.DataFrame([[arrayRa[self.SugarMill.index(pa)][self.process.index(k)][self.periodsN2.index(w)] for w in self.periodsN2] for pa in self.SugarMill for k in self.process],index=it.product(self.SugarMill,self.process),columns=self.periodsN2)
            self.Re = pd.DataFrame([[arrayRe[self.EthanolMill.index(pe)][self.process.index(k)][self.periodsN2.index(w)] for w in self.periodsN2] for pe in self.EthanolMill for k in self.process],index=it.product(self.EthanolMill,self.process),columns=self.periodsN2)
            self.Rm = pd.DataFrame([[arrayRm[self.MolassesMill.index(pm)][self.process.index(k)][self.periodsN2.index(w)] for w in self.periodsN2] for pm in self.MolassesMill for k in self.process],index=it.product(self.MolassesMill,self.process),columns=self.periodsN2)
            self.A = pd.concat([self.Ra, self.Re, self.Rm]).div(1000)
        #End of productive matrix values
        #Begin of industrial cost of TRS produced in process k in week w  calculation(Ck)
        arrayConvA = np.zeros((len(self.SugarMill),len(self.process),len(self.periodsN2)),np.float64)
        arrayConvD = np.zeros((len(self.EthanolMill),len(self.process),len(self.periodsN2)),np.float64)
        arrayConvM = np.zeros((len(self.MolassesMill),len(self.process),len(self.periodsN2)),np.float64)
        for k in self.process:
            for w in self.periodsN2:
                for pa in self.SugarMill:
                    arrayConvA[self.SugarMill.index(pa)][self.process.index(k)][self.periodsN2.index(w)] = self.Ra.at[(pa,k),w]*(self.dataProcess.at[k,'Pol']/100)*(1-self.dataProcess.at[k,'Umid']/100)*self.dataProcess.at[k,pa]
                for pe in self.EthanolMill:
                    arrayConvD[self.EthanolMill.index(pe)][self.process.index(k)][self.periodsN2.index(w)] = self.Re.at[(pe,k),w]*self.dataProdMill.at[pe,'Fator']*self.dataProcess.at[k,pe]
                for pm in self.MolassesMill:
                    arrayConvM[self.MolassesMill.index(pm)][self.process.index(k)][self.periodsN2.index(w)] = (self.Rm.at[(pm,k),w]*self.SmF/100)/0.95+self.Rm.at[(pm,k),w]*self.arm/100
        self.ConvA = pd.DataFrame([[arrayConvA[self.SugarMill.index(pa)][self.process.index(k)][self.periodsN2.index(w)] for w in self.periodsN2] for pa in self.SugarMill for k in self.process],index=it.product(self.SugarMill,self.process),columns=self.periodsN2)
        self.ConvD = pd.DataFrame([[arrayConvD[self.EthanolMill.index(pe)][self.process.index(k)][self.periodsN2.index(w)] for w in self.periodsN2] for pe in self.EthanolMill for k in self.process],index=it.product(self.EthanolMill,self.process),columns=self.periodsN2)
        self.ConvM = pd.DataFrame([[arrayConvM[self.MolassesMill.index(pm)][self.process.index(k)][self.periodsN2.index(w)] for w in self.periodsN2] for pm in self.MolassesMill for k in self.process],index=it.product(self.MolassesMill,self.process),columns=self.periodsN2)
        self.ConvAll = pd.concat([self.ConvA, self.ConvD, self.ConvM])
        arrayARTproc = np.zeros((len(self.process),len(self.periodsN2)),np.float64)
        arrayCustoEnsaque = np.zeros((len(self.process),len(self.periodsN2)),np.float64)
        for k in self.process:
            for w in self.periodsN2:
                arrayARTproc[self.process.index(k)][self.periodsN2.index(w)] = sum(self.ConvAll.at[(p,k),w]/self.dataProdMill.at[p,'divART'] for p in self.ProductsMill)
                arrayCustoEnsaque[self.process.index(k)][self.periodsN2.index(w)] = sum((self.Ra.at[(pa,k),w]*self.dataProcess.at[k,'CustoEnsaque'])/50.0 for pa in self.SugarMill)
        self.ARTproc = pd.DataFrame(arrayARTproc,index=self.process,columns=self.periodsN2)
        self.CustoEnsaque = pd.DataFrame(arrayCustoEnsaque,index=self.process,columns=self.periodsN2)
        self.Ck = self.ARTproc.mul(self.cproc).add(self.CustoEnsaque)
        #Begin of cost of TRS from supplier m in week w  calculation(Cm)
        arrayCm = np.zeros((len(self.sourceRaw),len(self.periodsN2)),np.float64)
        for m in self.sourceRaw:
            for w in self.periodsN2:
                arrayCm[self.sourceRaw.index(m)][self.periodsN2.index(w)] = self.atrMill.at[m,w]*self.dataHorizon.at[w,'vAtr']+self.dataSourceRaw.at[m,'delta']
        self.Cm = pd.DataFrame(arrayCm,index=self.sourceRaw,columns=self.periodsN2)
        #Input calculation for PASUC/L2 from outputs of PASUC/L1
        self.P1 = sum(modelL1.PU[p,u,t].X*modelL1.dataProd.at["Patr",p] for p in modelL1.products for u in modelL1.mills for t in modelL1.periods)
        self.P2 = sum(modelL1.PU[pa,u,t].X*modelL1.dataProd.at["ConvCrpaaa",pa] for pa in modelL1.sugar for u in modelL1.mills for t in modelL1.periods)
        self.P3 = sum(modelL1.PU[pd,u,t].X*modelL1.dataProd.at["ConvCrpaaa",pd] for pd in modelL1.ethanol for u in modelL1.mills for t in modelL1.periods)
        self.P4 = sum(modelL1.PU[pm,u,t].X*modelL1.dataMills.at["Matr",u] for pm in modelL1.molasses for u in modelL1.mills for t in modelL1.periods)
        self.R1 = sum(modelL1.vp.at[p,t]*modelL1.D[p,t].X for p in modelL1.products for t in modelL1.periods)
        self.R2 = sum(modelL1.vp.at[pa,t]*modelL1.D[pa,t].X for pa in modelL1.sugar for t in modelL1.periods)
        self.R3 = sum(modelL1.vp.at[pd,t]*modelL1.D[pd,t].X for pd in modelL1.ethanol for t in modelL1.periods)
        self.R4 = 0 #sum(modelL1.vp.at[pm,t]*modelL1.D[pm,t].X for pm in modelL1.molasses for t in modelL1.periods)
        arrayGoals = np.zeros((len(self.ProductsMill),len(self.periodsN2)),dtype=np.float64)
        #Goals for mill self.millName
        for i in self.periodsGoals.index:
            for p in self.ProductsMill:
                arrayGoals[self.ProductsMill.index(p)][self.periodsN2.index(self.periodsGoals.at[i,"Week"])] = modelL1.PU[p,self.millName,self.periodsGoals.at[i,"Month"]].X
        self.Goals = pd.DataFrame(arrayGoals,index=self.ProductsMill,columns=self.periodsN2)

    def createModelL2(self,name,modelL1):
        #This function create the PASUC/L2 model
        self.modelL2 = gp.Model(name)
        #Variables instances
        self.x = self.modelL2.addVars(self.process,self.periodsN2,name="Xkt",vtype=gp.GRB.BINARY)
        self.MMill = self.modelL2.addVars(self.periodsN2,name="Mt")
        self.M1 = self.modelL2.addVars(self.sourceRaw,self.periodsN2,name="M1mt")
        self.M2 = self.modelL2.addVars(self.transport,self.periodsN2,name="M2ft")
        self.M3 = self.modelL2.addVars(self.process,self.periodsN2,name="M3kt")
        self.PU2 = self.modelL2.addVars(self.ProductsMill,self.periodsN2,name="PU2")
        self.C = self.modelL2.addVars(self.sourceRaw,self.periodsN2,name="Cmt")
        self.S = self.modelL2.addVars(self.periodsN2,name="St")
        self.ImMill = self.modelL2.addVars(self.ProductsMill,self.periodsN2,name="I-pt")
        self.Ib = self.modelL2.addVars(self.periodsN2,name="Ibt")
        self.Mb = self.modelL2.addVars(self.periodsN2,name="Mbt")
        self.VAP = self.modelL2.addVars(self.periodsN2,name="VAPt")
        self.EG = self.modelL2.addVars(self.periodsN2,name="EGt")
        self.EE = self.modelL2.addVars(self.periodsN2,name="EEt")
        self.modelL2.update()
        self.calcParamsL2(modelL1)
        self.modelL2.addConstrs((gp.quicksum(self.A.at[(p,k),w]*self.M3[k,w] for k in self.process) == self.PU2[p,w]  for w in self.periodsN2 for p in self.ProductsMill),name="constPU")
        self.modelL2.addConstrs((gp.quicksum(self.PU2[p,wl] for wl in self.periodsN2 if self.periodsN2.index(wl)<=self.periodsN2.index(w))+self.ImMill[p,w]
                            -gp.quicksum(self.Goals.at[p,wl] for wl in self.periodsN2 if self.periodsN2.index(wl)<=self.periodsN2.index(w))>=0
                            for p in self.ProductsMill for w in self.periodsN2),name="const23")
        self.modelL2.addConstrs((gp.quicksum(self.M1[m,w] for m in self.sourceRaw) == self.MMill[w] for w in self.periodsN2),name="const24a")
        self.modelL2.addConstrs((gp.quicksum(self.M2[f,w] for f in self.transport) == self.MMill[w] for w in self.periodsN2),name="const24b")
        self.modelL2.addConstrs((gp.quicksum(self.M3[k,w] for k in self.process) == self.MMill[w] for w in self.periodsN2),name="const24c")
        self.modelL2.addConstrs((self.C[m,self.periodsN2[w]]==self.C[m,self.periodsN2[w-1]]-self.M1[m,self.periodsN2[w-1]] for m in self.sourceRaw for w in range(1,len(self.periodsN2))),name="const25")
        self.modelL2.addConstrs((self.C[m,self.periodsN2[0]]== self.dataSourceRaw.at[m,'C0'] for m in self.sourceRaw),name="const25_0")
        self.modelL2.addConstrs((self.C[m,w]>=self.M1[m,w] for m in self.sourceRaw for w in self.periodsN2),name="const25b")
        self.modelL2.addConstr((gp.quicksum(self.dataSourceRaw.at[m,'C0'] for m in self.sourceRaw)==gp.quicksum(self.M1[m,w] for m in self.sourceRaw for w in self.periodsN2)),name="const26")
        self.modelL2.addConstrs((gp.quicksum(self.C[mp,w] for mp in self.selfRaw)<=self.dataHorizon.at[w,'disp'] for w in self.periodsN2),name="const27")
        self.modelL2.addConstrs((gp.quicksum(self.M1[mf,w] for mf in self.externalRaw)<= self.dataHorizon.at[w,'alfa']*self.MMill[w] for w in self.periodsN2),name="const28")
        self.modelL2.addConstrs((self.MMill[w] >= modelL1.dataMills.at["Mmin",self.millName]*self.dataHorizon.at[w,'n']*modelL1.phi for w in self.periodsN2),name="constr29a")
        self.modelL2.addConstrs((self.MMill[w] <= modelL1.dataMills.at["Mmax",self.millName]*self.dataHorizon.at[w,'n']*modelL1.phi for w in self.periodsN2),name="constr29b")
        self.modelL2.addConstrs((self.S[self.periodsN2[w]]==self.S[self.periodsN2[w-1]]+(gp.quicksum((self.dataProdMill.at[p,'rca']+self.dataProdMill.at[p,'Patr']*self.rcb)*self.PU2[p,self.periodsN2[w]] for p in self.ProductsMill)+self.dataHorizon.at[self.periodsN2[w],'rcc'])
                        -(gp.quicksum(self.Gm.at[m,self.periodsN2[w]]*self.M1[m,self.periodsN2[w]] for m in self.sourceRaw)+gp.quicksum(self.Gf.at[f,self.periodsN2[w]]*self.M2[f,self.periodsN2[w]] for f in self.transport)
                        +gp.quicksum(self.Gk.at[k,self.periodsN2[w]]*self.M3[k,self.periodsN2[w]] for k in self.process)+self.dataHorizon.at[self.periodsN2[w],"Gfix"]) for w in range(1,len(self.periodsN2))),name="const30")
        self.modelL2.addConstr((self.S[self.periodsN2[0]] == self.initBudget+(gp.quicksum((self.dataProdMill.at[p,'rca']+self.dataProdMill.at[p,'Patr']*self.rcb)*self.PU2[p,self.periodsN2[0]] for p in self.ProductsMill)+self.dataHorizon.at[self.periodsN2[0],'rcc'])
                        -(gp.quicksum(self.Gm.at[m,self.periodsN2[0]]*self.M1[m,self.periodsN2[0]] for m in self.sourceRaw)+gp.quicksum(self.Gf.at[f,self.periodsN2[0]]*self.M2[f,self.periodsN2[0]] for f in self.transport)
                        +gp.quicksum(self.Gk.at[k,self.periodsN2[0]]*self.M3[k,self.periodsN2[0]] for k in self.process)+self.dataHorizon.at[self.periodsN2[0],"Gfix"])),name="const30_0")
        self.modelL2.addConstrs((self.M2[f,w]<=self.dataTransport.at[f,'fmax']*self.dataHorizon.at[w,'n']*self.beta.at[f,w] for f in self.transport for w in self.periodsN2),name="const31")
        self.modelL2.addConstrs((self.PU2[p,w]<=modelL1.cpmax.at[p,self.millName]*self.dataHorizon.at[w,'n'] for p in self.ProductsMill for w in self.periodsN2),name="const32")
        self.modelL2.addConstrs((gp.quicksum(self.PU2[pa,w] for pa in self.SugarMill) <= modelL1.dataMills.at['Cfmax',self.millName]*self.dataHorizon.at[w,'n'] for w in self.periodsN2),name="const33")
        self.modelL2.addConstrs((gp.quicksum(self.PU2[pe,w] for pe in self.EthanolMill) <= modelL1.dataMills.at['Cdmax',self.millName]*self.dataHorizon.at[w,'n'] for w in self.periodsN2),name="const34")
        self.modelL2.addConstrs((gp.quicksum(self.x[k,w] for k in self.process)==1 for w in self.periodsN2),name="const35")
        self.modelL2.addConstrs((self.M3[k,w]<=self.G*self.x[k,w] for k in self.process for w in self.periodsN2),name="const36")
        self.modelL2.addConstrs((self.Ib[self.periodsN2[w]]==self.Ib[self.periodsN2[w-1]]+gp.quicksum(self.fiber.at[m,self.periodsN2[w]]/(1-self.dataHorizon.at[self.periodsN2[w],'ub'])*self.M1[m,self.periodsN2[w]] for m in self.sourceRaw)-self.Mb[self.periodsN2[w]] for w in range(1,len(self.periodsN2))),name="const37")
        self.modelL2.addConstr((self.Ib[self.periodsN2[0]]==gp.quicksum(self.fiber.at[m,self.periodsN2[0]]/(1-self.dataHorizon.at[self.periodsN2[0],'ub'])*self.M1[m,self.periodsN2[0]] for m in self.sourceRaw)-self.Mb[self.periodsN2[0]]),name="const37_0")
        self.modelL2.addConstrs((self.Ib[w]>=gp.quicksum(self.fiber.at[m,w]/(1-self.dataHorizon.at[w,'ub'])*self.M1[m,w] for m in self.sourceRaw)*self.eb for w in self.periodsN2),name="const38")
        self.modelL2.addConstr((self.Ib[self.periodsN2[-1]]>=self.Epb),name="const39")
        self.modelL2.addConstrs((self.rc*self.Mb[w]==self.VAP[w] for w in self.periodsN2),name="const40")
        self.modelL2.addConstrs((self.VAP[w]>=gp.quicksum(self.M3[k,w] for k in self.process)*self.cfvap + (self.EG[w]/self.rcf) for w in self.periodsN2),name="const41a")
        self.modelL2.addConstrs((gp.quicksum(self.M3[k,w] for k in self.process)*self.cfvap + (self.EG[w]/self.rcf) >= gp.quicksum(self.dataProdMill.at[p,'cvap']*self.PU2[p,w] for p in self.ProductsMill) for w in self.periodsN2),name="const41b")
        self.modelL2.addConstrs((self.EG[w]-(gp.quicksum(self.cfe*self.M3[k,w] for k in self.process)+gp.quicksum(self.dataProdMill.at[p,'cve']*self.PU2[p,w] for p in self.ProductsMill))==self.EE[w] for w in self.periodsN2),name="const42")
        self.modelL2.addConstrs((self.VAP[w]<=self.vapmax*self.dataHorizon.at[w,'n'] for w in self.periodsN2),name="const43")
        self.modelL2.addConstrs((self.EG[w]<=self.egmax*self.dataHorizon.at[w,'n'] for w in self.periodsN2),name="const44")
        self.modelL2.setObjective((gp.quicksum(self.PU2[pa,w]/self.P2*self.R2*modelL1.dataProd.at['ConvCrpaaa',pa] for w in self.periodsN2 for pa in self.SugarMill)
                            +gp.quicksum(self.PU2[pe,w]/self.P3*self.R3*modelL1.dataProd.at['ConvCrpaaa',pe] for w in self.periodsN2 for pe in self.EthanolMill)
                            +gp.quicksum(self.PU2[pm,w]/self.P4*self.R4*modelL1.dataMills.at['Matr',self.millName] for w in self.periodsN2 for pm in self.MolassesMill)
                            +gp.quicksum(self.ve*self.EE[w] for w in self.periodsN2)
                            -gp.quicksum(self.Cm.at[m,w]*self.M1[m,w] for m in self.sourceRaw for w in self.periodsN2)
                            -gp.quicksum(self.Cf.at[f,w]*self.M2[f,w] for f in self.transport for w in self.periodsN2)
                            -gp.quicksum(self.Ck.at[k,w]*self.M3[k,w] for k in self.process for w in self.periodsN2)
                            -gp.quicksum(modelL1.ca*self.ImMill[p,w] for p in self.ProductsMill for w in self.periodsN2)),gp.GRB.MAXIMIZE)
        self.modelL2.update()
        return self.modelL2
