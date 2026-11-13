import pandas as pd
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.ticker as tk
import gurobipy_pandas as gppd
import gurobipy as gp
import itertools as it
import classPASUC as ps
#import classMPBRAP as mr

if __name__ == '__main__':

    originalPasuc1 = ps.PASUCL1()
    originalPasuc2 = ps.PASUCL2()

    originalPasuc1.readDataCoop("./data/data.xlsx")
    modelPasuc = originalPasuc1.createModelL1("PASUC/N1")
    modelPasuc.setParam("MIPGap",1e-10)
    modelPasuc.optimize()

    originalPasuc2.readDataMill("./data/dataMill2.xlsx")
    modelPasuc2 = originalPasuc2.createModelL2("MPBRAP/N2",originalPasuc1)
    modelPasuc2.setParam("MIPGap",1e-10)
    modelPasuc2.optimize()

    
    # pd.concat([    
    # #pd.DataFrame([[originalPasuc2.x[k,w].X for w in originalPasuc2.periodsN2]for k in originalPasuc2.process],columns=originalPasuc2.periodsN2, index=["X"+ k for k in originalPasuc2.process]),
    # pd.DataFrame([[sum((originalPasuc2.process.index(k)+1)*originalPasuc2.x[k,w].X for k in originalPasuc2.process) for w in originalPasuc2.periodsN2]],columns=originalPasuc2.periodsN2, index=["K"]),
    # #gp.quicksum(self.txvin.at[k]*self.A.at[(pd,k),self.periodsN2[w]]*self.M3[k,self.periodsN2[w]] for pd in self.EthanolMill  for k in self.process)
    # pd.DataFrame([[sum(originalPasuc2.txvin.at[k]*originalPasuc2.A.at[(pd,k),w]*originalPasuc2.M3[k,w].X for k in originalPasuc2.process for pd in originalPasuc2.EthanolMill) for w in originalPasuc2.periodsN2]],columns=originalPasuc2.periodsN2,index=["Xvin"]),
    # pd.DataFrame([[originalPasuc2.Ivin[w].X for w in originalPasuc2.periodsN2]],columns=originalPasuc2.periodsN2,index=["Ivin"]),
    # #pd.DataFrame([[originalPasuc2.Ovin[k,w].X for w in originalPasuc2.periodsN2]for k in originalPasuc2.process],columns=originalPasuc2.periodsN2, index=["Ovin"+ k for k in originalPasuc2.process]),
    # pd.DataFrame([[sum(originalPasuc2.Ovin[k,w].X for k in originalPasuc2.process) for w in originalPasuc2.periodsN2]],columns=originalPasuc2.periodsN2, index=["OvinTotal"]),
    # pd.DataFrame([[sum(originalPasuc2.OMvin[k]*originalPasuc2.Ovin[k,w].X for k in originalPasuc2.process) for w in originalPasuc2.periodsN2]],columns=originalPasuc2.periodsN2, index=["OMvin"]),
    # pd.DataFrame([[originalPasuc2.Fvin[w].X for w in originalPasuc2.periodsN2]],columns=originalPasuc2.periodsN2,index=["Fvin"]),
    # pd.DataFrame([[sum(originalPasuc2.txcake.at[k]*originalPasuc2.A.at[(pa,k),w]*originalPasuc2.M3[k,w].X for pa in originalPasuc2.SugarMill  for k in originalPasuc2.process) for w in originalPasuc2.periodsN2]],columns=originalPasuc2.periodsN2,index=["Xcake"]),
    # pd.DataFrame([[originalPasuc2.Icake[w].X for w in originalPasuc2.periodsN2]],columns=originalPasuc2.periodsN2,index=["Icake"]),
    # pd.DataFrame([[originalPasuc2.Ocake[w].X for w in originalPasuc2.periodsN2]],columns=originalPasuc2.periodsN2,index=["Ocake"]),
    # pd.DataFrame([[originalPasuc2.Ibiogas[w].X for w in originalPasuc2.periodsN2]],columns=originalPasuc2.periodsN2,index=["Ibiogas"]),
    # pd.DataFrame([[originalPasuc2.Otgen[w].X for w in originalPasuc2.periodsN2]],columns=originalPasuc2.periodsN2,index=["Otgen"]),
    # pd.DataFrame([[originalPasuc2.Obiogas[w].X for w in originalPasuc2.periodsN2]],columns=originalPasuc2.periodsN2,index=["Obiogas"]),
    # pd.DataFrame([[originalPasuc2.Ofleet[w].X for w in originalPasuc2.periodsN2]],columns=originalPasuc2.periodsN2,index=["Ofleet"]),
    # pd.DataFrame([[originalPasuc2.Ibiomet[w].X for w in originalPasuc2.periodsN2]],columns=originalPasuc2.periodsN2,index=["Obiomet"]),
    # pd.DataFrame([[originalPasuc2.Q[w].X for w in originalPasuc2.periodsN2]],columns=originalPasuc2.periodsN2,index=["Q"]),
    # pd.DataFrame([[originalPasuc2.Opipe[w].X for w in originalPasuc2.periodsN2]],columns=originalPasuc2.periodsN2,index=["Opipe"]),
    # pd.DataFrame([[originalPasuc2.Ilost[w].X for w in originalPasuc2.periodsN2]],columns=originalPasuc2.periodsN2,index=["Ilost"])
    # ]).to_excel("c:/Users/marko/source/repos/PASUC/Vars.xlsx")

