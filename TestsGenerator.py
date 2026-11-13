import pandas as pd
from openpyxl import load_workbook
import shutil
pd.options.mode.chained_assignment = None  # default='warn'

if __name__ == '__main__':

    
    testVar = pd.read_excel("./updatedData/dataTestVar.xlsx", sheet_name = "dataVar", index_col=0, header=0, converters={'IDTest': str})
    #print(testVar)
    for id in testVar.index:
        shutil.copy("./updatedData/dataMillBioP4C1.xlsx", "./dataGen/dataMillBio"+str(id)+".xlsx")
        #Ajuste da planilha Biodigester (parametros do biodigestor)
        dataBio = pd.read_excel("./dataGen/dataMillBio"+str(id)+".xlsx",sheet_name="Biodigester", index_col=0, header=None)
        dataBio[2] = ''#Sobrescreve a coluna de comentários do arquivo base com ''
        #print(dataBio)
        dataBio[1]["txBiogasVin"] *= testVar.at[id,"txBiogas"]
        dataBio[1]["txBiogasCake"] *= testVar.at[id,"txBiogas"]
        dataBio[1]["digMaxVin"] *= testVar.at[id,"digMax"]
        dataBio[1]["digMaxCake"] *= testVar.at[id,"digMax"]
        dataBio[1]["deltaOMvin"] *= testVar.at[id,"deltaOM"]
        dataBio[1]["deltacake"] *= testVar.at[id,"deltaOM"]
        dataBio[1]["tgmax"] *= testVar.at[id,"tgmax"]
        dataBio[1]["txbiomet"] = testVar.at[id,"txbiomet"]
        #dataBio[1]["rcftg"] *= testVar.at[id,"rcftg"]
        dataBio[1]["txConvCake"] *= testVar.at[id,"txConv"]
        dataBio[1]["EpbCake"] = testVar.at[id,"EpbCake"]*dataBio[1]["digMaxCake"]
        dataBio[1]["I0Cake"] = testVar.at[id,"I0Cake"]*dataBio[1]["digMaxCake"]
        dataBio[1]["biogasPlantCost"] =testVar.at[id,"biogasPlantCost"]
        #print(dataBio)
        writer = pd.ExcelWriter("./dataGen/dataMillBio"+str(id)+".xlsx",engine='openpyxl',mode='a', if_sheet_exists='overlay')
        dataBio.to_excel(writer, sheet_name='Biodigester',header=False)
        writer.close()
        #Ajuste da planilha Process (parametros de produção e carga organica de vinhaca e torta)
        generalN2 = pd.read_excel("./dataGen/dataMillBio"+str(id)+".xlsx", sheet_name = "General",index_col=0, header=None)
        numProc = generalN2[1]["numProc"]
        process = ['Proc'+str(i) for i in range(1,numProc+1)]
        dataProcess = pd.read_excel("./dataGen/dataMillBio"+str(id)+".xlsx",sheet_name="Process", index_col=0, header=0)
        for k in process:
            dataProcess["txvin"][k]*= testVar.at[id,"txprod"]
            dataProcess["txcake"][k]*= testVar.at[id,"txprod"]
            #dataProcess["OMvin"][k]*= testVar.at[id,"OMvin"]
            #dataProcess["Tc"][k]*= testVar.at[id,"Tc"]
        #print(dataProcess)
        writer = pd.ExcelWriter("./dataGen/dataMillBio"+str(id)+".xlsx",engine='openpyxl',mode='a', if_sheet_exists='overlay')
        dataProcess.to_excel(writer, sheet_name='Process')
        writer.close()
        numWeeks = generalN2[1]["Periods"]
        numSeason = generalN2[1]["Season"]
        periodsSeason = ['Week'+str(t) for t in range(1,numSeason+1)]
        periodsOffSeason = ['Week'+str(t) for t in range(numSeason+1,numWeeks+1)]
        #print(periodsSeason,periodsOffSeason)
        dataHorizon = pd.read_excel("./dataGen/dataMillBio"+str(id)+".xlsx",sheet_name="Horizon", index_col=0, header=0)
        for w in periodsSeason:
            dataHorizon["txConvVin"][w]*= testVar.at[id,"txConv"]
            dataHorizon["Dbmmin"][w]*= testVar.at[id,"Dbm"]
            dataHorizon["Dbmmax"][w]*= testVar.at[id,"Dbm"]
            dataHorizon["DbmPipe"][w]*= testVar.at[id,"Dbm"]
            dataHorizon["VbmPipe"][w]*= testVar.at[id,"Vbm"]
            dataHorizon["VbmDiesel"][w]*= testVar.at[id,"Vbm"]
            dataHorizon["Vebio"][w]*=testVar.at[id,"Vene"]
        for w in periodsOffSeason:
            dataHorizon["txConvVin"][w]*= testVar.at[id,"txConv"]
            dataHorizon["Dbmmin"][w]*= testVar.at[id,"Dbm"]
            dataHorizon["Dbmmax"][w]*= testVar.at[id,"Dbm"]
            dataHorizon["DbmPipe"][w]*= testVar.at[id,"Dbm"]
            dataHorizon["VbmPipe"][w]*= testVar.at[id,"Vbm"]
            dataHorizon["VbmDiesel"][w]*= testVar.at[id,"Vbm"]
            dataHorizon["Vebio"][w]*=testVar.at[id,"Vene"]
        #print(dataHorizon)
        writer = pd.ExcelWriter("./dataGen/dataMillBio"+str(id)+".xlsx",engine='openpyxl',mode='a', if_sheet_exists='overlay')
        dataHorizon.to_excel(writer, sheet_name='Horizon')
        writer.close()
        #dataGeneral = pd.read_excel("dataGen\dataMillBio"+str(id)+".xlsx",sheet_name="General", index_col=0, header=None)
        #dataGeneral[1]["Ve"]*=testVar.at[id,"Vene"]
        #writer = pd.ExcelWriter("dataGen\dataMillBio"+str(id)+".xlsx",engine='openpyxl',mode='a', if_sheet_exists='overlay')
        #dataGeneral.to_excel(writer, sheet_name='General',header=None,columns=None)
        #writer.close()
