''' 

Created on 07/31/2024
 
Author: Nicholas Hadas

This program contains the tool to analyze Annovar output data and identifies 
the pathogenic mutations. This program outputs a .xlsx file. Code has been 
updated from the first Annovar tool to allow for csv inputs. 

'''


def Annovar_Analysis_Sort(Input_File_Path, Output_File_Path):
        """
        This function iterates through the columns of an Annovar ouput data frame (.csv), and adds points to genes 
        for meeting the various thresholds of multiple pathogenicity predictors. The function sorts the 
        genes in increasing normalized score and returns the first ten sorted genes while outputting  the sorted 
        .csv file

        Annovar_Analysis_Sort: Str Str -> None 
        """
        import pandas as pd

        A_data = pd.read_csv(Input_File_Path) # Open File 
        A_data["N_Score"] = 0 # Add column to later lump the various scores 
        A_data["N_Score_D"] = 9 # Set normalizing column to 9 (number of predictors)
        for index, row in A_data.loc[:, ['CLNSIG']].iterrows(): # ClNSIG = Pathogenic 
                if "na" == row['CLNSIG'] or "." == row['CLNSIG']:
                        A_data.loc[index,'N_Score_D'] -= 1 
                elif "pathogenic" in row['CLNSIG'].lower(): 
                        A_data.loc[index,'N_Score'] += 1
        for index, row in A_data.loc[:, ['PopFreqMax']].iterrows(): # PopFreqMax <= 0.01
                if any(chr.isdigit() for chr in row['PopFreqMax']): 
                        if float(row['PopFreqMax']) <= 0.01: 
                                A_data.loc[index,'N_Score'] += 1
                elif "na"== row['PopFreqMax'] or "." == row['PopFreqMax']:
                        A_data.loc[index,'N_Score_D'] -= 1
        for index, row in A_data.loc[:, ['Polyphen2_HVAR_pred']].iterrows(): # Polyphen2_HVAR_pred = Damaging 
                if "na" in row['Polyphen2_HVAR_pred'] or "." == row['Polyphen2_HVAR_pred']:
                        A_data.loc[index,'N_Score_D'] -= 1 
                elif row['Polyphen2_HVAR_pred'] == "D":
                        A_data.loc[index,'N_Score'] += 1
        for index, row in A_data.loc[:, ['SIFT_pred']].iterrows(): # SIFT_pred = Deleterious 
                if "na" == row['SIFT_pred'] or "." == row['SIFT_pred']:
                        A_data.loc[index,'N_Score_D'] -= 1 
                elif row['SIFT_pred'] == "D":
                        A_data.loc[index,'N_Score'] += 1
        for index, row in A_data.loc[:, ['CADD_phred']].iterrows(): # CADD_phred >= 20, >= 30
                if any(chr.isdigit() for chr in row['CADD_phred']):
                        if 20 <= float(row['CADD_phred']) < 30:
                                A_data.loc[index,'N_Score'] += 0.5
                        elif float(row['CADD_phred']) >= 30:
                                A_data.loc[index,'N_Score'] += 1
                elif "na" == row['CADD_phred'] or "." == row['CADD_phred']:
                        A_data.loc[index,'N_Score_D'] -= 1
        for index, row in A_data.loc[:, ['DANN_score']].iterrows(): # DANN_score >= 0.96
                if any(chr.isdigit() for chr in row['DANN_score']):
                        if float(row['DANN_score']) >= 0.96:
                                A_data.loc[index,'N_Score'] += 1
                elif "na" == row['DANN_score'] or "." == row['DANN_score']:
                        A_data.loc[index,'N_Score_D'] -= 1
        for index, row in A_data.loc[:, ['MetaSVM_pred']].iterrows(): # MetaSVM_pred = Damaging
                if "na" == row['MetaSVM_pred'] or "." == row['MetaSVM_pred']:
                        A_data.loc[index,'N_Score_D'] -= 1 
                elif row['MetaSVM_pred'] == "D":
                        A_data.loc[index,'N_Score'] += 1
        for index, row in A_data.loc[:, ['MetaLR_pred']].iterrows(): # MetaLR_pred = Damaging
                if "na" == row['MetaLR_pred'] or "." == row['MetaLR_pred']:
                        A_data.loc[index,'N_Score_D'] -= 1 
                elif row['MetaLR_pred'] == "D":
                        A_data.loc[index,'N_Score'] += 1
        for index, row in A_data.loc[:, ['M-CAP_pred']].iterrows(): # M-CAP_pred = Deleterious
                if "na" == row['M-CAP_pred'] or "." == row['M-CAP_pred']:
                        A_data.loc[index,'N_Score_D'] -= 1 
                elif row['M-CAP_pred'] == "D":
                        A_data.loc[index,'N_Score'] += 1
        for index, row in A_data.loc[:, ['REVEL']].iterrows(): # REVEL >= 0.625
                if any(chr.isdigit() for chr in row['REVEL']):
                        if float(row['REVEL']) >= 0.625: 
                                A_data.loc[index,'N_Score'] += 1
                elif "na" == row['REVEL'] or "." == row['REVEL']:
                        A_data.loc[index,'N_Score_D'] -= 1
        A_data = A_data.loc[A_data['N_Score_D'] > 0] # Filter results that have no valid data
        A_data['N_Score_P'] = (A_data['N_Score'])/(A_data["N_Score_D"]) # Normalize N_Score in new column 
        A_data = A_data.sort_values(['N_Score_P', 'N_Score_D'], ascending=False) # Sort on decreasing N_Score_P, then decreasing N_Score_D
        A_data.to_excel(Output_File_Path, index=False) # Output csv File 
        print(A_data[['Chr', 'Start', 'End', 'Gene.refGene']].head(11)) # Print first ten most probable results 
        

functions = {
    "Annovar_Analysis_Sort": Annovar_Analysis_Sort
}

import argparse
parser = argparse.ArgumentParser(description = "Annovar Tool")
parser.add_argument('function', metavar = 'function', help = "Specifiy which function to use")
parser.add_argument('input_path', metavar = 'input_path', type = str,  help = "Enter Annovar File Path")
parser.add_argument('output_path', metavar = 'output_path', type = str,  help = "Enter Annovar File Path")
args = parser.parse_args()

name = args.function
functions[name](args.input_path, args.output_path)