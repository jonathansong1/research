
import pandas as pd

data = pd.read_csv("/Users/jonathansong/Research_Folder/clean_data.tsv", sep="\t")
transvar = pd.read_csv("/Users/jonathansong/Research_Folder/transvar_results_2/final_results.tsv", sep="\t")

#transvar_1 = pd.read_csv("/Users/jonathansong/Research_Folder/transvar_results_2/ensembl.tsv", sep="\t")
#transvar_2 = pd.read_csv("/Users/jonathansong/Research_Folder/transvar_results_2/CCDS.tsv", sep="\t")
#transvar_3 = pd.read_csv("/Users/jonathansong/Research_Folder/transvar_results_2/RefSeq.tsv", sep="\t")
#transvar_4 = pd.read_csv("/Users/jonathansong/Research_Folder/transvar_results_2/GENCODE.tsv", sep="\t")
#transvar_5 = pd.read_csv("/Users/jonathansong/Research_Folder/transvar_results_2/UCSC.tsv", sep="\t")
#transvar_6 = pd.read_csv("/Users/jonathansong/Research_Folder/transvar_results_2/AceView.tsv", sep="\t")
#
#temp = [transvar_1, transvar_2, transvar_3, transvar_4, transvar_5, transvar_6]
#to = pd.concat(temp)
#to.to_csv(sep='\t', path_or_buf='/Users/jonathansong/Research_Folder/transvar_results_2/final_results.tsv')


for idx,row in transvar.iterrows():
    try:
        if transvar.at[idx, "input"].split(":")[0] != "[wrap_exception] warning":
            transvar.at[idx, "gene"] = transvar.at[idx, "input"].split(":")[0]
        else:
            transvar.at[idx, "gene"] = transvar.at[idx, "input"].split("_")[3]
    except:
        pass

for idx,row in data.iterrows():
    if type(data.at[idx, "allele"]) != str:
        sub_transvar = transvar[transvar["gene"] == data.at[idx, "Gene_symbol"]]
        for idx_t, row_t in sub_transvar.iterrows():
            if str(sub_transvar.at[idx_t, "coordinates(gDNA/cDNA/protein)"])[:3] == "chr":
                gDNA = sub_transvar.at[idx_t, "coordinates(gDNA/cDNA/protein)"].split("/")[0]
                chrom,pos = gDNA.split("g.")
                newchrompos = chrom + pos[:-3]
                if newchrompos != data.at[idx, "Chrom_pos"]:
                    data.at[idx, "New_Chrom_pos"] = newchrompos
                aaMut = sub_transvar.at[idx_t, "coordinates(gDNA/cDNA/protein)"].split("/")[2]
                aaMut = aaMut.split(".")[1]
                mut = gDNA[-3:]
                loc = gDNA.split(":")[1]
                loc = loc[2:-3]
                temp_list = mut.split(">")
                if len(temp_list)==2:
                    if aaMut == data.at[idx, "AA_mutation"]: 
                        data.at[idx, "allele"] = temp_list[0] + "/" + temp_list[1]
                        data.at[idx, "mapped_gene_works"] = 0
                    else:
                        data.at[idx, "Suggested_AA_mutation"] = aaMut
    
# data.to_csv(sep='\t', path_or_buf='/Users/jonathansong/Research_Folder/clean_data.tsv')
# len(data) - data.count()