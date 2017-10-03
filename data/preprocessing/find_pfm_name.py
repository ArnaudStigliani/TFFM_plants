# -*-coding:Utf-8 -*
import requests, sys

list_name="total_list"


TF=str()
with open(list_name,"r") as f1:
    TF=f1.read()
    

TF=TF.split("\n")
#print(TF)

 
api_root = "http://hfaistos.uio.no:8000/api/v1/"
# TF=list()
# TF.append("ERF105")
with open("totaal_output","w") as f1:
    for elt in TF:
        
        query = "matrix?format=json&search="+elt

        result = requests.get(api_root+query)
    
        if not result.ok:
            print(elt)
            r.raise_for_status()
            sys.exit()
        elif (len(result.json()["results"]) == 0):
            f1.write(elt+"\n")
        else :
            decoded = result.json()
            for results_decoded in decoded["results"]:
                if((results_decoded['name']).lower())==elt.lower():
                    line=elt + "\t" + results_decoded["base_id"]+ "\t" + results_decoded["version"]+"\n"
            f1.write(line)
