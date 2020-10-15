import csv
import os
from collections import defaultdict
from collections import Counter
import pandas as pd

from reportfunk.funks import io_functions as qcfunk

def configure_input_query(config):

    query_metadata = config["query"]
    patient_id_col = config["patient_id_col"]
    path_to_dir = config["outdir"]

    new_query = os.path.join(path_to_dir, "query_with_groups.csv") 

    patient_to_group = {}
    count = 1
    single_line = False
    with open(query_metadata) as f:
        reader = csv.DictReader(f)
        fieldnames=reader.fieldnames
        data = [r for r in reader]
        for seq in data:
            
            if len(fieldnames) == 1:
                patient = "Patient A"
                single_line = True
            else:
                patient = seq[patient_id_col]

            if patient not in patient_to_group.keys():
                patient_to_group[patient] = count
                count += 1

    with open(query_metadata) as csv_input:
        with open(new_query, "w") as csv_output:
            reader = csv.DictReader(csv_input)
            data = [r for r in reader]
            
            if single_line:
                fieldnames.append("patient_id")
            fieldnames.append("group")
            
            writer = csv.DictWriter(csv_output, fieldnames = fieldnames)
            writer.writeheader()
            
            for row in data:
                if single_line:
                    row["patient_id"] = patient
                    row["group"] = patient_to_group[patient]
                else:
                    row["group"] = patient_to_group[row[patient_id_col]]
                
                writer.writerow(row)

    config["query"] = new_query

def make_patient_files(config):

    path_to_dir = config["tempdir"]
    query = config["query"]
    name_col = config["input_column"]
    patient_col = config["patient_id_col"]

    try:
        os.mkdir(os.path.join(path_to_dir,"patient_data"))
    except FileExistsError:
        pass

    patient_dict = defaultdict(list)

    with open(query) as f:
        reader = csv.DictReader(f)
        data = [r for r in reader]
        for row in data:
            name = row[name_col]
            patient = row[patient_col]

            patient_dict[patient].append(name)

    for k in patient_dict.keys():
        fw = open(os.path.join(path_to_dir,"patient_data",f"{k}_taxa.txt"), 'w')
        fw.write('name\n')
        for l in patient_dict[k]:
            fw.write(l + "\n")
        fw.close()

def reinfection_args_to_config(args, config, defaultdict):

    patient_id_col = qcfunk.check_arg_config_default("patient_id_col", args.patient_id_col, config, defaultdict)
    config["patient_id_col"] = patient_id_col

def sort_into_patients(query_dict):

    patient_to_taxa = defaultdict(list)
    patient_to_tree = defaultdict(set)
    patient_list = set()
    
    for tax in query_dict.values():
        
        patient = tax.attribute_dict["patient"]
        patient_to_taxa[patient].append(tax)
        patient_to_tree[patient].add(tax.tree)
        patient_list.add(patient)


    patient_list = list(patient_list)

    return patient_list, patient_to_taxa, patient_to_tree

def make_summary_table(focal_patient, patient_to_taxa, tree_to_tips, full_taxon_dict):

    df_dict = defaultdict(list)
    query_taxa = patient_to_taxa[focal_patient]

    for taxa in query_taxa:
        same_adm2 = 0
        same_epiweek = 0
        tree = taxa.tree
        tips = tree_to_tips[taxa.tree]
        focal_epiweek = taxa.epiweek
        focal_adm2 = taxa.attribute_dict["adm2"]

        total_tips = len(tips) - len(query_taxa)

        for tip in tips:
            if tip not in query_taxa:
                if tip in full_taxon_dict.keys():
                    taxon_obj = full_taxon_dict[tip]
                    adm2 = taxon_obj.attribute_dict['adm2']
                    if focal_adm2 == "NA" or focal_adm2 == "":
                        same_adm2 = "NA"
                    else:
                        if adm2 == focal_adm2 and adm2 != "" and adm2 != "NA":
                            same_adm2 += 1
                    if type(focal_epiweek) == str:
                        same_epiweek = "NA"
                    else:
                        if type(taxon_obj.epiweek) != str:
                            if taxon_obj.epiweek == focal_epiweek:
                                same_epiweek += 1
        
        if type(same_adm2) != str:
            adm2_perc = (same_adm2/total_tips)*100
        else:
            adm2_perc = "NA"
        if type(same_epiweek) != str:
            epiweek_perc = (same_epiweek/total_tips)*100
        else:
            epiweek_perc = "NA"

        df_dict["Query sequence"].append(taxa.name)
        df_dict["Adm2"].append(focal_adm2)
        df_dict["Epiweek"].append(focal_epiweek)
        df_dict["Tree"].append(taxa.tree)
        df_dict["Percentage of tips in tree with the same adm2"].append(str(adm2_perc))
        df_dict["Percentage of tips in tree with the same epiweek"].append(str(epiweek_perc))

    df = pd.DataFrame(df_dict)
    df.set_index("Query sequence", inplace=True)

    return df




            

