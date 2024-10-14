# import sys
# sys.path.append( '/Users/adam/miniforge3/envs/test_env/lib/python3.9/site-packages')
# from ccdc.search import TextNumericSearch
from ccdc.search import TextNumericSearch
import json
import copy
from rcsbsearchapi import rcsb_attributes as attrs
from rcsbsearchapi.search import ChemSimilarityQuery
import requests
from urllib.parse import quote
def search_name(str, count, total_count):
    print(f"searching for {str}, #{count}/{total_count}")
    text_numeric_search = TextNumericSearch()
    text_numeric_search.add_compound_name(str)
    results = text_numeric_search.search()
    print(f"found {len(results)} results")
    return results

def get_ccdc():
    
    with open('../data/collatedSynonyms.json', "r") as file:
        data = json.load(file)

    data_copy = copy.deepcopy(data)

    total_count = 0
    for i in range(len(data)):
        total_count += 1
        total_count += len(data[i]["synonyms"])

    current_count = 0
    # for result in searchName('2-Hydroxy-1-naphthaldehyde', 1):
        # print(result.identifier)
    for i in range(len(data)):
        current_count += 1
        datum = data[i]
        hits = []
        # print(datum["drugName"])
        # hits.append( {datum["drugName"]: searchName(datum["drugName"])})
        results = search_name(datum["drugName"], current_count, total_count)
        identifiers = []
        for result in results:
            identifiers.append(result.identifier)
        hits.append({datum["drugName"]: identifiers})
        for synonym in datum["synonyms"]:
            current_count += 1
            # print(synonym)
            identifiers = []
            results = search_name(synonym, current_count, total_count)
            for result in results: 
                print(result.identifier)
                identifiers.append(result.identifier)
            hits.append({synonym: identifiers})
        data_copy[i]["hits"] = hits


    with open("../data/collatedSynonyms_withHits.json", "w") as file:
        json.dump(data_copy, file)

def search_pdb(search_str, smiles_str):
    print(f'searching pdb for {search_str}')
    # q1 = AttributeQuery(attribute ="chem_comp.name", operator="contains_phrase", value=search_str, service="text_chem")
    # q2 = AttributeQuery(attribute ="rcsb_chem_comp_synonyms.name", operator="contains_phrase", value=search_str, service="text_chem")
    # q3 = AttributeQuery(attribute ="rcsb_id", operator="exact_match", value=search_str, service="text_chem")
    # q4 = AttributeQuery(attribute="struct.title", operator="contains_words", value=search_str, service="text")
    # # q2 = attrs.rcsb_chem_comp_related.resource_accession_code =="7L8"
    # query = q1 | q2 | q3 | q4
    query = ChemSimilarityQuery(value=smiles_str, query_type="descriptor", descriptor_type="SMILES", match_type="fingerprint-similarity")
    results = list(query())
    print(f'found {len(results)} results')
    print(results[:3])
    return results[:3]
    
def get_pdb():
    
    with open("../data/collatedSynonyms.json", "r") as file:
        data = json.load(file)

    current_count = 0
    for datum in data:
        if datum.get("smiles", "none") == "none":
            continue
        current_count += 1
        print(f"{current_count}/{len(data)}")
        datum["pdbHits"] = search_pdb(datum["drugName"], datum["smiles"])

    with open("../data/collatedSynonyms_withPDBHits.json", "w") as file:
        json.dump(data, file)
        file.close()
    

get_pdb()
