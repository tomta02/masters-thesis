# some functions for processing and comparing yaml files

# function "yaml_compare": comparing content of 2 yaml files
# comparing the content of both yaml files to make sure all imgs underwent both thresholding and binarization and did not forget anything :-)

# function "filter_yaml_3rd_layer":
# not all markers that were used for CODEX imaging worked, therefore, some markers were thresholded but not binarized, because we decided to drop them in the end.
# Before comparing the thresholding and binarization yaml files, subset them by markers that actually worked and that we are actually interested in - otherwise warning is raised for markers present in thresholding yaml not present in binarization
# yaml that we are not interested in at all anymore. Protein markers are given in third layer of yaml file, preserve first and second level structure, altering third lvl only

import yaml

def load_yaml(pat):
    '''load yaml file as dict'''
    with open(pat, "r") as f:
        return yaml.safe_load(f)


def flt_yaml_3rd_layer(yamldict, markers2keep):
    '''
    In our yml files: 3rd layer = protein markers
    remove protein markers from yaml that were disregarded due to their quality.
    
    In: dictionary from yaml file (from load_yaml), list with markers to keep
    Out: yaml dict with only the markers of interest
    '''
    if not isinstance(yamldict, dict):
        return yamldict
    flt_first_layer = {}
    for top_key, second_lvl_dict in yamldict.items():
        if isinstance(second_lvl_dict, dict):
            flt_first_layer[top_key] = {}
            for samplename, markers_dict in second_lvl_dict.items():
                flt_first_layer[top_key][samplename] = {
                    k: v for k, v in markers_dict.items() if k in markers2keep}

    if all(not v for v in flt_first_layer.values()):
        print("[!] subdicts empty")
        
    return flt_first_layer


def yaml_compare(dict1, dict2, pat=""):
    '''
    compare 2 yaml files (thresh yaml and binarization yaml) recursively and
    1.) warn if keys are missing in either yaml (should be the same but just in case)
    2.) warn if one has a value while other has missing value at given position

    In: two yaml files to compare, loaded as dicts
    '''
    warn = []
    empty_vals = (None, '')
    
    if isinstance(dict1, dict) and isinstance(dict2, dict):
        allkeys = set(dict1.keys()) | set(dict2.keys())
        for k in allkeys:
            npat = f"{pat}.{k}" if pat else k
            if k not in dict1:
                warn.append(f"[!] key '{npat}' is missing in yaml1")
            elif k not in dict2:
                warn.append(f"[!] key '{npat}' is missing in yalm2")
            else: 
                v1 = dict1.get(k)
                v2 = dict2.get(k)
                warn.extend(yaml_compare(v1, v2, npat))
  
    else:
        if (dict1 in empty_vals) != (dict2 in empty_vals):
            warn.append(f"[!] mismatch at '{pat}': yaml1 has '{dict1}', yaml2 has '{dict2}'")

    return warn
            

if __name__ == "__main__":
    yaml_t = load_yaml("/g/saka/Tatjana/analysis/codex/yamlfiles/thresholds.yaml")
    yaml_b = load_yaml("/g/saka/Tatjana/analysis/codex/yamlfiles/binarization.yaml")

    # get list of markers of interes
    with open("/g/saka/Tatjana/analysis/codex/codex_markers.txt", "r") as f:
        mk = f.read().splitlines()

    yaml_t_flt = flt_yaml_3rd_layer(yaml_t, mk)
    yaml_b_flt = flt_yaml_3rd_layer(yaml_b, mk)

    
    diffs = yaml_compare(yaml_t_flt, yaml_b_flt)
    if diffs:
        print("\n".join(diffs))
    else:
        print("both yml files match")


