#!/usr/bin/env python3

import argparse
from bisect import bisect_left

def parse_pseudoalignments(mapping_file, contig_map, rnames):
    num_found = 0
    num_mapped = 0
    num_missed = 0
    num_positive = 0
    total_mappings = 0
    num_spurious = 0
    is_themisto = len(rnames) > 0
    target_idx = 1 if is_themisto else 2

    with open(mapping_file, 'r') as ifile:
        for i,l in enumerate(ifile):
            if i % 500000 == 0 and i > 0:
                print(f"parsed {i} pseudoalignments")
            num_mapped += 1
            ntargets = 0
            toks = None
            if is_themisto:
                toks = l.split(' ')
                ref_key = rnames[int(toks[0])]
                ntargets = len(toks) - 1
            else:
                toks = l.split('\t')
                ref_key = toks[0].split(':')[-2]
                ntargets = int(toks[1])
            total_mappings += ntargets
            ref_id = contig_map.get(ref_key)
            if ref_id is not None:
                num_positive += 1
                if ntargets > 1:
                    targets = list(map(int, toks[target_idx:]))
                    i = bisect_left(targets, ref_id)
                    if i != len(targets) and targets[i] == ref_id:
                        num_found += 1
                    else:
                        num_missed += 1
                else:
                    num_missed += 1
            else:
                # if not ref_key.startswith("NC"):
                #     print(ref_key)
                if ntargets > 0:
                    num_spurious += 1
                # else:
                #     true_neg += 1
    return {'num_found' : num_found, 'observed_reads' : num_mapped, 
            'num_missed' : num_missed, 'total_mappings' : total_mappings, 
            'num_positive' : num_positive, 'num_spurious' : num_spurious,
            "TPR" : float(num_found/num_positive), "FPR": float(num_spurious/(num_mapped-num_positive)),
            "Precision": float(num_found/(num_found+num_spurious)),
            "Hits Per Read": float(total_mappings/num_mapped)}

def parse_read_names(rnf):
    rn = []
    with open(rnf, 'r') as ifile:
        for l in ifile:
            rn.append(l.rstrip().split(':')[0])
    return rn

def build_contig_map(contig_map_file):
    cm = {}
    with open(contig_map_file, 'r') as ifile:
        for l in ifile:
            toks = l.strip().rstrip().split(':')
            ctg = toks[-1]
            ref_id = int(toks[1].split('_')[0]) - 1 
            cm[ctg] = ref_id
    return cm

def evaluate_mappings(args):
    import json
    contig_map = build_contig_map(args.contig_map)
    print(f"Read contig map, contains {len(contig_map)} contigs")
    rnames = []
    if args.read_names is not None:
        print(f"Reading read name map {args.read_names} (to process themisto output)")
        rnames = parse_read_names(args.read_names)
        print(f"read {len(rnames)} names")

    res = parse_pseudoalignments(args.mapping_file, contig_map, rnames)
    print(f"results : {res}")
    print(f"recall : {res['num_found']/float(res['num_positive'])}")
    print(f"hits-per-read : {res['total_mappings']/float(res['observed_reads'])}")
    
    res_json = json.dumps(res, indent=2)
    with open(args.output, 'w') as ofile:
        ofile.write(res_json)

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("contig_map")
    parser.add_argument("mapping_file")
    parser.add_argument("--names", dest="read_names", help="optional file of read names (for themisto)")
    parser.add_argument("-o", dest="output", help="where output json dictionary should be written")
    args = parser.parse_args()
    evaluate_mappings(args)

