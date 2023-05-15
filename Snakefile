from os.path import join
import os
from norm_guilio_newton import norm_single_contig, df_list
# configfile:
#     "config.json"
configfile:
    config["cuttlefish_json"]


os.environ['TMPDIR'] = '/mnt/scratch4'

inp_path = config["file_path"]
out_fq_path = config["out_fq_path"]
out_psa_path = config["out_psa_path"]
non_samp_path = config["non_samp_path"]
extension=""
if "inp_file_exten" in config.keys():
    extension = config["inp_file_exten"]

inp_query_file = join(inp_path, "filenames_list_{nquery}k_"+extension+".txt")
if inp_path.endswith("tsv") or inp_path.endswith("txt"):
    inp_query_file = inp_path
print(inp_query_file)

num_reads = config["reads_per_reference"]

nqueries = [5,10,50,100,150]
nqueries = [50]
# ratios = [5,90]
ratios = [5,10,15,30,45,60,75,90]
# ratios = [5,10,15,25,30,45,50,60,75,90]

out_file_fq = join(out_fq_path, "sim_reads_{nquery}_ratio={ratio}.fq")
out_seeds = join(out_fq_path, "seeds_{nquery}_ratio={ratio}")

map_psa_k_out = join(out_psa_path, "map_kallisto_nquery={nquery}_ratio={ratio}.out")
map_psa_skip_out = join(out_psa_path, "map_skip_nquery={nquery}_ratio={ratio}.out")
map_psa_hybrid_out = join(out_psa_path, "map_hybrid_thresh={hyb_thresh}_nquery={nquery}_ratio={ratio}.out")
map_psa_full_int_out = join(out_psa_path, "map_fullint_nquery={nquery}_ratio={ratio}.out")

map_psa_k_json = join(out_psa_path, "map_kallisto_nquery={nquery}_ratio={ratio}.json")
map_psa_skip_json = join(out_psa_path, "map_skip_nquery={nquery}_ratio={ratio}.json")
map_psa_hybrid_json = join(out_psa_path, "map_hybrid_thresh={hyb_thresh}_nquery={nquery}_ratio={ratio}.json")
map_psa_full_int_json = join(out_psa_path, "map_fullint_nquery={nquery}_ratio={ratio}.json")

psa_bin_path = config["psa_path"]
cuttlefish_index_path = config["cuttlefish_index_path"]
hybdrid_thresh  =config["hybrid_thresh"]

rule all:
    input:
        inp_fq = expand(out_file_fq, nquery = nqueries, ratio = ratios),
        inp_seed = expand(out_seeds, nquery = nqueries, ratio = ratios),
        # map_contig = config['contig_map'],
        psa_kjson = expand(map_psa_k_json, nquery = nqueries, ratio = ratios),
        psa_sjson = expand(map_psa_skip_json, nquery = nqueries, ratio = ratios),
        psa_hjson = expand(map_psa_hybrid_json, nquery = nqueries, ratio = ratios, hyb_thresh = hybdrid_thresh),
        psa_fjson = expand(map_psa_full_int_json, nquery = nqueries, ratio = ratios)

rule run_psa_eval:
    input:
        psa_kjson = expand(map_psa_k_json, nquery = nqueries, ratio = ratios),
        psa_sjson = expand(map_psa_skip_json, nquery = nqueries, ratio = ratios),
        psa_hjson = expand(map_psa_hybrid_json, nquery = nqueries, ratio = ratios, hyb_thresh = hybdrid_thresh),
        psa_fjson = expand(map_psa_full_int_json, nquery = nqueries, ratio = ratios)

rule _run_spsa:
    input:
        out_file_fq
    output:
        map_psa_skip_json,
    params:
        map_out = map_psa_skip_out,
        psa_bin_path = psa_bin_path,
        cind = cuttlefish_index_path,
        map_contig = config['contig_map'],
        python_eval_script = config["eval_script_path"]
    shell:
        """
            {params.psa_bin_path} -i {params.cind} -q {input} -t 4 -o {params.map_out} --skipping
            python {params.python_eval_script} {params.map_contig} {params.map_out} -o {output}
            rm {params.map_out}
        """
        
rule _run_hpsa:
    input:
        out_file_fq
    output:
        map_psa_hybrid_json
    params:
        map_out = map_psa_hybrid_out,
        psa_bin_path = psa_bin_path,
        cind = cuttlefish_index_path,
        map_contig = config['contig_map'],
        thresh = hybdrid_thresh,
        python_eval_script = config["eval_script_path"]
    shell:
        """
            {params.psa_bin_path} -i {params.cind} -q {input} -t 4 -o {params.map_out} --threshold {params.thresh}
            python {params.python_eval_script} {params.map_contig} {params.map_out} -o {output}
            rm {params.map_out}
        """

rule _run_fpsa:
    input:
        out_file_fq
    output:
        map_psa_full_int_json
    params:
        map_out = map_psa_full_int_out,
        psa_bin_path = psa_bin_path,
        cind = cuttlefish_index_path,
        map_contig = config['contig_map'],
        python_eval_script = config["eval_script_path"]
    shell:
        """
            {params.psa_bin_path} -i {params.cind} -q {input} -t 4 -o {params.map_out}
            python {params.python_eval_script} {params.map_contig} {params.map_out} -o {output}
            rm {params.map_out}
        """

rule _run_kpsa:
    input:
        out_file_fq
    output:
        map_psa_k_json,
    params:
        map_out = map_psa_k_out,
        psa_bin_path = psa_bin_path,
        cind = cuttlefish_index_path,
        map_contig = config['contig_map'],
        python_eval_script = config["eval_script_path"]
    shell:
        """
            {params.psa_bin_path} -i {params.cind} -q {input} -t 4 -o {params.map_out} --skipping-kallisto
            python {params.python_eval_script} {params.map_contig} {params.map_out} -o {output}
            rm {params.map_out}
        """

ref_files = config['parameters info']['input']
ref_files = ref_files.split(", ")
print(len(ref_files))
# # assert(len(ref_files)==50000)
print(ref_files[0:10])

if config["guilio_ind"]:
    ref_files = list(map(norm_single_contig, ref_files))
    print(ref_files[0:10])

rule map_contig:
    input:
        join(config['out_path'],'ref_files')
    output:
        config['contig_map']
    params:
        out_path=config['out_path']
    shell:
        "bash map_contig.sh {output} {input} {params.out_path}"

checkpoint ref_files_create:
    input:ref_files[0]
    output:join(config['out_path'],'ref_files')
    run:
        with open(output[0], "w") as f:
            for file in ref_files:
                f.write(f"{file}\n")

rule run_sim:
    input:
        inp_fq = expand(out_file_fq, nquery = nqueries, ratio = ratios),
        inp_seed = expand(out_seeds, nquery = nqueries, ratio = ratios)
rule _run_sim:
    input:
        inp_query_file
    output:
        fq = out_file_fq,
        seed = out_seeds
    params:
        non_samp_path = non_samp_path,
        ratio = lambda wildcards: wildcards.ratio,
        out_path = out_fq_path,
        nquery =  lambda wildcards: wildcards.nquery,
        num_reads = num_reads

    shell:
        """
            bash sim_reads.sh {input} {params.non_samp_path} {params.ratio} {params.out_path} {params.nquery} {params.num_reads}
        """

# rule inp_files:
#     input:expand(inp_query_file, nquery = nqueries)