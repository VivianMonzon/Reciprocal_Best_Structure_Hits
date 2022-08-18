#!/usr/bin/env python3

###################
#
# @author T. Paysan-Lafosse
#
# @brief This script generates Reciprocal Best Structural Hits
# Requires foldseek (https://github.com/steineggerlab/foldseek) to be installed in /bin of the project root
#
# usage: python find_RBSH.py config.ini --runfoldseek [yes, no]
#
###################

import argparse
import os, re, sys
from configparser import ConfigParser
import subprocess
import shutil

def run(command_list):
    stdout = ""
    out_log = subprocess.run(command_list, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, shell=True)
    stdout = str(out_log.stdout, "utf-8")

    return stdout


def get_model(outputdir, link):
    script_path = os.getcwd()
    filename=link.split('/')[-1]
    tarfilepath=os.path.join(outputdir,filename)
    proteome_id=filename.split('.')[0]
    
    #verify is foldseek db exists for proteome id, if not download data from AlphaFold and create it
    proteome_db=os.path.join(outputdir, f"{proteome_id}_db")
    proteomedir=os.path.join(outputdir,proteome_id)
    if not os.path.isfile(proteome_db) or os.path.getsize(proteome_db) == 0:
        if not os.path.isdir(proteomedir) or not os.listdir(proteomedir):
            os.makedirs(proteomedir, exist_ok=True)
            if not os.path.isfile(tarfilepath): #download models from AF
                cmd = f"""wget -O {tarfilepath} {link};
                        tar xf {tarfilepath} -C {proteomedir};
                    """
                print(f"download AF proteome {proteome_id}")
                run(cmd)
        
            os.remove(tarfilepath)
            delete_cif(proteomedir)
        else:
            print(f"models already downloaded for {proteome_id}")
        
        if os.path.isdir(proteomedir) and os.listdir(proteomedir):
            print(f"generate foldseek db for {proteome_id}")
            cmd = [os.path.join(script_path,f"bin/foldseek createdb {proteomedir}/ {proteome_db}")]
            run(cmd)
    else:
        print(f"foldseek db already exists for {proteome_id}")
    
    print(f"Number of models for {proteome_id}: {len(os.listdir(proteomedir))}")
    return proteome_id


def delete_cif(filedir):
    pattern=r'\.cif\.gz'
    for f in os.listdir(filedir):
        if re.search(pattern, f):
            os.remove(os.path.join(filedir, f))


def get_foldseek(script_dir):
    bindir=os.path.join(script_dir,"bin")
    if not os.path.isdir(bindir):
        print("download foldseek")
        foldseekdir=os.path.join(script_dir, "foldseek/bin")
        tarfile=os.path.join(script_dir, "foldseek-linux-sse41.tar.gz")
        cmd= f"""wget -O {tarfile} --no-check-certificate https://mmseqs.com/foldseek/foldseek-linux-sse41.tar.gz; 
                tar xzf {tarfile}; 
            """
        run(cmd)

        shutil.move(foldseekdir, bindir)
        os.remove(tarfile)
        shutil.rmtree(os.path.join(script_dir, "foldseek"))
    sys.path.append(f"{bindir}/")


def run_foldseek(outputdir, script_path, proteome_id1, proteome_id2):
    print(f"Running foldseek for {proteome_id1}/{proteome_id2}")
    foldseek = os.path.join(script_path, "bin/foldseek")
    m8_file = os.path.join(outputdir, f"{proteome_id1}_{proteome_id2}.m8")
    if not os.path.isfile(m8_file) or os.path.getsize(m8_file) == 0:
        input_dir = os.path.join(outputdir, proteome_id1)
        db_file = os.path.join(outputdir, f"{proteome_id2}_db")
        tmp_dir = os.path.join(outputdir, "tmp")
        os.makedirs(tmp_dir, exist_ok=True)
        
        output_format = "query,target,qlen,tlen,fident,alnlen,mismatch,qstart,qend,tstart,tend,evalue,bits"

        cmd = f"{foldseek} easy-search --format-output {output_format} {input_dir} {db_file} {m8_file} {tmp_dir}"
        run(cmd)

    return m8_file


def get_relevant_matches(outputdir, matches):
    all_protein_matches = dict()
    if os.path.isfile(matches) and os.path.getsize(matches) > 0:
        match = os.path.basename(matches)

        species = match.split('.')[0]
        output_best_hits = os.path.join(outputdir, f"{species}_best_hits.csv")
        output_all_hits = os.path.join(outputdir, f"{species}_all_dom_hits.csv")
        try:
            os.remove(output_all_hits)
        except FileNotFoundError:
            pass
        with open(matches, "r") as f, open(output_best_hits, "w") as fout:
            af = ""
            list_matches = []
            count_best_hits = 0
            same_bits_score = 0

            for line in f.readlines():
                line = line.strip()
                values = line.split()
                aftmp = values[0][:-4]

                if not af:
                    af = aftmp
                if aftmp != af:
                    protein_matches, text, count_best, same_bits = get_best_match(output_all_hits, list_matches)
                    count_best_hits += count_best
                    same_bits_score += same_bits
                    all_protein_matches |= protein_matches
                    fout.write(text)

                    af = aftmp
                    list_matches = []

                list_matches.append(values)

            protein_matches, text, count_best, same_bits = get_best_match(output_all_hits, list_matches)
            count_best_hits += count_best
            same_bits_score += same_bits
            all_protein_matches |= protein_matches

            fout.write(text)
    
    print(f"best hits {species}: {count_best_hits}")
    print(f"same bit score {species}: {same_bits_score}")
    return all_protein_matches


def get_best_match(output_all_hits, list_matches):
    keep = dict()
    pattern = r"[0-9]+\.pdbnum\.pdb_[A-Za-z0-9]+"
    with open(output_all_hits, "a") as fout:
        for line in list_matches:
            if re.search(pattern, line[1]):  # non existing pdb file => ignore match ## long proteins, how to download all the pdb files generated by AF?
                continue
            qstart = int(line[7])
            qend = int(line[8])
            midpoint = (qend + qstart) / 2

            if not keep:
                keep[f"{qstart}-{qend}"] = [line]
            else:
                # group results by domains
                new = True
                for key in keep.keys():
                    saved_s, saved_e = key.split("-")
                    midpoint2 = (int(saved_s) + int(saved_e)) / 2
                    if (midpoint > int(saved_s) and midpoint < int(saved_e)) or (
                        midpoint2 > qstart and midpoint2 < qend
                    ):  # same domain
                        keep[key].append(line)
                        new = False
                if new:
                    keep[f"{qstart}-{qend}"] = [line]
                    new = False

        prot_prot = dict()
        text = ""
        count_best_hits = 0
        same_bits_score = 0
    
        for key, lines in keep.items():
            keeplines = []
            evalue_saved = ""
            bits_saved = 0
            
            for match in lines:
                evalue = float(match[11])
                prot1 = match[0].split('-')[1]
                prot2 = match[1].split('-')[1]
                bits = int(match[12])

                fout.write(f"{prot1},{prot2},{int(line[7])}-{int(line[8])},{int(line[9])}-{int(line[10])},{evalue},{bits}\n")
                if evalue < 1.0e-04 and not keeplines:
                    evalue_saved = evalue
                    bits_saved = bits
                    keeplines.append(match)
                else:
                    if evalue < 1.0e-04 and evalue <= evalue_saved:
                        if evalue == evalue_saved:

                            if bits > bits_saved: #filter same e-value results using the bitscore
                                keeplines = [match]
                            elif bits == bits_saved:
                                same_bits_score +=1
                                keeplines.append(match)
                        elif evalue < evalue_saved:
                            keeplines = [match]
                        evalue_saved = evalue

            for keepline in keeplines:
                count_best_hits += 1
                # AF-P50993-F1-model_v2.pdb.gz => P50993'
                prot1 = keepline[0].split('-')[1]
                prot2 = keepline[1].split('-')[1]
                vals = ', '.join(keepline[2:])

                if prot1 in prot_prot:
                    if prot2 in prot_prot[prot1]:
                        prot_prot[prot1][prot2].append(vals)
                    else:
                        prot_prot[prot1][prot2] = [vals]
                else:
                    prot_prot[prot1] = {prot2: [vals]}
                text+=f"{prot1},{prot2},{int(keepline[7])}-{int(keepline[8])},{int(keepline[9])}-{int(keepline[10])},{float(keepline[11])},{int(keepline[12])}\n"

    return prot_prot, text, count_best_hits, same_bits_score


def get_orthologs(orthologs_file, matches1, matches2):
    with open(orthologs_file, "w") as fout:
        orthologs = []
        for match1_prot1, match1_prots2 in matches1.items():
            for match1_prot2 in match1_prots2:
                for match2_prot1, match2_prots2 in matches2.items():
                    for match2_prot2 in match2_prots2:
                        if match1_prot1 == match2_prot2 and match1_prot2 == match2_prot1:
                            orthologs.append(f"{match1_prot1}-{match1_prot2}")
                            fout.write(f"{match1_prot1},{match1_prot2},{match1_prots2[match1_prot2]},{match2_prots2[match2_prot2]}\n")
        print(f"Number of orthologs found: {len(orthologs)}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("config", metavar="CONFIG_FILE", help="configuration file")
    parser.add_argument(
        "-r",
        "--runfoldseek",
        help="Start analysis by running Foldseek or not",
        choices=["yes", "no"],
        required=True,
    )
    args = parser.parse_args()

    if not os.path.isfile(args.config):
        parser.error(f"Cannot open '{args.config}': " f"no such file or directory")

    config = ConfigParser()
    config.read(args.config)
    outputdir = config["misc"]["basedir"]
    os.makedirs(outputdir, exist_ok=True)
    script_path = os.getcwd()

    if args.runfoldseek == "yes":
        print("Analysis with Foldseek run")
        link1=config["misc"]["af_link1"]
        link2=config["misc"]["af_link2"]

        get_foldseek(script_path)

        proteome_id1 = get_model(outputdir, link1)
        proteome_id2 = get_model(outputdir, link2)

        os.chdir(script_path)

        outputfile1to2 = run_foldseek(outputdir, script_path, proteome_id1, proteome_id2)
        outputfile2to1 = run_foldseek(outputdir, script_path, proteome_id2, proteome_id1)

        protein_matches1 = get_relevant_matches(outputdir, outputfile1to2)
        protein_matches2 = get_relevant_matches(outputdir, outputfile2to1)

        orthologs_file = os.path.join(outputdir, f"{proteome_id1}_{proteome_id2}.csv")

    else:
        # part ran if m8 file already generated in different directory
        print("Analysis of Foldseek files only")
        os.chdir(script_path)
        outputfile1to2 = config["misc"]["foldseek_file1"]
        outputfile2to1 = config["misc"]["foldseek_file2"]
        csvfile = config["misc"]["csvfile"]

        protein_matches1 = get_relevant_matches(outputdir, outputfile1to2)
        protein_matches2 = get_relevant_matches(outputdir, outputfile2to1)

        orthologs_file = os.path.join(outputdir, csvfile)
        
    get_orthologs(orthologs_file, protein_matches1, protein_matches2)
