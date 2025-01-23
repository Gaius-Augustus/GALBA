#!/usr/bin/env python3
"""
startAlign.py

This script splits a genome file into single sequences (or sequence parts) and a protein file
according to contigIDs, then runs alignment programs (exonerate, spaln, or gth) for each
contigID/sequence. If position/list files aren't provided, the entire protein file is used.
After the alignments, the script adjusts coordinate outputs and combines them into a final file.

Ported to Python 3 from the original Perl script: startAlign.pl
Authors of original Perl script: Katharina Hoff & Simone Lange

Usage examples:
  startAlign.py --genome genome.fa --prot db.fa --prg gth
  startAlign.py --genome genome.fa --prot db.fa --prg exonerate --list BLAST.hit.list --pos dna.pos
  startAlign.py --help

Required arguments:
  --genome=FILE    Fasta file with DNA sequences
  --prot=FILE      Fasta file with protein sequences
  --prg=STRING     Alignment program: 'gth', 'exonerate', or 'spaln'

Optional arguments:
  --help           Print this help message
  --CPU=n          Number of CPUs (parallel processes) to use (default: 1)
  --dir=DIR        Working directory path (default: current directory)
  --list=FILE      File with contigID proteinID pairs
  --pos=FILE       File with contig positions: contigID nr_of_prots_mapped start end strand chrID
  --log=FILE       Log file (default: startAlign_<prg>.log in working dir)
  --maxintron=n    For exonerate: maximal intron length (default: 30000)
  --reg            Use region-based alignment (split by positions), else entire contigs
  --offset=n       Number of bp to add before/after cutout coordinates (default: 10000)
  --ALIGNMENT_TOOL_PATH=PATH  Directory containing exonerate/spaln/gth binary
  --args=STRING    Additional command line parameters for the alignment tool
  --nice           Execute alignment calls with 'nice'
"""

import argparse
import os
import sys
import shutil
import subprocess
import tempfile
from concurrent.futures import ProcessPoolExecutor, as_completed
from textwrap import dedent
import time

VERSION = "1.0"


def main():
    parser = argparse.ArgumentParser(add_help=False)
    parser.add_argument("--genome", help="Fasta file with DNA sequences.")
    parser.add_argument("--prot", help="Fasta file with protein sequences.")
    parser.add_argument("--prg", choices=["exonerate", "spaln", "gth"],
                        help="Alignment program: 'exonerate', 'spaln', or 'gth'.")
    parser.add_argument("--CPU", type=int, default=1,
                        help="Number of CPUs (parallel processes). Default 1.")
    parser.add_argument("--dir", help="Working directory. Default is current directory.")
    parser.add_argument("--list", help="File listing contigID proteinID pairs.")
    parser.add_argument("--pos", help="File with contig positions.")
    parser.add_argument("--log", help="Log file (default: startAlign_<prg>.log in working dir).")
    parser.add_argument("--maxintron", type=int, default=30000,
                        help="Max intron length for exonerate (default: 30000).")
    parser.add_argument("--offset", type=int, default=10000,
                        help="Number of bp to add before/after coordinates. Default 10000.")
    parser.add_argument("--reg", action="store_true",
                        help="Use region parts, not entire contigs (like Perl's --reg).")
    parser.add_argument("--help", action="store_true", help="Print this help message and exit.")
    parser.add_argument("--nice", action="store_true", help="Run alignment calls with 'nice'.")
    parser.add_argument("--ALIGNMENT_TOOL_PATH",
                        help="Path to alignment tool binary (exonerate/spaln/gth).")
    parser.add_argument("--args",
                        help="Additional command line parameters for the alignment tool.")
    args = parser.parse_args()

    if args.help:
        print(__doc__)
        parser.print_help(sys.stdout)
        sys.exit(0)

    # Basic checks
    if not args.genome or not args.prot or not args.prg:
        sys.stderr.write("ERROR: --genome, --prot, and --prg are required.\n")
        parser.print_help(sys.stderr)
        sys.exit(1)

    # Determine working directory
    if not args.dir:
        work_dir = os.getcwd()
    else:
        work_dir = os.path.abspath(args.dir)

    # Log file
    if not args.log:
        args.log = os.path.join(work_dir, f"startAlign_{args.prg}.log")

    # Open the log
    try:
        log_handle = open(args.log, "w")
    except OSError:
        sys.stderr.write(f"ERROR: Could not open log file {args.log}\n")
        sys.exit(1)

    # Write function to log
    def log(msg):
        now = time.strftime("%Y-%m-%d %H:%M:%S")
        log_handle.write(f"# {now}: {msg}\n")

    # If no list file or pos file => use entire protein file
    use_whole_prot = False
    if not args.list or not args.pos:
        use_whole_prot = True

    # Check existence of files
    if not os.path.exists(args.genome):
        sys.stderr.write(f"ERROR: Genome file {args.genome} does not exist.\n")
        sys.exit(1)

    if not os.path.exists(args.prot):
        sys.stderr.write(f"ERROR: Protein file {args.prot} does not exist.\n")
        sys.exit(1)

    pos_file_abs = None
    list_file_abs = None
    if not use_whole_prot:
        if not os.path.exists(args.pos):
            sys.stderr.write(f"ERROR: Position file {args.pos} does not exist.\n")
            sys.exit(1)
        else:
            pos_file_abs = os.path.abspath(args.pos)

        if not os.path.exists(args.list):
            sys.stderr.write(f"ERROR: List file {args.list} does not exist.\n")
            sys.exit(1)
        else:
            list_file_abs = os.path.abspath(args.list)

    # Check environment for spaln
    if args.prg == "spaln":
        if "ALN_DBS" not in os.environ:
            sys.stderr.write(
                "ERROR: The environment variable ALN_DBS for spaln2 is not defined.\n"
                "       Please set e.g. 'export ALN_DBS=/path/to/spaln2/seqdb'.\n"
            )
            sys.exit(1)
        if "ALN_TAB" not in os.environ:
            sys.stderr.write(
                "ERROR: The environment variable ALN_TAB for spaln2 is not defined.\n"
                "       Please set e.g. 'export ALN_TAB=/path/to/spaln2/table'.\n"
            )
            sys.exit(1)

    # Adjust for "spalnErrAdj" in the Perl script
    spalnErrAdj = 80 if args.prg == "spaln" else 0

    # Create or clean temporary folder
    tmp_dir = os.path.join(work_dir, f"tmp_{args.prg}")
    if os.path.isdir(tmp_dir):
        log(f"Removing existing temp directory {tmp_dir}")
        shutil.rmtree(tmp_dir)
    os.makedirs(tmp_dir, exist_ok=True)
    log(f"Created temp directory {tmp_dir}")

    # Create or clean align output directory
    align_dir = os.path.join(work_dir, f"align_{args.prg}")
    if os.path.isdir(align_dir):
        log(f"Removing existing alignment directory {align_dir}")
        shutil.rmtree(align_dir)
    os.makedirs(align_dir, exist_ok=True)
    log(f"Created alignment directory {align_dir}")

    # Link genome_file and protein_file into tmp_dir (avoid multiple dots for spaln)
    genome_file_abs = os.path.abspath(args.genome)
    genome_link_name = link_in_tmp(genome_file_abs, tmp_dir, log)
    prot_file_abs = os.path.abspath(args.prot)
    prot_link_name = link_in_tmp(prot_file_abs, tmp_dir, log)

    # We'll define the names the same as the Perl script does
    # The "prot_file_base" is the name of the protein file without the path:
    prot_file_base = os.path.basename(prot_link_name)
    # The prot_addstop_file is like "X_addstop.Y" if X.Y was the original
    prot_addstop_file = None
    split_basename = prot_file_base.split(".")
    if len(split_basename) == 2:
        prot_addstop_file = f"{split_basename[0]}_addstop.{split_basename[1]}"
    else:
        sys.stderr.write(
            "ERROR: Could not define prot_addstop_file name properly.\n"
            "       Protein file name has unexpected format.\n"
        )
        sys.exit(1)

    # "contigIDs" => { contigID -> {"start": x, "end": y, "seq": someName} }
    # "protIDs" => { protID -> [ contigIDs... ] } or [ seqNames... ]
    contigIDs = {}
    protIDs = {}

    # If not using whole protein, read pos and list
    if not use_whole_prot:
        read_pos_file(pos_file_abs, contigIDs, log)
        read_list_file(list_file_abs, protIDs, contigIDs, args.reg, log)

    # Next, read protein sequences and write them out
    # each protein may map to multiple contigIDs or seqNames
    # "counterW" is used in Perl to log "added stop amino acid"
    # but we replicate the same logic in code
    # We'll store how many times we wrote a protein with stop.
    counterW = [0]  # keep inside a list so we can mutate in function

    if not write_proteins(
        prot_link_name,
        tmp_dir,
        prot_addstop_file,
        use_whole_prot,
        protIDs,
        contigIDs,
        args.reg,
        args.prg,
        log,
        counterW
    ):
        log("ERROR while writing proteins.")
        sys.exit(1)

    # If CPU>1 or (CPU=1 and not protWhole), we read genome sequences from genome_link_name
    seq_dict = {}
    if args.CPU > 1 or (args.CPU == 1 and not use_whole_prot):
        read_genome_into_dict(genome_link_name, seq_dict, log)

    # Now we call start_align in the Perl script
    # We'll replicate that in a function:
    start_align(
        tmp_dir,
        align_dir,
        seq_dict,
        args,
        contigIDs,
        protIDs,
        prot_addstop_file,
        spalnErrAdj,
        log
    )

    # Clean up empty files, remove tmp_dir
    clean_up_empty_files(align_dir, log)
    log(f"Removing temp directory {tmp_dir}")
    shutil.rmtree(tmp_dir, ignore_errors=True)

    log_handle.close()


def link_in_tmp(original_path, tmp_dir, log_fn):
    """
    Create a symlink of 'original_path' into 'tmp_dir'.
    If the file has multiple dots, try to rename it so that for Spaln it has at most one dot.
    Return the link name path.
    """
    base_name = os.path.basename(original_path)
    # If more than one dot in the file name, replace extra ones with '_'
    # e.g. "my.file.fa" => "my_file.fa"
    dot_count = base_name.count(".")
    if dot_count > 1:
        # replace all but the last dot
        parts = base_name.rsplit(".", 1)
        left = parts[0].replace(".", "_")
        base_name = f"{left}.{parts[1]}"

    link_name = os.path.join(tmp_dir, base_name)
    if not os.path.exists(link_name):
        cmd = f"ln -s {original_path} {link_name}"
        log_fn(f"Creating symlink: {cmd}")
        ret = os.system(cmd)
        if ret != 0:
            sys.stderr.write(f"ERROR: failed to run: {cmd}\n")
            sys.exit(1)
    else:
        log_fn(f"Symlink {link_name} already exists. Assuming it's valid.")
    return link_name


def read_pos_file(pos_file, contigIDs, log_fn):
    """Read the --pos file, fill contigIDs dict: { contigID -> { start, end, seq } }."""
    log_fn(f"Reading position data from {pos_file}")
    with open(pos_file, "r") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            parts = line.split()
            # contigID nr_of_prots_mapped start end strand chrID
            if len(parts) == 6:
                contig = parts[0]
                start_ = int(parts[2]) - 1
                end_ = int(parts[3]) - 1
                seq_name = parts[5]
                contigIDs[contig] = {
                    "start": start_,
                    "end": end_,
                    "seq": seq_name
                }


def read_list_file(list_file, protIDs, contigIDs, use_reg, log_fn):
    """Read the --list file, fill protIDs dict: { protID -> [ contigIDs or seqNames ] }."""
    log_fn(f"Reading list data from {list_file}")
    with open(list_file, "r") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            parts = line.split()
            if len(parts) == 2:
                contigID = parts[0]
                protID = parts[1]
                if use_reg:
                    # store each contigID
                    protIDs.setdefault(protID, []).append(contigID)
                else:
                    # store seqName from contigIDs
                    if contigID in contigIDs and "seq" in contigIDs[contigID]:
                        seqname = contigIDs[contigID]["seq"]
                        # avoid duplicates
                        if seqname not in protIDs.setdefault(protID, []):
                            protIDs[protID].append(seqname)


def read_genome_into_dict(genome_file, seq_dict, log_fn):
    """Read entire genome fasta into a dictionary: { seqname: sequenceString }."""
    log_fn(f"Reading genome data from {genome_file}")
    with open(genome_file, "r") as f:
        data = f.read().strip().split(">")
    for chunk in data:
        if not chunk.strip():
            continue
        lines = chunk.split("\n", 1)
        header = lines[0].strip()
        seqname = header.split()[0]
        sequence = ""
        if len(lines) > 1:
            sequence = lines[1].replace("\n", "")
        seq_dict[seqname] = sequence


def write_proteins(
    prot_file, tmp_dir, prot_addstop_file,
    use_whole_prot, protIDs, contigIDs, reg, prg, log_fn, counterW
):
    """
    Read protein sequences from 'prot_file' and write them to either:
      1) one single file 'prot_addstop_file' if use_whole_prot is True
      2) multiple smaller files if use_whole_prot is False
    If prg == 'gth', append '*' to the end of each protein sequence.
    'counterW' is a mutable list so we can increment it to track logging the addition of '*'.
    """
    try:
        pf = open(prot_file, "r")
    except OSError:
        return False

    log_fn(f"Reading protein data from {prot_file}")

    seqname = None
    seqbuf = []
    first_line = True

    def flush_protein(seq_name, seq_data):
        # flush one protein record
        if not seq_name:
            return
        joined_seq = "".join(seq_data)
        if use_whole_prot:
            # single output: tmp_dir/prot_addstop_file
            # "counterW[0]" tracks if we appended '*'
            append_prot(
                os.path.join(tmp_dir, prot_addstop_file),
                seq_name,
                joined_seq,
                prg,
                log_fn,
                counterW
            )
        else:
            # For each contig (region or entire seq) that references seq_name
            if seq_name in protIDs:
                for cID in protIDs[seq_name]:
                    if reg:
                        # region mode -> cID is contig ID, but we only write a file if contigIDs have an entry
                        if cID in contigIDs:
                            # e.g. "$tmpDir/@{ $protIDs{ $seqname }}[$i].fa"
                            out_file = os.path.join(tmp_dir, f"{cID}.fa")
                            append_prot(out_file, seq_name, joined_seq, prg, log_fn, counterW)
                    else:
                        # not region mode => cID is actually the "seq" from contigIDs
                        out_file = os.path.join(tmp_dir, f"prot{cID}.fa")
                        append_prot(out_file, seq_name, joined_seq, prg, log_fn, counterW)

    for line in pf:
        line = line.rstrip("\n")
        if line.startswith(">"):
            # flush previous
            if not first_line:
                flush_protein(seqname, seqbuf)
            first_line = False
            seqbuf = []
            seqname = line[1:].strip()
        else:
            seqbuf.append(line)

    # flush the last one
    flush_protein(seqname, seqbuf)
    pf.close()
    return True


def append_prot(outfile, seq_name, seq_data, prg, log_fn, counterW):
    """
    Write a single protein entry to 'outfile'. If prg=='gth', append '*' to the sequence.
    Log the addition the first time we do it.
    """
    mode = "a" if os.path.exists(outfile) else "w"
    added_stop = (prg == "gth")
    with open(outfile, mode) as out:
        out.write(f">{seq_name}\n")
        if added_stop:
            out.write(seq_data + "*\n")
            if counterW[0] == 0:
                log_fn("Adding stop amino acid ('*') to protein sequences")
            counterW[0] += 1
        else:
            out.write(seq_data + "\n")


def start_align(tmp_dir, align_dir, seq_dict, args,
                contigIDs, protIDs, prot_addstop_file,
                spalnErrAdj, log_fn):
    """
    Core alignment logic. Creates tasks per contig or per sequence, calls the
    alignment tool in parallel if CPU>1, merges output into 'prgsrc.concat.aln',
    and adjusts coordinate files if needed.
    """
    prgsrc = args.prg
    whole_pred_file = os.path.join(align_dir, f"{prgsrc}.concat.aln")
    # Remove if exists
    if os.path.exists(whole_pred_file):
        os.remove(whole_pred_file)

    # We replicate the parallel approach from Perl with Parallel::ForkManager
    # We'll define tasks, run them with ProcessPoolExecutor
    tasks = []
    results = []
    n_tasks = 0

    # We'll define a helper that does the alignment calls
    # Return the path to the final alignment chunk to be appended
    def align_one_target(
        genome_file, protein_file, out_file, err_file,
        prg, offset, CPU_local
    ):
        # Actually call the alignment programs
        if prg == "exonerate":
            call_exonerate(genome_file, protein_file, out_file, err_file,
                           args.maxintron, args.nice, args.ALIGNMENT_TOOL_PATH, args.args)
            # Adjust exonerate output
            if not args.reg:
                # If the entire sequence is used (no offset shifting for each contig),
                # we skip coordinate shifting. The original script does coordinate
                # shifting only for region-based approach. But let's carefully check:
                # The Perl script does "adjust_exonerate" only if $reg=1. Indeed, it does that in sub adjust_exonerate($outfile, $ID).
                # We'll replicate that logic. If reg=1 => do it, else skip.
                return out_file
            else:
                # We don't know the contig ID here. The original script passes $ID.
                # We'll handle that in the function that calls align_one_target,
                # so that it can pass a contigID for exonerate adjusting.
                return out_file
        elif prg == "spaln":
            call_spaln(genome_file, protein_file, out_file, err_file, CPU_local,
                       args.nice, args.ALIGNMENT_TOOL_PATH, args.args)
            return out_file
        elif prg == "gth":
            call_gth(genome_file, protein_file, out_file, err_file,
                     args.nice, args.ALIGNMENT_TOOL_PATH, args.args)
            return out_file
        else:
            return None

    executor = ProcessPoolExecutor(max_workers=args.CPU)

    # If reg=1 => for each ID in contigIDs => create tasks
    # else => for each seq in seq_dict => create tasks
    # or handle the single-case if CPU=1 but we have a large genome, etc.
    # We'll replicate the Perl logic as closely as possible.
    if args.reg:
        # region-based approach
        for contigID, cdata in contigIDs.items():
            # Create target file if not present
            target_file = os.path.join(
                tmp_dir, f"{cdata['seq']}{contigID}.fa"
            )  # e.g. chr1contig123.fa
            if not os.path.exists(target_file):
                length_ = cdata["end"] - cdata["start"] + 1 + (2 * args.offset)
                subseq_start = cdata["start"] - args.offset
                if subseq_start < 0:
                    subseq_start = 0
                if cdata["seq"] in seq_dict:
                    subseq = seq_dict[cdata["seq"]][subseq_start: subseq_start + length_]
                else:
                    subseq = ""

                # Write to file
                write_fasta(target_file, cdata["seq"], subseq)

            # The protein file for this contig is something like f"{contigID}.fa"
            protein_file = os.path.join(tmp_dir, f"{contigID}.fa")
            if not os.path.exists(protein_file):
                # If the user never had a matching protein for that contig, it might not exist
                # Then skip alignment
                continue

            err_file = os.path.join(align_dir, f"{cdata['seq']}.{contigID}.{args.prg}.stderr")
            out_file = os.path.join(align_dir, f"{cdata['seq']}.{contigID}.{args.prg}.aln")

            # We'll define a function wrapper for the exonerate coordinate shift,
            # because the Perl script calls "adjust_exonerate($outfile, $ID)" with $ID=contigID
            def run_task():
                partial = align_one_target(
                    target_file, protein_file, out_file, err_file,
                    args.prg, args.offset, 1  # parallel=1 for the actual tool?
                )
                # If exonerate => adjust coords only if reg=1
                # If prg=spaln/gth => run "adjust" sub. Then cat partial => concat file => remove partial
                # We'll replicate: sub adjustExonerate => shift coords by start-offset
                if args.prg == "exonerate":
                    adj = adjust_exonerate(
                        partial, contigID, contigIDs, args.offset
                    )
                    os.remove(partial)
                    return adj
                elif args.prg in ["spaln", "gth"]:
                    # sub adjust => shift coords by start-offset + spalnErrAdj
                    adj = adjust_spaln_gth(
                        partial, contigID, contigIDs, args.offset, spalnErrAdj
                    )
                    os.remove(partial)
                    return adj
                else:
                    return partial

            future = executor.submit(run_task)
            tasks.append((future, out_file))
            n_tasks += 1

    else:
        # no region => entire genome. If CPU>1 or (CPU=1 and not protWhole)
        # we do one alignment per seq in seq_dict. Else we do a single call with the entire genome_file
        if args.CPU > 1 or (args.CPU == 1 and not is_whole_prot(protIDs)):
            # For each seq in seq_dict => a task
            for seqname in seq_dict.keys():
                target_file = os.path.join(tmp_dir, f"genome{seqname}.fa")
                if not os.path.exists(target_file):
                    # write the entire chromosome as a single file
                    write_fasta(target_file, seqname, seq_dict[seqname])

                # The protein file is "prot{seqname}.fa" if not use_whole_prot
                # or if using the entire DB => "prot_addstop_file"
                if not is_whole_prot(protIDs):
                    protein_file = os.path.join(tmp_dir, f"prot{seqname}.fa")
                    if not os.path.exists(protein_file):
                        # skip if no relevant proteins
                        continue
                else:
                    protein_file = os.path.join(tmp_dir, prot_addstop_file)

                err_file = os.path.join(align_dir, f"{seqname}.{args.prg}.stderr")
                out_file = os.path.join(align_dir, f"{seqname}.{args.prg}.aln")

                def run_task():
                    partial = align_one_target(
                        target_file, protein_file, out_file, err_file,
                        args.prg, args.offset, args.CPU if args.prg == "spaln" else 1
                    )
                    # For exonerate/spaln/gth in no region mode, the Perl script does no coordinate shift
                    return partial  # we just return the alignment
                future = executor.submit(run_task)
                tasks.append((future, out_file))
                n_tasks += 1
        else:
            # Single call with the entire genome_file + entire prot DB
            # exonerate => calls with entire genome_file, entire prot DB
            # Then cat => prgsrc.concat.aln
            # We'll do it outside parallel
            # The original script does call_gth(...) or call_exonerate(...) or call_spaln(...)
            err_file = os.path.join(align_dir, f"{args.prg}.stderr")
            out_file = os.path.join(align_dir, f"{args.prg}.aln")

            partial = align_one_target(
                genome_link_name,
                os.path.join(tmp_dir, prot_addstop_file),
                out_file, err_file,
                args.prg, args.offset, 1
            )
            # Then cat out_file => prgsrc.concat.aln, remove out_file
            with open(os.path.join(align_dir, f"{args.prg}.concat.aln"), "ab") as wpf:
                with open(partial, "rb") as rf:
                    wpf.write(rf.read())
            os.remove(partial)
            return  # done

    # Wait for parallel tasks to finish
    with open(whole_pred_file, "ab") as wpf:
        for future, out_file in as_completed(dict(tasks)):
            adj_file = future.result()
            # cat adj_file => whole_pred_file
            # remove adj_file
            with open(adj_file, "rb") as rf:
                wpf.write(rf.read())
            os.remove(adj_file)


def call_exonerate(
    genome_file, prot_file, stdout_file, stderr_file,
    max_intron, use_nice, tool_path, extra_args
):
    cmd_list = []
    if use_nice:
        cmd_list.append("nice")

    if tool_path:
        exo_bin = os.path.join(tool_path, "exonerate")
    else:
        exo_bin = "exonerate"

    cmd_list.extend([
        exo_bin,
        "--model", "protein2genome",
        f"--maxintron", str(max_intron),
        "--showtargetgff",
        "--showalignment", "no",
        "--query", prot_file,
        "--target", genome_file
    ])
    if extra_args:
        cmd_list.extend(extra_args.split())

    cmd_str = " ".join(cmd_list) + f" > {stdout_file} 2> {stderr_file}"
    run_command(cmd_str)


def call_gth(genome_file, prot_file, stdout_file, stderr_file,
             use_nice, tool_path, extra_args):
    # We replicate the logic: create a subfolder for this combination, link the protein file
    # but the parallel approach might cause collisions if done similarly. We'll mimic it carefully.
    # The script does: "mkdir <protFileFor_...>", "ln -s ../<prot_file>", call gth, then rm -r
    # We'll do that in a temp directory unique to the genome_file's name.
    base_g = os.path.basename(genome_file)
    safe_subdir = os.path.join(os.path.dirname(prot_file), f"protFileFor_{base_g}_{os.getpid()}")
    os.makedirs(safe_subdir, exist_ok=True)
    link_basename = os.path.basename(prot_file)
    link_path = os.path.join(safe_subdir, link_basename)
    if not os.path.exists(link_path):
        # e.g. ln -s ../<prot_file> <link_path>
        cmd = f"ln -s ../{link_basename} {link_path}"
        # but we have to be sure it's from the perspective of the subdir
        os.chdir(os.path.dirname(prot_file))
        run_command(cmd)
        os.chdir(safe_subdir)

    # Build GTH command
    cmd_list = []
    if use_nice:
        cmd_list.append("nice")

    if tool_path:
        gth_bin = os.path.join(tool_path, "gth")
    else:
        gth_bin = "gth"

    if not extra_args:
        # default in the Perl script:
        # -prseedlength 20 -prminmatchlen 20 -prhdist 2
        extra_args = "-prseedlength 20 -prminmatchlen 20 -prhdist 2"

    # Hard-coded from the Perl script:
    # -gcmincoverage 80, plus extra_args, plus -skipalignmentout -paralogs -gff3out
    cmd_list.extend([
        gth_bin,
        "-genomic", genome_file,
        "-protein", os.path.join(safe_subdir, link_basename),
        "-gff3out",
        "-skipalignmentout",
        "-paralogs",
        "-gcmincoverage", "80"
    ])
    # add user-specified extra_args
    cmd_list.extend(extra_args.split())
    cmd_list.extend(["-o", stdout_file])
    final_cmd = " ".join(cmd_list) + f" 2> {stderr_file}"

    # run
    run_command(final_cmd)

    # remove the subdir
    os.chdir(os.path.dirname(prot_file))  # return to the tmp_dir
    shutil.rmtree(safe_subdir, ignore_errors=True)


def call_spaln(genome_file, prot_file, stdout_file, stderr_file,
               cpus, use_nice, tool_path, extra_args):
    # For Spaln, we might need to run 'makdbs -KP' / 'perl makblk.pl ...' if .idx or .bkp do not exist
    # We'll replicate the logic. We'll store them next to the .fa in the same directory.
    # The script references $split[0].idx, so let's do similar.
    # We'll define a helper function check_and_create_spaln_db(...) below
    check_and_create_spaln_db(genome_file, "P", use_nice, tool_path)
    check_and_create_spaln_blk(genome_file, "P", use_nice)

    check_and_create_spaln_db(prot_file, "A", use_nice, tool_path)
    check_and_create_spaln_blk(prot_file, "A", use_nice)

    cmd_list = []
    if use_nice:
        cmd_list.append("nice")
    if tool_path:
        # we also ensure that environment PATH includes tool_path for 'makdbs', but let's do directly:
        spaln_bin = os.path.join(tool_path, "spaln")
    else:
        spaln_bin = "spaln"

    cmd_list.append(spaln_bin)
    # from the Perl script:
    #   spaln -Q7 -O0 [ -t[$cpus] ] $tFile $qFile > $stdoutfile 2>$errorfile
    if cpus > 1:
        cmd_list.append(f"-t[{cpus}]")
    cmd_list.extend(["-Q7", "-O0"])
    if extra_args:
        cmd_list.extend(extra_args.split())

    cmd_list.extend([genome_file, prot_file])

    final_cmd = " ".join(cmd_list) + f" > {stdout_file} 2> {stderr_file}"
    run_command(final_cmd)


def check_and_create_spaln_db(fasta_path, mode, use_nice, tool_path):
    """
    If .idx for the fasta doesn't exist, run 'makdbs -K<mode> fasta_path'.
    """
    # e.g. "mycontig.fa" => "mycontig" => check if "mycontig.idx" exists
    base, ext = os.path.splitext(fasta_path)
    idxfile = f"{base}.idx"
    if not os.path.isfile(idxfile):
        cmd_list = []
        if use_nice:
            cmd_list.append("nice")
        if tool_path:
            mkdbs_bin = os.path.join(tool_path, "makdbs")
        else:
            mkdbs_bin = "makdbs"
        cmd_list.extend([mkdbs_bin, f"-K{mode}", fasta_path])
        run_command(" ".join(cmd_list))


def check_and_create_spaln_blk(fasta_path, mode, use_nice):
    """
    If .bkp or .bka for the fasta doesn't exist, run
    'perl $ENV{'ALN_DBS'}/makblk.pl -Wxxx.bkp -K<mode> fasta_path'.
    """
    # e.g. "mycontig.fa" => "mycontig.bkp"
    base, ext = os.path.splitext(fasta_path)
    ext_map = {"P": "bkp", "A": "bka"}
    if mode not in ext_map:
        return
    bfile = f"{base}.{ext_map[mode]}"
    if not os.path.isfile(bfile):
        # e.g. perl $ENV{'ALN_DBS'}/makblk.pl -Wmycontig.bkp -KP mycontig.fa
        cmd_list = []
        if use_nice:
            cmd_list.append("nice")
        perl_script = os.path.join(os.environ["ALN_DBS"], "makblk.pl")
        cmd_list.extend([
            "perl", perl_script,
            f"-W{bfile}",
            f"-K{mode}{base[-1]}",  # but the original uses e.g. -KP or -KA
        ])
        # Actually, the Perl code: `-W$split[0].bkp -KP $tFile`.
        # If mode='P' => "-KP", if mode='A' => "-KA"
        cmd_list[-1] = f"-K{mode} {fasta_path}"
        # We'll do partial approach:
        # fix: we do `-K{mode} {fasta_path}`
        # The script is slightly tricky, let's do it literally:
        # for a genomic file => -KP
        # for a protein file => -KA
        # We'll just do a direct approach:
        cmd_list[-1] = f"-K{mode} {fasta_path}"
        # so the final is: perl <...>/makblk.pl -Wmycontig.bkp -KP mycontig.fa
        run_command(" ".join(cmd_list))


def run_command(cmd):
    """Run a shell command, raise if fails."""
    ret = os.system(cmd)
    if ret != 0:
        raise RuntimeError(f"Failed to execute: {cmd}")


def write_fasta(filepath, header, sequence, line_width=80):
    """Write a FASTA with given header and sequence, wrapping at line_width."""
    with open(filepath, "w") as f:
        f.write(f">{header}\n")
        start = 0
        n = len(sequence)
        while start < n:
            f.write(sequence[start:start+line_width] + "\n")
            start += line_width


def is_whole_prot(protIDs):
    """Return True if there are no specific references in protIDs (i.e., empty => entire DB?)."""
    return len(protIDs) == 0


def adjust_spaln_gth(output_file, contigID, contigIDs, offset, spalnErrAdj):
    """
    Adjust the GFF3 coordinates for spaln/gth by start-offset + spalnErrAdj.
    Return path to the adjusted file (output_file + '.adj').
    """
    adj_file = f"{output_file}.adj"
    start_ = contigIDs[contigID]["start"] - offset
    with open(output_file, "r") as inp, open(adj_file, "w") as outp:
        for line in inp:
            line = line.rstrip("\n")
            if line.startswith("##sequence-region"):
                # e.g. "##sequence-region <id> start end"
                parts = line.split()
                if len(parts) >= 4:
                    try:
                        # the last two are start..end
                        old_start = int(parts[-2])
                        old_end = int(parts[-1])
                        new_start = old_start + start_ + spalnErrAdj
                        new_end = old_end + start_ + spalnErrAdj
                        outp.write(f"##sequence-region\t{parts[1]} {new_start} {new_end}\n")
                    except ValueError:
                        outp.write(line + "\n")
                else:
                    outp.write(line + "\n")
            else:
                f = line.split()
                if len(f) >= 5:
                    # typical GFF line
                    # f[3], f[4] => start, end
                    try:
                        old_s = int(f[3])
                        old_e = int(f[4])
                        new_s = old_s + start_ + spalnErrAdj
                        new_e = old_e + start_ + spalnErrAdj
                        f[3] = str(new_s)
                        f[4] = str(new_e)
                        outp.write("\t".join(f) + "\n")
                    except ValueError:
                        outp.write(line + "\n")
                else:
                    outp.write(line + "\n")
    return adj_file


def adjust_exonerate(output_file, contigID, contigIDs, offset):
    """
    Adjust exonerate output. The Perl script modifies lines containing
    'Target range:', 'vulgar:', or GFF lines. It's complicated layout.
    We'll replicate the logic. Then return the new file name.
    """
    adj_file = f"{output_file}.adj"
    start_ = contigIDs[contigID]["start"] - offset

    with open(output_file, "r") as inp, open(adj_file, "w") as outp:
        add_space = ""
        for line in inp:
            line_nl = line.rstrip("\n")
            if "  Target range: " in line_nl:
                # e.g. "  Target range: 123 -> 456"
                # We'll parse out the two coords
                import re
                m = re.search(r"Target range: (\d+) -> (\d+)", line_nl)
                if m:
                    st = int(m.group(1)) + start_
                    en = int(m.group(2)) + start_
                    diff = len(str(en)) - len(m.group(2))
                    add_space = " " * diff
                    outp.write(f"  Target range: {st} -> {en}\n")
                else:
                    outp.write(line_nl + "\n")
            elif line_nl.strip().startswith(("vulgar:")):
                # 'vulgar:' lines => shift fields 6,7 by start_
                parts = line_nl.split()
                if len(parts) >= 8:
                    # parts[6], parts[7]
                    parts[6] = str(int(parts[6]) + start_)
                    parts[7] = str(int(parts[7]) + start_)
                outp.write(" ".join(parts) + "\n")
            elif "\t" in line_nl:
                # Possibly a GFF line with 8 or 9 fields
                f = line_nl.split("\t")
                if len(f) >= 5:
                    try:
                        s_ = int(f[3]) + start_
                        e_ = int(f[4]) + start_
                        f[3] = str(s_)
                        f[4] = str(e_)
                        outp.write("\t".join(f) + "\n")
                    except ValueError:
                        outp.write(line_nl + "\n")
                else:
                    outp.write(line_nl + "\n")
            else:
                # Possibly lines with the format: "   123 : something : 345"
                # or lines with weird alignment
                # The Perl script tries to preserve spacing. We'll do partial.
                # We'll adopt the minimal approach: if there's a pattern: ^(\s+)(\d+) : ([^:]+) : (\d+)
                import re
                m = re.match(r"^(\s+)(\d+) : ([^:]+) : (\d+)", line_nl)
                if m:
                    s_ = int(m.group(2)) + start_
                    e_ = int(m.group(4)) + start_
                    prefix = m.group(1)
                    middle = m.group(3)
                    outp.write(f"{prefix}{s_} : {middle} : {e_}\n")
                else:
                    # lines with ^(\s+)(\d+) : ([a-zA-Z\.\-]+) or such
                    # The original script tries to preserve spacing. We'll just replicate it as is.
                    # The Perl script has some logic to add add_space. We'll do a minimal approach:
                    if line_nl.strip().startswith(tuple(["|", ".", "!", ":"])):
                        outp.write(add_space + line_nl + "\n")
                    else:
                        outp.write(line_nl + "\n")

    return adj_file


def clean_up_empty_files(align_dir, log_fn):
    """
    Find empty files in align_dir, remove them. Then log the step.
    """
    log_fn("Deleting empty files in alignment directory.")
    for root, dirs, files in os.walk(align_dir):
        for f in files:
            path = os.path.join(root, f)
            if os.path.getsize(path) == 0:
                log_fn(f"Removing empty file: {path}")
                os.remove(path)


if __name__ == "__main__":
    main()
