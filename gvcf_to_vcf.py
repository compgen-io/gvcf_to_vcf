#!/usr/bin/env python3

import os
import sys
import gzip

try:
    import pysam
except:
    sys.stderr.write("Missing pysam library!\n")
    sys.exit(1)


def help(msg=None):
    if msg:
        sys.stderr.write("ERROR: %s\n" % msg)
    else:
        sys.stderr.write("""
This program will convert a single-patient GVCF file to a VCF file. This program will 
re-call GT values using the below options. For positions that are part of a GVCF block, 
by default any homozygous ref call will also be outputted, provided it meets the minimum 
depth (MIN_DP) threshold. Any FORMAT fields that are based on the number of alt calls 
(type A or R) will be truncated to remove <NON_REF> values. Any FORMAT field with a type "G"
will be removed.

The GVCF file is assumed to contain at least the following FORMAT fields:
  DP     - read depth
  AD     - depth for ref and alt alleles
  MIN_DP - minumum read depth across a GVCF block
                         
Additionally, the END INFO field is also expected to be present. For reference calls that are
part of GVCF blocks (contiguous blocks of homozygous ref calls), AD and DP values will be
inferred from the MIN_DP value for the block.
                         
""")
    sys.stderr.write("""
Usage: gvcf_to_vcf.py {opts} input.gvcf genome.fa

  Options:
    --min-dp val            Minimum DP to report (default: 10)

    --recalc-gt             Recalculate the GT field based on AD values (see below)
    --min-gt-count val      Minimum reads to support a GT call (default: 1)
    --min-gt-af val         Minumum allele-frequency for a GT call (default: 0.1)
                     
    --skip-ambiguous        Skip writing reference bases that are "N"
    --only-alt              Only write out positions with a variant call
    --only-passing          Only write out positions that meet minumum DP threshold
                            (default: write out lines with FILTERs set)
                     
""")
    sys.exit(1)


def convert_gvcf(fname, ref_fname, min_dp, recalc_gt, min_gt_count, min_gt_af, only_passing, only_alt, skip_ambiguous):
    
    if fname == '-':
        f = sys.stdin
    else:
        f = open(fname, 'rb')
        magic = f.read(2)
        f.close()
        if magic[0] == 0x1f and magic[1] == 0x8b:
            f = gzip.open(fname, 'rt')
        else:
            f = open(fname, 'rt')

    fasta = pysam.FastaFile(ref_fname)

    format_keys_to_remove = []
    format_keys_to_adj = []
    info_keys_to_adj = []
    found_min_dp = False
    found_dp = False
    found_ad = False

    last_chrom = None
    chrom_fail = False

    for line in f:
        if line[0] == '#':
            # TODO: add checks for presence of required FORMAT/INFO fields
            if line[:8] == '##INFO=<':
                tmp = line[8:].strip()
                if tmp[-1] == '>':
                    tmp = tmp[:-1]

                id = ""
                num = ""
                for el in csv_split(tmp):
                    kv = el.split('=')
                    if len(kv) == 2:
                        k = kv[0]
                        v = kv[1]
                        if k == 'ID':
                            id = v
                        elif k == 'Number':
                            num = v
                if num.upper() in ["A","."]:
                    info_keys_to_adj.append(id)

            elif line[:10] == '##FORMAT=<':
                tmp = line[10:].strip()
                if tmp[-1] == '>':
                    tmp = tmp[:-1]

                id = ""
                num = ""

                for el in csv_split(tmp):
                    kv = el.split('=')
                    if len(kv) == 2:
                        k = kv[0]
                        v = kv[1]
                        if k == 'ID':
                            id = v
                        elif k == 'Number':
                            num = v

                if id == 'MIN_DP':
                    found_min_dp = True
                if id == 'DP':
                    found_dp = True
                if id == 'AD':
                    found_ad = True
                
                if num.upper() in ["A","."]:
                    format_keys_to_adj.append(id)
                if num.upper() == "G":
                    format_keys_to_remove.append(id)

            elif line[:10] == '##contig=<':
                tmp = line[10:].strip()
                if tmp[-1] == '>':
                    tmp = tmp[:-1]

                id = ""

                for el in csv_split(tmp):
                    kv = el.split('=')
                    if len(kv) == 2:
                        k = kv[0]
                        v = kv[1]
                        if k == 'ID':
                            id = v
                
                if id and id not in fasta.references:
                    # don't write out CONTIG lines that aren't present in the FASTA ref
                    continue



            if line[:6] == '#CHROM':
                sys.stdout.write('##FILTER=<ID=LowDepth,Description="Low sequencing depth">\n')
                sys.stdout.write('##INFO=<ID=DELETION,Number=0,Type=Flag,Description="Is the position part of a deleted region?">\n')

                if not found_min_dp or not found_dp or not found_ad:
                    sys.stderr.write("Missing a required FORMAT element!\n")
                    sys.exit(1)

                sys.stdout.write("##%s\n" % ",".join(info_keys_to_adj))
                sys.stdout.write("##%s\n" % ",".join(format_keys_to_adj))

            sys.stdout.write(line)
            continue

        cols = line.strip('\n').split('\t')

        chrom = cols[0]
        pos = int(cols[1])
        id = cols[2]
        refbase = cols[3]
        alts = cols[4]
        qual = cols[5]
        filter = cols[6]
        info = cols[7]
        format = cols[8]
        format_vals = cols[9]

        format_kv = {}
        info_kv = {}

        # check chrom is in fasta
        if chrom != last_chrom:
            last_chrom = chrom
            if chrom in fasta.references:
                chrom_fail = False
            else:
                sys.stderr.write("Skipping ref: %s (missing from FASTA)\n" % chrom)
                chrom_fail = True
        
        if chrom_fail:
            continue

        gvcf_block_end = -1

        for val in info.split(';'):
            if '=' in val:
                k, v = val.split('=')
                if k == 'END':
                    gvcf_block_end = int(v)

        for val in info.split(';'):
            if '=' in val:
                kv = val.split('=',1)
                info_kv[kv[0]] = kv[1]
            else:
                info_kv[val] = True

        for i, k in enumerate(format.split(':')):
            format_kv[k] = format_vals.split(':')[i]

        if gvcf_block_end > 0:
            # we are in a homozygous ref block.
            if not only_alt:
                # if the MIN_DP is above DP threshold, write it

                if not only_passing or int(format_kv['MIN_DP']) >= min_dp:
                    write_ref_block(chrom, pos, gvcf_block_end, int(format_kv['DP']), int(format_kv['MIN_DP']), fasta, skip_ambiguous, int(format_kv['MIN_DP']) < min_dp)
        else:
            # We aren't in a homozygous ref block...

            # find the depth
            dp = -1
            if 'DP' in format_kv:
                dp = int(format_kv['DP'])
            elif 'DP' in info_kv:
                dp = int(info_kv['DP'])

            # add LowDepth filter (or continue)
            if dp < min_dp:
                if only_passing:
                    continue
                elif filter in ['PASS', '.']:
                    filter = 'LowDepth'
                else:
                    filter = '%s;LowDepth' % filter

            # remove bad format vals first
            if format != '.':
                new_vals = []
                new_format = []
                for i, val in enumerate(format.split(':')):
                    if val not in format_keys_to_remove:
                        new_format.append(val)
                        new_vals.append(format_kv[val])

                format = ':'.join(new_format)
                format_vals = ':'.join(new_vals)

            # alter A/R format vals (and recalc GT if necessary)
            new_format_vals = []
            for i, k in enumerate(format.split(':')):
                if k == 'GT' and alts != '<NON_REF>' and recalc_gt:
                    # only recalculate GT for positions with a variant

                    ad_vals = format_kv['AD'].split(',')[:-1] # removes the NON_REF count

                    counts = []
                    gt = []
                    total = 0

                    for i, count in enumerate(ad_vals):
                        counts.append((int(count), i))
                        total += int(count)

                    counts = sorted(counts, reverse=True)

                    alt1_num = counts[0][1]
                    alt1_count = counts[0][0]

                    if alt1_count >= min_gt_count and alt1_count / total >= min_gt_af:
                        gt.append(alt1_num)

                    if len(counts) > 1:
                        alt2_num = counts[1][1]
                        alt2_count = counts[1][0]

                        if alt2_count >= min_gt_count and alt2_count / total >= min_gt_af:
                            gt.append(alt2_num)

                            if len(counts) > 2 and counts[2][0] == alt2_count:
                                # allele 1 and allele 2 have the same count (and both are valid) -- skip this (ref is 0)
                                gt = []
                    

                    if len(gt) == 0:
                        new_format_vals.append('./.') # this is strange and shouldn't happen
                    else:
                        gt = sorted(gt)
                        if len(gt) < 2:
                            new_format_vals.append('%s/%s' % (gt[0], gt[0]))
                        else:
                            new_format_vals.append('%s/%s' % (gt[0], gt[1]))
    
                elif k not in format_keys_to_adj:
                    new_format_vals.append(format_kv[k])

                else:
                    # remove the final value for a FORMAT field that includes the "<NON_REF>" alt
                    spl = format_kv[k].split(',')
                    new_format_vals.append(','.join(spl[:-1]))

            format_vals = ':'.join(new_format_vals)

            # alter A/R INFO  TODO: FIXME
            if info != '.':
                new_info = []
                for val in info.split(';'):
                    if '=' in val:
                        kv = val.split('=')
                        if kv[0] not in info_keys_to_adj:
                            new_info.append(val)
                        else:
                            spl = kv[1].split(',')
                            new_info.append('%s=%s' % (kv[0], ','.join(spl[:-1])))
                    else:
                        new_info.append(val)        
                info = ';'.join(new_info)

            if alts == '*,<NON_REF>':
                # this is a deletion... this is covered in an earlier VCF record, but not necessarily as easily queryable
                # however, this isn't part of the VCF spec, so we'll add an INFO flag and leave it

                if not info or info == '.':
                    info="DELETION"
                else:
                    info = '%s;DELETION' % info

                if not only_alt:
                    write_ref_single(chrom, pos, id, refbase, qual, filter, info, format, format_vals)

            elif alts == '<NON_REF>':
                # we have a single-base REF
                if not only_alt:
                    write_ref_single(chrom, pos, id, refbase, qual, filter, info, format, format_vals)

            else:
                if len(refbase) > 1:
                    # this is the start of a deletion
                    if not info or info == '.':
                        info="DELETION"
                    else:
                        info = '%s;DELETION' % info

                write_alt_single(chrom, pos, id, refbase, alts, qual, filter, info, format, format_vals)

    fasta.close()
    if fname:
        f.close()


    pass


def csv_split(s):
    in_quote = False
    vals = []

    buf = ""

    for c in s:
        if c == '"':
            if in_quote:
                in_quote = False
            else:
                in_quote = True
            buf += c
        elif c == ',':
            if in_quote:
                buf += c
            else:
                vals.append(buf)
                buf = ""
        else:
            buf += c

    if buf:
        vals.append(buf)

    return vals


def write_alt_single(chrom, pos, id, refbase, alts, qual, filter, info, format, format_vals):

    # Note: GT should be replaced before this call

    new_alts = []
    for alt in alts.split(','):
        if alt != '<NON_REF>':
            new_alts.append(alt)

    outcols = [chrom, pos, id, refbase, ','.join(new_alts), qual, filter, info, format, format_vals]
    sys.stdout.write('%s\n' % '\t'.join([str(x) for x in outcols]))


def write_ref_single(chrom, pos, id, refbase, qual, filter, info, format, format_vals):
    # replace whatever GT call there is with "0/0" for a homozygous ref call
    format_vals_spl = format_vals.split(":")
    for i, k in enumerate(format.split(":")):
        if k == "GT":
            format_vals_spl[i] = "0/0"

    outcols = [chrom, pos, id, refbase, '.', qual, filter, info, format, ":".join(format_vals_spl)]
    sys.stdout.write('%s\n' % '\t'.join([str(x) for x in outcols]))


def write_ref_block(chrom, pos, end, dp, min_dp, fasta, skip_ambiguous, is_lowdepth):
    for i in range(pos, end+1):
        refbase = fasta.fetch(chrom, i-1, i)
        if skip_ambiguous and refbase in ['N', 'n', '.']:
            continue

        filter = '.'
        if is_lowdepth:
            filter = 'LowDepth'

        # hard-code GT as "0/0", and add MIN_DP to signify that the DP here is at least MIN_DP
        if pos == i and dp != '':
            outcols = [chrom, i, '.', refbase, '.', '.', filter, '.', 'GT:DP:MIN_DP', '0/0:%s:%s' % (dp, min_dp)]
        else:
            outcols = [chrom, i, '.', refbase, '.', '.', filter, '.', 'GT:MIN_DP', '0/0:%s' % min_dp]
        sys.stdout.write('%s\n' % '\t'.join([str(x) for x in outcols]))


if __name__ == '__main__':
    min_dp = 10
    min_gt_count = 1
    min_gt_af = 0.1
    only_passing = False
    only_alt = False
    skip_ambiguous = False
    recalc_gt = False

    vcf_fname = None
    fasta_fname = None
    
    last = None
    for arg in sys.argv[1:]:
        if last == '--min-dp':
            min_dp = int(arg)
            last = None
        elif last == '--min-gt-count':
            min_gt_count = int(arg)
            last = None
        elif last == '--min-gt-af':
            min_gt_af = int(arg)
            last = None
        elif last:
            help("Unknown option: %s" % last)
        elif arg == '--only-passing':
            only_passing = True
        elif arg == '--only-alt':
            only_alt = True
        elif arg == '--recalc-gt':
            recalc_gt = True
        elif arg == '--skip-ambiguous':
            skip_ambiguous = True
        elif arg in ['--min-dp', '--min-gt-count', '--min-gt-af']:
            last = arg
        elif arg in ['--help', '-h', '-help']:
            help()
        elif not vcf_fname and (arg == '-' or os.path.exists(arg)):
            vcf_fname = arg
        elif not fasta_fname and os.path.exists(arg):
            fasta_fname = arg
    
    if not vcf_fname or not fasta_fname:
        help("Missing input file! (use - for stdin)")

    convert_gvcf(vcf_fname, fasta_fname, min_dp, recalc_gt, min_gt_count, min_gt_af, only_passing, only_alt, skip_ambiguous)
