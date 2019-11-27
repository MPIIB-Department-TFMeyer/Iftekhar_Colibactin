import sys, os

ifile = sys.stdin

per_sample_files = {}
for line in ifile:
    fasta_file = line.strip("\n\r")
    run = os.path.split(os.path.dirname(fasta_file))[-1]
    ff = os.path.basename(fasta_file).replace("_f.fastq.gz","").replace(".fastq.gz","")
    ff_fields = ff.split("_")
    sample = "_".join(ff_fields[0:2])
    run = ff_fields[2]
    index = ff_fields[3]
    sampleID = ff_fields[4]
    lane = ff_fields[5]
    pair_end = ff_fields[6]
    part = ff_fields[7]

    fkey = (sample)
    complete_key = (sample,lane,part)

    if not fkey in per_sample_files:
        per_sample_files[fkey] = {"R1": {}, "R2": {}}

    per_sample_files[fkey][pair_end][complete_key] = fasta_file

#print per_sample_files

for s in sorted(per_sample_files.keys()):
    reps = sorted(per_sample_files[s]["R1"].keys())
    r1 = []
    r2 = []
    for r in reps:
        r1.append(per_sample_files[s]["R1"][r])
        #r2.append(per_sample_files[s]["R2"][r])
    print "%s\t%s\t%s" % (s, ",".join(r1), ",".join(r2) )
    #ofile.write("%s\t%s\t%s\n" % ("_".join(s), ",".join(r1), ",".join(r2) ) )

