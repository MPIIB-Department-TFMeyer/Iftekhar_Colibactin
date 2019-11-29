#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Fri May  5 10:23:32 2017

@author: hilmar
"""
import vcf
import re

# order taken from http://snpeff.sourceforge.net/VCFannotationformat_v1.0.pdf
annotype_sort_order = """chromosome_number_variation
exon_loss_variant
frameshift_variant
stop_gained
stop_lost
start_lost
splice_acceptor_variant
splice_donor_variant
rare_amino_acid_variant
missense_variant
inframe_insertion
disruptive_inframe_insertion
inframe_deletion
disruptive_inframe_deletion
5_prime_UTR_truncation+exon_loss_variant
3_prime_UTR_truncation+exon_loss
splice_branch_variant
splice_region_variant
splice_branch_variant
stop_retained_variant
initiator_codon_variant
synonymous_variant
initiator_codon_variant+non_canonical_start_codon
stop_retained_variant
coding_sequence_variant
5_prime_UTR_variant
3_prime_UTR_variant
5_prime_UTR_premature_start_codon_gain_variant
upstream_gene_variant
downstream_gene_variant
TF_binding_site_variant
regulatory_region_variant
miRNA
custom
sequence_feature
conserved_intron_variant
intron_variant
intragenic_variant
conserved_intergenic_variant
intergenic_region
coding_sequence_variant
non_coding_exon_variant
nc_transcript_variant
gene_variant
chromosome"""

annotype_order_hash = {}
tmp = annotype_sort_order.strip().split("\n")
for ii in tmp:
    annotype_order_hash[ii]=tmp.index(ii)
annotype_order_hash_max = max(annotype_order_hash.values())+1

impact_order = """HIGH
MODERATE
LOW
MODIFIER
"""

impact_order_hash = {}
tmp = impact_order.strip().split("\n")
for ii in tmp:
    impact_order_hash[ii]=tmp.index(ii)


def vcf_to_table(ifile, ofile):
    vcf_reader = vcf.Reader(ifile)
    
    ann_fields = vcf_reader.infos['ANN']
    af_list =  ann_fields.desc.replace("Functional annotations:","").strip(" '").split(" | ")
    af_ind_hash = {}
    for a in af_list:
        af_ind_hash[a] = af_list.index(a)

    rec_cnt = 0
    
    # global overview on relevant variants per gene
    per_gene_variants = {'LOF': {}, 'NMD': {}, 'ANN': {}}

    unknown_annotation_types = {}
    
    for record in vcf_reader:
        if rec_cnt == 0:
            #ofile.write("CHROM\tPOS\tREF\tTYPE\tFILTER\tALT\tSAMPLE\tDP\tAD\tGQ\tFT\tGene\tImpact\tAnnotationType\tAluSite\tdbSNP\tExonJunction")
            ofile.write("CHROM\tPOS\tREF\tTYPE\tFILTER\tALT\tSAMPLE\tDP\tAD\tGQ\tFT\tGene\tImpact\tAnnotationType\tdbSNP")
            ofile.write("\n")

#        if rec_cnt < 5:
#            print record

        chrom = record.CHROM
        pos = record.POS
        #rid = record.ID
        ref = record.REF
        alt = record.ALT
        #qual = record.QUAL

        # Replace missing FILTER value with empty string instead of None
        filter_stat = record.FILTER
        if not filter_stat is None:
            filter_str = "PASS" if len(filter_stat)==0 else ";".join(filter_stat)
        else:
            filter_str = "."
        
        #print record.INFO
        # Replace missing QUAL value with empty string instead of None
        #q_str = str(qual) if not qual is None else ""
        
        # the default columns of the variant
        out_rec = ("%s\t%d\t%s\t%s\t%s" % (chrom, pos, ref, record.var_type, filter_str))
        #print out_rec
        all_alleles = [ref]
        all_alleles.extend([str(e) for e in alt])
        
        # map alleles to samples
        variant_alleles_per_sample = {}
        for sample in record.samples:
            #print sample
            sn = sample.sample
            
            # exclude duplicate reference samples
            if sn[1]==":":
                continue
            
            # Strelka does not call the genotype for somativ variants
            if sample['DP']==".":
                gt_raw="0/0"
            else:
                gt_raw="1/1"
            #gt_raw = sample['GT']
            if gt_raw is None:
                var_allele_ind = None
            else:
                var_allele_ind = [int(e) for e in re.split("/|\|", gt_raw) if e != "0"]
            
            if not var_allele_ind is None and len(var_allele_ind)>0:
                if not sn in variant_alleles_per_sample:
                    variant_alleles_per_sample[sn] = {}

                for vi in var_allele_ind:
                    try:
                        depth = sample['DP'] 
                    except:
                        depth = sample['DPI']
                    #variant_alleles_per_sample[sn][all_alleles[vi]] = (depth, sample['AD'][vi], sample['GQ'], sample['FT']) 
                    variant_alleles_per_sample[sn][all_alleles[vi]] = (depth, -1, -1, -1) 
        
           
        # treat the SnpEff info parameters
        # Loss of function info
        """ From MacArthur et al Science 2012 p. 823-828
            We adopted a definition for LoF variants expected to correlate with complete loss of function of the affected transcripts: 
            stop codon–introducing (nonsense) or 
            splice site–disrupting single-nucleotide variants (SNVs), 
            insertion/deletion (indel) variants predicted to disrupt a transcript’s reading frame,
            or larger deletions removing either the first exon or more than 50% of the protein-coding sequence of the affected transcript. 
            """
        
        if 'LOF' in record.INFO:
            lof_record = record.INFO['LOF']
            for rec_str in lof_record:
                gg = rec_str.strip("()").split("|", -1)
                gene_symbol = gg[0]
                lof_result = (gg[1], gg[2], gg[3]) # gene id, number of transcripts affected, % of tx affected
                if not gene_symbol in per_gene_variants['LOF']:
                    per_gene_variants['LOF'][gene_symbol]= [lof_result]
                else:
                    per_gene_variants['LOF'][gene_symbol].append(lof_result)
        
        # nonsense-mediated decay
        # for definition see Nat Rev Mol Cell Biol. 2004 Feb;5(2):89-99.
        if 'NMD' in record.INFO:
            nmd_record = record.INFO['NMD']
            for rec_str in nmd_record:
                gg = rec_str.strip("()").split("|", -1)
                gene_symbol = gg[0]
                nmd_result = (gg[1], gg[2], gg[3]) # gene id, number of transcripts affected, % of tx affected
                if not gene_symbol in per_gene_variants['NMD']:
                    per_gene_variants['NMD'][gene_symbol]= [nmd_result]
                else:
                    per_gene_variants['NMD'][gene_symbol].append(nmd_result)

        # complete annotation
        if 'ANN' in record.INFO:
            ann_record = record.INFO['ANN']
            curr_anno = {}
            curr_anno2 = {}
            for rec_str in ann_record:
                ar = rec_str.split("|", -1)
                allele = ar[af_ind_hash['Allele']]
                atype = ar[af_ind_hash['Annotation']]
                impact = ar[af_ind_hash['Annotation_Impact']]
                gsymbol = ar[af_ind_hash['Gene_Name']]
                
                if not atype in annotype_order_hash:
                    if not atype in unknown_annotation_types:
                         unknown_annotation_types[atype] = 1
                    else:
                         unknown_annotation_types[atype] += 1
                
                if not allele in curr_anno:
                    curr_anno[allele] = {}
                
                if not gsymbol in curr_anno[allele]:
                    curr_anno[allele][gsymbol]={}
                
                if not impact in curr_anno[allele][gsymbol]:
                    curr_anno[allele][gsymbol][impact] = {}
                
                if not atype in curr_anno[allele][gsymbol][impact]:
                    curr_anno[allele][gsymbol][impact][atype] = []
                
                curr_anno[allele][gsymbol][impact][atype].append(ar)
                
                #############################################
                impact_id = (impact, atype)
                if not allele in curr_anno2:
                    curr_anno2[allele] = {}
                
                if not gsymbol in curr_anno2[allele]:
                    curr_anno2[allele][gsymbol]={}
                
                if not impact_id in curr_anno2[allele][gsymbol]:
                    curr_anno2[allele][gsymbol][impact_id] = []

                curr_anno2[allele][gsymbol][impact_id].append(ar)

        polymorphism_site="."
        if not record.ID is None and record.ID[0:2]=="rs":
            polymorphism_site = record.ID

                
#        print "-----------"
#        print curr_anno2
#        print "-----------"
#        print variant_alleles_per_sample
        # identify highest impact for each allele and gene


        for c_allele in curr_anno2.keys():
            for c_gene in curr_anno2[c_allele].keys():
                sorted_impacts = sorted(curr_anno2[c_allele][c_gene].keys(), key=lambda x: (impact_order_hash[x[0]], annotype_order_hash[x[1]] if x[1] in annotype_order_hash else annotype_order_hash_max))
                #print c_allele, c_gene, sorted_impacts

                highest_impact = sorted_impacts[0]
                gene_symbol = c_gene
        
                if not gene_symbol in per_gene_variants['ANN']:
                    per_gene_variants['ANN'][gene_symbol]= {}
                
                if not highest_impact in per_gene_variants['ANN'][gene_symbol]:
                    per_gene_variants['ANN'][gene_symbol][highest_impact] = []
                
                per_gene_variants['ANN'][gene_symbol][highest_impact].append((chrom, pos, alt))
                
                for s in variant_alleles_per_sample:
                    if c_allele in variant_alleles_per_sample[s]:
                        sample_gt = variant_alleles_per_sample[s][c_allele]
                        ofile.write(out_rec + "\t%s\t%s" % (c_allele, s) + ("\t%s\t%s\t%s\t%s" %  (sample_gt) ) + ("\t%s\t%s\t%s\t%s" % (c_gene, highest_impact[0], highest_impact[1], polymorphism_site  )) + "\n")
                       
        rec_cnt += 1
#        if rec_cnt > 5:
#            break

        #print(unknown_annotation_types)

if __name__ == '__main__':

    import argparse, bz2file, gzip, fileinput, sys
    import time
    parser = argparse.ArgumentParser(description='Convert SnpEff-annotated VCF to tall-skinny table format')
    parser.add_argument('input_vcf', type=str, help='Input vcf')

#    parser.add_argument('--exclude-uncultured', action='store_true', default=False,
#                       help='Exclude uncultured and unclassified bacteria as far as possible.')
    parser.add_argument('-o','--output-file', action='store', type=str, default=None,
                       help='File where output should go to. Defaults to STDOUT')


    args = parser.parse_args()

    if args.input_vcf.strip()=="-":
        ifile=sys.stdin
    elif args.input_vcf.endswith(".bz2"):
        try:        
            ifile=bz2file.BZ2File(args.input_vcf,"r", buffering=0)  
        except Exception, e:
            raise e
    elif args.input_vcf.endswith(".gz") or args.input_vcf.endswith(".bgz"):
        try:
            ifile=gzip.GzipFile(args.input_vcf,"r")
        except Exception, e:
            raise e
    else:
        try:
            ifile=open(args.input_vcf,"r")
        except Exception, e:
            raise e

    if args.output_file is None:    
        ofile = sys.stdout
    else:
        ofile = open(args.output_file,"w")
       
#

#    start = time.clock()
    vcf_to_table(ifile, ofile)
 #   end = time.clock()
 #   diff = (end-start)

 

