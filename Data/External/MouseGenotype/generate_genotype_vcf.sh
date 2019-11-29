GENOTYPE_FOLDER=.
bcftools view -s 129P2/OlaHsd ${GENOTYPE_FOLDER}/SC_MOUSE_GENOMES.genotype.vcf.gz | bcftools view -o 129P2.vcf.gz -O z --exclude 'GT="0/0" | GT="./."' 
tabix -p vcf 129P2.vcf.gz
bcftools view -s C57BL/6NJ ${GENOTYPE_FOLDER}/SC_MOUSE_GENOMES.genotype.vcf.gz | bcftools view -o C57BL6.vcf.gz -O z --exclude 'GT="0/0" | GT="./."'
tabix -p vcf C57BL6.vcf.gz
bcftools merge -o 129P2_and_C57BL6_calls.vcf.gz -O z 129P2.vcf.gz C57BL6.vcf.gz
tabix -p vcf 129P2_and_C57BL6_calls.vcf.gz

