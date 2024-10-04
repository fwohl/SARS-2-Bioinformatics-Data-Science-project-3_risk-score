#Download data
wget --no-check-certificate https://osf.io/download/vd7gj/ -O vd7gj.fasta
#download reference
wget --no-check-certificate https://osf.io/download/kt4dq/ -O nCoV-2019.reference.fasta
# we also download an index for later usage
wget --no-check-certificate https://osf.io/download/5eunp/ -O nCoV-2019.reference.fasta.fai
#activate environment
conda activate /storage/mi/transat/workshop/workshop/envs/workshop/



#map reads and generate sam
minimap2 -x sr -t 4 -a -o vd7gj.sam nCoV-2019.reference.fasta vd7gj.fasta
#convert to bam
samtools view -bS vd7gj.sam > vd7gj.bam
#sort
samtools sort vd7gj.bam > vd7gj.sorted.bam
#index
samtools index vd7gj.sorted.bam

#make merged vcf file (we dont work with it but do it anyways)
freebayes -f nCoV-2019.reference.fasta vd7gj.sorted.bam > merge_vcf.vcf



#annotate lineage with singularity

#
conda install singularity
singularity pull --name nanozoo-pangolin-v4-4.3--1.30.img docker://nanozoo/pangolin-v4:4.3--1.30
singularity run nanozoo-pangolin-v4-4.3--1.30.img
pangolin -t 4 vd7gj.fasta
#exit singularity
exit
#make tsv
tr ',' '\t' <lineage_report.csv >lineage_report.tsv

#calculate the number of ambigous nucleotides per sample
awk -vRS=">" -vORS="\n" -vFS="\n" -vOFS="\t" 'BEGIN {NR>1; FS = "\n" }; {print $1"\t"gsub(/N/,"")}' vd7gj.fasta >number_of_ns.txt

#make file of seqnames
awk -vRS=">" -vORS="\n" -vFS="\n" -vOFS="\t" 'BEGIN {NR>1; FS = "\n" }; /seq*/ {print $1}' vd7gj.fasta >seqnames.txt

#make array of seqnames
while IFS= read -r line; do   arr+=("$line"); done < seqnames.txt

#covsonar

git clone https://github.com/rki-mf1/covsonar.git
mamba env create -p envs/sonar -f covsonar/sonar.env.yml
conda activate ~/workshop/envs/sonar
#make database
~/covsonar/sonar.py add -f vd7gj.fasta --db mydb --cpus 4

#add metadata
~/covsonar/sonar.py update --tsv lineage_report.tsv --fields accession=taxon lineage=lineage --db mydb
#make one big file for all variants
~/covsonar/sonar.py match --db mydb > all_variants.txt

#extract lineage infomration
for value in "${arr[@]}"; do ~/covsonar/sonar.py match --db mydb --acc "$value" | awk -v FS=" " -v value="$value" '{print $0}'| awk -v FS="," -v value="$value" 'NR>1 {print value "\t" $18}' ;done >lineages.txt
for value in "${arr[@]}"; do ~/covsonar/sonar.py match --db mydb --acc "$value" | awk -v FS=" " -v value="$value" '{print $0}'| awk -v FS="," -v value="$value" 'NR>1 {print $18}'|tr "\n" " "; done | xargs -n1 | sort |uniq > unique_lineages.txt


#make vcf files for each sequence
for value in "${arr[@]}"; do ~/covsonar/sonar.py var2vcf --db mydb --acc "$value" -o ${value}merge.vcf  ; done
#do the stats
for value in "${arr[@]}"; do bcftools stats  ${value}merge.vcf.gz  > ${value}.txt ; done

#count deletions
for value in "${arr[@]}"; do bcftools stats  ${value}merge.vcf.gz | awk -v FS=" " -v value="$value" '/.*number of no-ALTs:.*/ {print value"\t" $6}'>>no_deletions.txt; done
#count number of SNPs
#
for value in "${arr[@]}"; do bcftools stats  ${value}merge.vcf.gz | awk -v FS=" " -v value="$value" '/.*number of SNPs:.*/ {print value"\t" $6}'>>no_snps.txt; done
#get all mutations on dns profile
#

for value in "${arr[@]}"; do ~/covsonar/sonar.py match --db mydb --acc "$value" | awk -v FS=" " -v value="$value" '{print $0}'| awk -v FS="," -v value="$value" 'NR>1 {print value "\t" $20}'>> allmutations.txt; done

#get all aa mutations
for value in "${arr[@]}"; do ~/covsonar/sonar.py match --db mydb --acc "$value" | awk -v FS=" " -v value="$value" '{print $0}'| awk -v FS="," -v value="$value" 'NR>1 {print value "\t" $21}'>> all_aa_mutations.txt; done
#get all deletions

for value in "${arr[@]}"; do ~/covsonar/sonar.py match --db mydb --acc "$value" | awk -v FS=" " -v value="$value" '{print $0}'| awk -v FS="," -v value="$value" 'NR>1 {print value "\t" $22}'>> all_deletions.txt; done



#list of unique mutations
for value in "${arr[@]}"; do ~/covsonar/sonar.py match --db mydb --acc "$value" | awk -v FS=" " -v value="$value" '{print $0}'| awk -v FS="," -v value="$value" 'NR>1 {print $20}'|tr "\n" " "; done | xargs -n1 | sort |uniq  >>list_of_unique_mutations.txt
#same for spike mutations
for value in "${arr[@]}"; do ~/covsonar/sonar.py match --db mydb --acc "$value" | awk -v FS=" " -v value="$value" '{print $0}'| awk -v FS="," -v value="$value" 'NR>1 {print $20}'|tr "\n" " "; done | xargs -n1 | sort |uniq | awk ' /*2156[3-9]|216[0-9]|2[2-4][0-9]{3}|25[0-2][0-9]{2}|253[0-7][0-9]|2538[0-4]*/ {print $0}'>>list_of_unique_mutations_sspike.txt


#count
for value in "${arr[@]}"; do ~/covsonar/sonar.py match --db mydb --acc "$value" | awk -v FS=" " -v value="$value" '{print $0}'| awk -v FS="," -v value="$value" 'NR>1 {print $20}'|tr "\n" " "; done | xargs -n1 | sort |uniq | wc -w 

for value in "${arr[@]}"; do ~/covsonar/sonar.py match --db mydb --acc "$value" | awk -v FS=" " -v value="$value" '{print $0}'| awk -v FS="," -v value="$value" 'NR>1 {print $20}'|tr "\n" " "; done | xargs -n1 | sort |uniq | awk ' /*2156[3-9]|216[0-9]|2[2-4][0-9]{3}|25[0-2][0-9]{2}|253[0-7][0-9]|2538[0-4]*/ {print $0}'|wc -w

#"" deletions
for value in "${arr[@]}"; do ~/covsonar/sonar.py match --db mydb --acc "$value" | awk -v FS=" " -v value="$value" '{print $0}'| awk -v FS="," -v value="$value" 'NR>1 {print $20}'|tr "\n" " "; done | xargs -n1 | sort |uniq  |awk -v RS="\n" '/del:*/ {print $0}'>list_of_unique_deletions.txt


for value in "${arr[@]}"; do ~/covsonar/sonar.py match --db mydb --acc "$value" | awk -v FS=" " -v value="$value" '{print $0}'| awk -v FS="," -v value="$value" 'NR>1 {print $20}'|tr "\n" " "; done | xargs -n1 | sort |uniq  |awk -v RS="\n" '/del:*/ {print $0}'| wc -w



#same for spike
for value in "${arr[@]}"; do ~/covsonar/sonar.py match --db mydb --acc "$value" | awk -v FS=" " -v value="$value" '{print $0}'| awk -v FS="," -v value="$value" 'NR>1 {print $22}'|tr "\n" " "; done | xargs -n1 | sort |uniq | awk -v RS="\n" '/del:*/ {print $0}' | awk ' /*2156[3-9]|216[0-9]|2[2-4][0-9]{3}|25[0-2][0-9]{2}|253[0-7][0-9]|2538[0-4]*/ {print $0}'>>list_of_unique_deletions_spike.txt



#count
#
for value in "${arr[@]}"; do ~/covsonar/sonar.py match --db mydb --acc "$value" | awk -v FS=" " -v value="$value" '{print $0}'| awk -v FS="," -v value="$value" 'NR>1 {print $22}'|tr "\n" " "; done | xargs -n1 | sort |uniq | wc -w

for value in "${arr[@]}"; do ~/covsonar/sonar.py match --db mydb --acc "$value" | awk -v FS=" " -v value="$value" '{print $0}'| awk -v FS="," -v value="$value" 'NR>1 {print $22}'|tr "\n" " "; done | xargs -n1 | sort |uniq | awk -v RS="\n" '/del:*/ {print $0}' | awk ' /*2156[3-9]|216[0-9]|2[2-4][0-9]{3}|25[0-2][0-9]{2}|253[0-7][0-9]|2538[0-4]*/ {print $0}'| wc -w

	
#""aa mutations
#
for value in "${arr[@]}"; do ~/covsonar/sonar.py match --db mydb --acc "$value" | awk -v FS=" " -v value="$value" '{print $0}'| awk -v FS="," -v value="$value" 'NR>1 {print $21}'|tr "\n" " "; done| xargs -n1 | sort |uniq >>listofunique_aa_mutations.txt


#count
for value in "${arr[@]}"; do ~/covsonar/sonar.py match --db mydb --acc "$value" | awk -v FS=" " -v value="$value" '{print $0}'| awk -v FS="," -v value="$value" 'NR>1 {print $21}'|tr "\n" " "; done| xargs -n1 | sort |uniq | wc -w
#prepare matrix
for value in "${arr[@]}"; do ~/covsonar/sonar.py match --db mydb --acc "$value" | awk -v FS=" " -v value="$value" '{print $0}'| awk -v FS="," -v value="$value" 'NR>1 {print value"\t" $20}'|tr "\n" " "| xargs -n1 | sort > ${value}_mutations.txt ; done

unset matrix
declare -A matrix
num_rows=46
num_columns=443


rm outtest.txt

for ((i=1;i<=num_rows;i++)) do         value=$(sed -n "${i}{p;q;}" seqnames.txt );     for ((j=1;j<=num_columns;j++)) do         name=$(sed -n "${j}{p;q;}" list_of_unique_mutations.txt);         if grep -Fxq "$name"  "$value"_mutations.txt;                then echo "1" >> outtest.txt;         else                echo "0" >> outtest.txt;         fi;     done;  done



#echo "${matrix[@]}" > matrix.txt

#same for spike
#
unset matrix
declare -A matrix
num_rows=46
num_columns=443



rm outtest2.txt
for ((i=1;i<=num_rows;i++)) do         value=$(sed -n "${i}{p;q;}" seqnames.txt );     for ((j=1;j<=num_columns;j++)) do         name=$(sed -n "${j}{p;q;}" list_of_unique_mutations_sspike.txt);         if grep -Fxq "$name"  "$value"_mutations.txt;                then echo "1" >> outtest2.txt;         else                echo "0" >> outtest2.txt;         fi;     done;  done



#statistics mutations
#seq level
for ((i=1;i<=num_rows;i++)) do         value=$(sed -n "${i}{p;q;}" seqnames.txt );count=0;     for ((j=1;j<=73;j++)) do         name=$(sed -n "${j}{p;q;}" allmutationsspike.txt );         if grep -Fxq "$name"  "$value"_mutations.txt;                then count=$(( count+1 ));          fi;done; echo "$value" "$count"  >>stats_mut_spike.txt; done

for ((i=1;i<=num_rows;i++)) do         value=$(sed -n "${i}{p;q;}" seqnames.txt );count=0;     for ((j=1;j<=num_columns;j++)) do         name=$(sed -n "${j}{p;q;}" listofuniquemutations.txt );         if grep -Fxq "$name"  "$value"_mutations.txt;                then count=$(( count+1 ));          fi;done; echo "$value" "$count"  >>stats_mut.txt; done

#for deletions
for ((i=1;i<=num_rows;i++)) do         value=$(sed -n "${i}{p;q;}" seqnames.txt );count=0;     for ((j=1;j<=num_columns;j++)) do         name=$(sed -n "${j}{p;q;}" listofuniquedeletions.txt );         if grep -Fxq "$name"  "$value"_mutations.txt;                then count=$(( count+1 ));          fi;done; echo "$value" "$count" >>stats_del.txt; done

#spike
for ((i=1;i<=num_rows;i++)) do         value=$(sed -n "${i}{p;q;}" seqnames.txt );count=0;     for ((j=1;j<=num_columns;j++)) do         name=$(sed -n "${j}{p;q;}" alldeletionsspike.txt );         if grep -Fxq "$name"  "$value"_mutations.txt;                then count=$(( count+1 ));          fi;done; echo "$value" "$count" >>stats_del_spike.txt; done


#genome level
for ((i=1;i<=num_columns;i++)) do        name=$(sed -n "${i}{p;q;}" listofuniquemutations.txt )  ;count=0;     for ((j=1;j<=num_rows;j++)) do         value=$(sed -n "${j}{p;q;}" seqnames.txt );         if grep -Fxq "$name"  "$value"_mutations.txt;                then count=$(( count+1 ));          fi;done; echo "$name" "$count"  >>stats_seq.txt; done

for ((i=1;i<=73;i++)) do        name=$(sed -n "${i}{p;q;}" allmutationsspike.txt )  ;count=0;     for ((j=1;j<=num_rows;j++)) do         value=$(sed -n "${j}{p;q;}" seqnames.txt );         if grep -Fxq "$name"  "$value"_mutations.txt;                then count=$(( count+1 ));          fi;done; echo "$name" "$count"  >>stats_seq_spike.txt; done


for ((i=1;i<=num_columns;i++)) do        name=$(sed -n "${i}{p;q;}" listofuniquedeletions.txt )  ;count=0;     for ((j=1;j<=num_rows;j++)) do         value=$(sed -n "${j}{p;q;}" seqnames.txt );         if grep -Fxq "$name"  "$value"_mutations.txt;                then count=$(( count+1 ));          fi;done; echo "$name" "$count"  >>stats_seq_del.txt; done

for ((i=1;i<=73;i++)) do        name=$(sed -n "${i}{p;q;}" alldeletionsspike.txt )  ;count=0;     for ((j=1;j<=num_rows;j++)) do         value=$(sed -n "${j}{p;q;}" seqnames.txt );         if grep -Fxq "$name"  "$value"_mutations.txt;                then count=$(( count+1 ));          fi;done; echo "$name" "$count"  >>stats_seq_del_spike.txt; done

