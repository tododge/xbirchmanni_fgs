NUMBER=20
genome=/scratch/groups/schumer/tris/falsegravid/gwas/xbir_WT/ref/xbir-COAC-16-VIII-22-M_v2023.1.fa

DIR=/scratch/groups/schumer/tris/falsegravid/gwas/xbir_WT/data/

ls $DIR/*read_1*.gz > R1list
ls $DIR/*read_2*.gz > R2list
paste R1list R2list > allfileslist

split allfileslist allfileslist. -l $NUMBER

ls allfileslist.* > listoflists

for line in `cat listoflists`
do
    touch $line.mapping.sh
    cat /home/groups/schumer/lab_member_folders/tris/sbatch_header_12hrs.txt | while read n
    do
          echo $n >> $line.mapping.sh
    done

  echo 'perl run_samtools_only_parental_v2.pl' $line $genome 'PE' >> $line.mapping.sh
done

#for file in `ls *mapping.sh`
#do
#    sbatch $file
#done
