MASON="mason_simulator"

# SAMPLED="pseudoaligner_paper/ENA2018-bacteria-661k-filenames/filenames_list_10k_norm.txt"
SAMPLED=$1
echo $SAMPLED
# NON_SAMPLED_PATH="/mnt/scratch5/wabi2023/chr19.fa"
NON_SAMPLED_PATH=$2
RATIO=$3
OUT_PATH=$4
OUT_EXTEND=$5
NUM_READS=$6
FILE_EXTEND=$OUT_EXTEND"_ratio="$RATIO
SEQKIT=/home/noor/bin/seqkit

NUM_SAMP=`wc -l $SAMPLED | cut -d' ' -f1`
NUM_NONSAMP=`python3 -c "print(int(${NUM_SAMP}* float(${RATIO})/(100-${RATIO})))"`

echo $NUM_NONSAMP

TF=`tempfile`
echo $SAM
echo "creating $NUM_NONSAMP samples of $NON_SAMPLED_PATH in $TF"
for i in `seq $NUM_NONSAMP`;
do 
  echo $NON_SAMPLED_PATH >> $TF;
done

echo "cat $SAMPLED $TF > TEMP_GENOMES.txt"
cat $SAMPLED $TF > TEMP_GENOMES.txt

rm $TF

for ref in `shuf -n 500 TEMP_GENOMES.txt`; do
  
  SEED=10
  while [ $SEED != 1 ]
  do
    genome_fa=genome_ratio=${RATIO}_nquery=${OUT_EXTEND}.fa
    $SEQKIT sort -l -r $ref | $SEQKIT head -n10 | $SEQKIT shuffle -s $SEED | $SEQKIT head -n 1 > ${genome_fa}
    read_prefix=`head -1 ${genome_fa} | awk '{print substr($1,2)}'`
    $MASON --read-name-prefix "$read_prefix:" -ir ${genome_fa} -n $NUM_READS -o tmp.fq
    if [ $? -eq 0 ]
    then
      cat tmp.fq >> $OUT_PATH/sim_reads_$FILE_EXTEND.fq
      echo $SEED $ref >> $OUT_PATH/seeds_$FILE_EXTEND
      SEED=1
      
    else
      rm tmp.fq
      echo $SEED
      SEED=$((SEED+10))
    fi
  rm ${genome_fa}
  rm ${genome_fa}.fai
  
  done
done

rm TEMP_GENOMES.txt
rm tmp.fq