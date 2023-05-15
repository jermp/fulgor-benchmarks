SEQKIT=/home/noor/bin/seqkit
map_contig=$1
ref_path=$3
ref_contigs=${ref_path}/ref_contigs
i=1
while read ref_file;do
    $SEQKIT fx2tab -n -i -H -l -q $ref_file | awk -F '\t' '{print $1}' | sed '1d' > $ref_contigs
    sed -i "s/^/Reference:${i}_Sequence:/g" $ref_contigs
    cat $ref_contigs >> $map_contig
    if [[ $((i%100)) -eq 0 ]]; then
        echo $i
    fi
    i=$((i+1))
    rm $ref_contigs
done < $2