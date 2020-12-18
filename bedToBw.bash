function bedToBw {

file=$1
scale=$2
chromSizeFile=$3

cat $file | awk '$1 != "chrEBV" {print $0}' > ${file%.bed}".tmp"
bedSort ${file%.bed}".tmp" $file
bedtools genomecov -scale $scale -i $file -g $chromSizeFile -bg -strand + > ${file%.bed}"-FW.bedGraph"
bedtools genomecov -scale $scale -i $file -g $chromSizeFile -bg -strand - > ${file%.bed}"-RV.bedGraph"
bedSort ${file%.bed}"-FW.bedGraph" ${file%.bed}"-FW.bedGraph"
bedSort ${file%.bed}"-RV.bedGraph" ${file%.bed}"-RV.bedGraph"
bedGraphToBigWig ${file%.bed}"-FW.bedGraph" $chromSizeFile ${file%.bed}"-FW.bw"
bedGraphToBigWig ${file%.bed}"-RV.bedGraph" $chromSizeFile ${file%.bed}"-RV.bw"
rm ${file%.bed}"-RV.bedGraph" ${file%.bed}"-FW.bedGraph" ${file%.bed}".tmp"

}

genome="/media/genomes/combinedgenome/combinedhg38_moth.chrom.sizes"

bedToBw Sample21-dedup.bed 0.876971972 $genome &
bedToBw Sample22-dedup.bed 1.29952819 $genome &
bedToBw Sample23-dedup.bed 1.825712356 $genome &
bedToBw Sample24-dedup.bed 2.771789905 $genome &
bedToBw Sample25-dedup.bed 1.273699439 $genome &
bedToBw Sample26-dedup.bed 1.18758729 $genome &
bedToBw Sample27-dedup.bed 0.601445391 $genome &
bedToBw Sample28-dedup.bed 0.799385596 $genome &
bedToBw Sample29-dedup.bed 1.157654349 $genome &
bedToBw Sample30-dedup.bed 1.479850181 $genome &
bedToBw Sample31-dedup.bed 0.726676314 $genome &
bedToBw Sample32-dedup.bed 0.557766921 $genome &
bedToBw Sample33-dedup.bed 0.872404903 $genome &
bedToBw Sample34-dedup.bed 1.731372368 $genome &
bedToBw Sample35-dedup.bed 1.375630735 $genome &
bedToBw Sample36-dedup.bed 2.587498112 $genome &
bedToBw Sample37-dedup.bed 0.679189383 $genome &
bedToBw Sample38-dedup.bed 0.706914671 $genome &
bedToBw Sample39-dedup.bed 0.863500812 $genome &
bedToBw Sample40-dedup.bed 0.953032069 $genome &
