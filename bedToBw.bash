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


bedToBw 2020-04-23-WT_6h_dTAG_Flavo_Exp3_New-dedup-hg38.bed 1.0535 /media/genomes/hg38/hg38.chrom.sizes &
bedToBw 2020-04-23-WT_6h_dTAG_6h_PFA_Flavo_Exp3_New-dedup-hg38.bed 0.9517 /media/genomes/hg38/hg38.chrom.sizes &

bedToBw 2020-04-23-UL87F_No_dTAG_No_Flavo_Exp3_New-dedup-hg38.bed 1.1031 /media/genomes/hg38/hg38.chrom.sizes &
bedToBw 2020-04-23-UL87F_6h_dTAG_No_Flavo_Exp3_New-dedup-hg38.bed 0.8685 /media/genomes/hg38/hg38.chrom.sizes &

bedToBw 2020-04-23-UL87F_No_dTAG_Flavo_Exp3_New-dedup-hg38.bed 1.1302 /media/genomes/hg38/hg38.chrom.sizes &
bedToBw 2020-04-23-UL87F_6h_dTAG_Flavo_Exp3_New-dedup-hg38.bed 1.0419 /media/genomes/hg38/hg38.chrom.sizes &

bedToBw 2020-04-23-UL87F_6h_dTAG_6h_PFA_No_Flavo_Exp3_New-dedup-hg38.bed 0.8822 /media/genomes/hg38/hg38.chrom.sizes &
bedToBw 2020-04-23-UL87F_6h_dTAG_6h_PFA_Flavo_Exp3_New-dedup-hg38.bed  1.0770 /media/genomes/hg38/hg38.chrom.sizes &

bedToBw 2020-04-23-UL87F_No_dTAG_6h_PFA_No_Flavo_Exp3_New-dedup-hg38.bed 1.2171 /media/genomes/hg38/hg38.chrom.sizes &
bedToBw 2020-04-23-UL87F_No_dTAG_6h_PFA_Flavo_Exp3_New-dedup-hg38.bed 0.8176 /media/genomes/hg38/hg38.chrom.sizes 



