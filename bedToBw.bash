function bedToBw {

file=$1
scale=$2
chromSizeFile=$3

cat $file | awk '$1 != "chrEBV" {print $0}' > ${file%.bed}".tmp"
bedSort ${file%.bed}".tmp" $file
bedtools genomecov -scale $scale -i $file -g $chromSizeFile -bg -strand + > ${file%.bed}"-FW.bedGraph"
bedtools genomecov -scale $scale -i $file -g $chromSizeFile -bg -strand - > ${file%.bed}"-RV.bedGraph"
bedtools genomecov -scale $scale -i $file -g $chromSizeFile -bg -strand + -5 > ${file%.bed}"-FW-5.bedGraph"
bedtools genomecov -scale $scale -i $file -g $chromSizeFile -bg -strand + -3 > ${file%.bed}"-FW-3.bedGraph"
bedtools genomecov -scale $scale -i $file -g $chromSizeFile -bg -strand - -5 > ${file%.bed}"-RV-5.bedGraph"
bedtools genomecov -scale $scale -i $file -g $chromSizeFile -bg -strand - -3 > ${file%.bed}"-RV-3.bedGraph"
bedSort ${file%.bed}"-FW.bedGraph" ${file%.bed}"-FW.bedGraph"
bedSort ${file%.bed}"-RV.bedGraph" ${file%.bed}"-RV.bedGraph"
bedSort ${file%.bed}"-FW-5.bedGraph" ${file%.bed}"-FW-5.bedGraph"
bedSort ${file%.bed}"-FW-3.bedGraph" ${file%.bed}"-FW-3.bedGraph"
bedSort ${file%.bed}"-RV-5.bedGraph" ${file%.bed}"-RV-5.bedGraph"
bedSort ${file%.bed}"-RV-3.bedGraph" ${file%.bed}"-RV-3.bedGraph"
bedGraphToBigWig ${file%.bed}"-FW.bedGraph" $chromSizeFile ${file%.bed}"-FW.bw"
bedGraphToBigWig ${file%.bed}"-RV.bedGraph" $chromSizeFile ${file%.bed}"-RV.bw"
bedGraphToBigWig ${file%.bed}"-FW-5.bedGraph" $chromSizeFile ${file%.bed}"-FW-5.bw"
bedGraphToBigWig ${file%.bed}"-FW-3.bedGraph" $chromSizeFile ${file%.bed}"-FW-3.bw"
bedGraphToBigWig ${file%.bed}"-RV-5.bedGraph" $chromSizeFile ${file%.bed}"-RV-5.bw"
bedGraphToBigWig ${file%.bed}"-RV-3.bedGraph" $chromSizeFile ${file%.bed}"-RV-3.bw"
rm ${file%.bed}"-RV.bedGraph" ${file%.bed}"-FW.bedGraph" ${file%.bed}".tmp"
rm ${file%.bed}"-RV-5.bedGraph" ${file%.bed}"-FW-5.bedGraph"
rm ${file%.bed}"-RV-3.bedGraph" ${file%.bed}"-FW-3.bedGraph"

}

genome="/media/genomes/combinedgenome/combinedhg38_moth.chrom.sizes"

bedToBw TAF1-DMSO-Rep1-combined.bed 1.705910314 $genome
bedToBw TAF1-DMSO-Rep2-combined.bed 1.51747873 $genome
bedToBw TAF1-VHL-Rep1-combined.bed 1.08858964 $genome
bedToBw TAF1-VHL-Rep2-combined.bed 1.315846809 $genome
bedToBw TAF4-DMSO-Rep1-combined.bed 1.307716169 $genome
bedToBw TAF4-VHL-Rep1-combined.bed 1.711831066 $genome
bedToBw TBP-DMSO-Rep1-combined.bed 0.969276449 $genome
bedToBw TBP-DMSO-Rep2-combined.bed 0.795498868 $genome
bedToBw TBP-VHL-Rep1-combined.bed 0.759328409 $genome
bedToBw TBP-VHL-Rep2-combined.bed 0.826110089 $genome
bedToBw TFIIB-DMSO-Rep1-combined.bed 1.278398188 $genome
bedToBw TFIIB-DMSO-Rep2-combined.bed 1.093492191 $genome
bedToBw TFIIB-VHL-Rep1-combined.bed 0.652399366 $genome
bedToBw TFIIB-VHL-Rep2-combined.bed 0.43608327 $genome
bedToBw XPB-DMSO-Rep1-combined.bed 1.180759755 $genome
bedToBw XPB-DMSO-Rep2-combined.bed 1.115740976 $genome
bedToBw XPB-VHL-Rep1-combined.bed 1.07467902 $genome
bedToBw XPB-VHL-Rep2-combined.bed 0.982029653 $genome


