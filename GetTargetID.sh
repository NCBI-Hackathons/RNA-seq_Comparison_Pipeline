# Arguments 
# -- name of the abundance file
# Transcript ID 

export LC_ALL="C"
N=100
Gene2Ref=../data/gene2refseq.sdx
#Loop through the different directories and extract transcript ID and Counts
for d in */ ; do
head -$N $d/$1 | cut -f1 |sed -e 's/gi\(.*\)ref|/\l/'|cut -d'|' -f1>$d/TID.txt
head -$N $d/$1 | cut -f4>$d/Count.txt
#Now get gene ID from gene2refseq correspond to transcripts ID 
echo "geneId">$d/GID.txt;
while read line; do grep  "$line"  $Gene2Ref|cut -f2>>$d/GID.txt;  done <$d/TID.txt
#Now concatenate the values that we gathered in to a single file 
pr -mts  $d/TID.txt $d/Count.txt $d/GID.txt > $d/Comp.txt
#For each unique gene ID extract transcripts corresponds to the gene and find total Counts
echo 0> $d/Score.txt
for t in `uniq  $d/GID.txt` ; do
export Sum=$(grep $t $d/Comp.txt|cut -f2|paste -sd+ | bc)
NumElem=$(uniq  $d/GID.txt|wc -l)
Mean=$(echo $Sum/$NumElem|bc -l)
Sum2=$(echo 0|bc -l)
for e in `grep $t $d/Comp.txt|cut -f2` ; do
Sum2=$(echo $Sum2+$e*$e-2*$e*$Mean+$Mean*$Mean |bc -l)
done #grep $t $d/Comp.txt|cut -f2
Vari=$(echo $Sum2/$NumElem|bc -l)
echo $t $Mean $Vari >> $d/Score.txt
done #uniq  $d/GID.txt

done # loop over files */
echo ""> CountForSample.txt
for d in */ ; do
cp $d/Comp.txt ./Comp-$d.txt
done
