#Given Kallisto's abundance calculations per SRA run per transcript, grouped into sets originating from the same tissue.
#These grouped calculations are the output of fetchAndRun.sh, symlinked locally named thyroid.all and heart.all for two tissues

#We compute the standard deviation between runs originating from the same tissue for each transcript.
#Heuristically, we descard the lowest quintile ranked on (std_dev(expression)/expression) to discard the 
# transcripts whose reported abundance is most variable between runs of similar origin.  Subsequent analysis 
# can be done both on all transcripts and on the union of cleanest 80% of both.

#Then we compute fold change normalized by (sqrt(std_dev_1*std_dev_1 + std_dev_2*std_dev_2)).  

perl -ane '{$sum=0; $sum2=0; for ($i=1;$i<=$#F;$i++){$sum += $F[$i]; $sum2 += $F[$i]*$F[$i]; } $stdev = sqrt( ($sum2-$sum* $sum/($#F))/($#F) ); print $F[0],"\t",$sum/($#F),"\t",$stdev,"\n"}' thyroid.all > thyroid.sum
perl -ane '{$sum=0; $sum2=0; for ($i=1;$i<=$#F;$i++){$sum += $F[$i]; $sum2 += $F[$i]*$F[$i]; } $stdev = sqrt( ($sum2-$sum* $sum/($#F))/($#F) ); print $F[0],"\t",$sum/($#F),"\t",$stdev,"\n"}' heart.all > heart.sum
sort heart.sum > heart.sum.s
sort thyroid.sum > thyroid.sum.s
 join heart,sum.s thyroid.sum.s > heart_thyroid.sum 

 perl -ane '{print $F[0],"\t",$F[2]/(0.0001+$F[1]),"\t",$F[1],"\n"}' heart_thyroid.sum > heart.spread
 sort -k 2,2nr heart.spread | sed -n '1,32623p' | grep NM_ | sort -k 3,3n
 sort -k 2,2nr heart.spread | sed -n '32623,$p' | sort > heart.clean80
 perl -ane '{print $F[0],"\t",$F[4]/(0.0001+$F[3]),"\t",$F[3],"\n"}' heart_thyroid.sum > thymus.spread
 sort -k 2,2nr thymus.spread | sed -n '32623,$p' | sort > thymus.clean80
 join -0 1.1 thymus.clean80 heart.clean80 
 join -o 1.1 thymus.clean80 heart.clean80 
 join -o 1.1 thymus.clean80 heart.clean80 > both.clean80
 ls
 head heart_thyroid.ratio 
 perl -ane '{print $F[0],"\t",log((0.000001+$F[1])/(0.000001+$F[3])),"\t", sqrt($F[2]*$F[2]+$F[4]*$F[4]),"\t",$F[1],"\t",$F[3],"\n"}' heart_thyroid.sum > heart_thyroid.ratio
 sort heart_thyroid.sum | join - both.80 > heart_thyroid.sum80
 ls
 sort heart_thyroid.sum | join - both.clean80 > heart_thyroid.sum80
 perl -ane '{print $F[0],"\t",log((0.000001+$F[1])/(0.000001+$F[3])),"\t", sqrt($F[2]*$F[2]+$F[4]*$F[4]),"\t",$F[1],"\t",$F[3],"\n"}' heart_thyroid.sum80 > heart_thyroid.ratio80
