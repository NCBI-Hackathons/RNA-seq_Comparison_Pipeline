#This script process the RNAseq tissue type and project name you input
#fastqdump to fasta
#alignment by kallisto
#combine the abundance.txt
#you need to combine all the abundance files, here we saved it as count.all.6.tissue
#Then you make the MDSplot and heatmap

########### input region
tissue_type='heart'
project_name='ERP003613'
###########

#grep all rows that contains the input criteria
grep ${project_name} sra_150715.tab| grep $tissue_type > ${tissue_type}.txt
awk -F'\t' '{print $3}' ${tissue_type}.txt > run.id.txt

mkdir Result_${project_name}
cd Result_${project_name}

mkdir Result_${tissue_type}
cd Result_${tissue_type}

mv ../../run.id.txt .
mv ../../${tissue_type}.txt .
cp ../../python.fastqdump.kalisto.gen.py .

#generate the fastadump and kalisto command
python python.fastqdump.kalisto.gen.py
source fastqdump.kallisto.sh

#let's save some space!
rm *.fasta

#Change the head line of the count.all files
sed -e '1s/\/home\/ubuntu\/Folder_Chen\/Indentify_samples\/Result_ERP003613\/Result_//g' count.all > count.all.new
sed -e '1s/\/output//g' count.all.new > count.all.new.new
sed -e '1s/\/abundance.txt//g' count.all.new.new > count.all
rm *.new
sed '1d' count.all > temp
mv temp > count.all

#MDS plot and heatmap of top 200 std transcript
Rscript command_MDSplot_heatmap.R

