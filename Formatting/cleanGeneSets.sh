#Script to re-format a gmx file with multiple gene sets
#Usage: bash cleanGeneSets.sh batchGeneSetFile
#Usage Ex: bash cleanGeneSets.sh ~/PfrenderLab/PA42_v4.1/genesets.gmx

#Determine the number of gene sets
numSets=$(head -1 ~/PfrenderLab/PA42_v4.1/genesets.gmx | sed "s/[[:space:]]/\n/g" | sed '/^$/d' | wc -l | sed "s/[[:space:]]//g")

#Loop over each gene set and output to separate files
for i in $( eval echo {1..$numSets} ); do
	setName=$(cat ~/PfrenderLab/PA42_v4.1/genesets.gmx | cut -f $i | head -1)
	cat ~/PfrenderLab/PA42_v4.1/genesets.gmx | cut -f $i | sed '/^$/d' | tail -n+3 > ~/PfrenderLab/PA42_v4.1/geneSet"$i"_"$setName".txt
done