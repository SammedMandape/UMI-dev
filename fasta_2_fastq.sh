# fasta to fastq, only works if you already have sequence per record in one line. 
awk 'BEGIN {
		RS=">";
		FS="\n"}
		NR>1
		{print "@"$1"\n"$2"\n+"; 
		for(c=0;c<length($2);c++) 
		printf "~"; 
		printf "\n"}' 