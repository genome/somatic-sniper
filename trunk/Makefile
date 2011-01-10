all:
	gcc -g -Wall -O2 -I/gsc/pkg/bio/samtools/samtools-0.1.6/ sniper_pileup.c sniper_maqcns.c somatic_sniper.c -o bam-somaticsniper -Wl,-Bstatic -lm -L/gsc/pkg/bio/samtools/samtools-0.1.6/ -lbam -lz -Wl,-Bdynamic
	#gcc -g -Wall -O2 -I/gsc/pkg/bio/samtools/samtools-0.1.6/ sniper_pileup.c sniper_maqcns.c somatic_sniper.c -o bam-somaticsniper -lm -lz -L/gsc/pkg/bio/samtools/samtools-0.1.6/ -lbam
