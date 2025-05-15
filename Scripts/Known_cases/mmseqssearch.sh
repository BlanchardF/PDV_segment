#!/bin/bash

#SBATCH --partition=normal
#SBATCH --time=5:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G

mkdir -p /beegfs/project/horizon/PDV_Segments/Results/Known_cases/run_mmseqs2/Genes_Cotesia_congregata_Genome_Microplitis_demolitor/
mkdir -p /beegfs/project/horizon/PDV_Segments/Results/Known_cases/run_mmseqs2/Genes_Microplitis_demolitor_Genome_Cotesia_congregata/
mkdir -p /beegfs/project/horizon/PDV_Segments/Results/Known_cases/run_mmseqs2/Genes_Cotesia_congregata_Genome_Hyposoter_didymator/
mkdir -p /beegfs/project/horizon/PDV_Segments/Results/Known_cases/run_mmseqs2/Genes_Microplitis_demolitor_Genome_Hyposoter_didymator/
mkdir -p /beegfs/project/horizon/PDV_Segments/Results/Known_cases/run_mmseqs2/Genes_Hyposoter_didymator_Genome_Microplitis_demolitor/
mkdir -p /beegfs/project/horizon/PDV_Segments/Results/Known_cases/run_mmseqs2/Genes_Hyposoter_didymator_Genome_Cotesia_congregata/


/beegfs/project/horizon/bin/miniconda3/bin/mmseqs search \
 /beegfs/project/horizon/PDV_Segments/Results/Known_cases/Genome_db/Cotesia_congregata_mmseqs2_db  \
              /beegfs/project/horizon/PDV_Segments/Results/Known_cases/Genes_db/Genes_Microplitis_demolitor \
              /beegfs/project/horizon/PDV_Segments/Results/Known_cases/run_mmseqs2/Genes_Microplitis_demolitor_Genome_Cotesia_congregata/_result_mmseqs2  \
	     /beegfs/project/horizon/PDV_Segments/Results/Known_cases/run_mmseqs2/Genes_Microplitis_demolitor_Genome_Cotesia_congregata/tmp/  \
                -a -s 7.5 -e 0.0001  --search-type 3 --threads 4 



/beegfs/project/horizon/bin/miniconda3/bin/mmseqs search \
	       /beegfs/project/horizon/PDV_Segments/Results/Known_cases/Genome_db/Hyposoter_didymator_mmseqs2_db  \
              /beegfs/project/horizon/PDV_Segments/Results/Known_cases/Genes_db/Genes_Cotesia_congregata \
              /beegfs/project/horizon/PDV_Segments/Results/Known_cases/run_mmseqs2/Genes_Cotesia_congregata_Genome_Hyposoter_didymator/_result_mmseqs2  \
	     /beegfs/project/horizon/PDV_Segments/Results/Known_cases/run_mmseqs2/Genes_Cotesia_congregata_Genome_Hyposoter_didymator/tmp/  \
                -a -s 7.5 -e 0.0001  --search-type 3 --threads 4 

/beegfs/project/horizon/bin/miniconda3/bin/mmseqs search \
              /beegfs/project/horizon/PDV_Segments/Results/Known_cases/Genome_db/Hyposoter_didymator_mmseqs2_db  \
              /beegfs/project/horizon/PDV_Segments/Results/Known_cases/Genes_db/Genes_Microplitis_demolitor \
              /beegfs/project/horizon/PDV_Segments/Results/Known_cases/run_mmseqs2/Genes_Microplitis_demolitor_Genome_Hyposoter_didymator/_result_mmseqs2  \
	     /beegfs/project/horizon/PDV_Segments/Results/Known_cases/run_mmseqs2/Genes_Microplitis_demolitor_Genome_Hyposoter_didymator/tmp/  \
                -a -s 7.5 -e 0.0001  --search-type 3 --threads 4 



/beegfs/project/horizon/bin/miniconda3/bin/mmseqs search \
              /beegfs/project/horizon/PDV_Segments/Results/Known_cases/Genome_db/Microplitis_demolitor_mmseqs2_db  \
              /beegfs/project/horizon/PDV_Segments/Results/Known_cases/Genes_db/Genes_Hyposoter_didymator \
              /beegfs/project/horizon/PDV_Segments/Results/Known_cases/run_mmseqs2/Genes_Hyposoter_didymator_Genome_Microplitis_demolitor/_result_mmseqs2  \
	     /beegfs/project/horizon/PDV_Segments/Results/Known_cases/run_mmseqs2/Genes_Hyposoter_didymator_Genome_Microplitis_demolitor/tmp/  \
                -a -s 7.5 -e 0.0001  --search-type 3 --threads 4 


/beegfs/project/horizon/bin/miniconda3/bin/mmseqs search \
              /beegfs/project/horizon/PDV_Segments/Results/Known_cases/Genome_db/Cotesia_congregata_mmseqs2_db  \
              /beegfs/project/horizon/PDV_Segments/Results/Known_cases/Genes_db/Genes_Hyposoter_didymator \
              /beegfs/project/horizon/PDV_Segments/Results/Known_cases/run_mmseqs2/Genes_Hyposoter_didymator_Genome_Cotesia_congregata/_result_mmseqs2  \
	     /beegfs/project/horizon/PDV_Segments/Results/Known_cases/run_mmseqs2/Genes_Hyposoter_didymator_Genome_Cotesia_congregata/tmp/  \
                -a -s 7.5 -e 0.0001  --search-type 3 --threads 4 










/beegfs/project/horizon/bin/miniconda3/bin/mmseqs convertalis \
	--format-output 'query,qlen,tlen,target,pident,alnlen,mismatch,gapopen,qstart,qend,tstart,tend,evalue,bits,qaln,tcov,qcov' \
              /beegfs/project/horizon/PDV_Segments/Results/Known_cases/Genome_db/Microplitis_demolitor_mmseqs2_db  \
              /beegfs/project/horizon/PDV_Segments/Results/Known_cases/Genes_db/Genes_Cotesia_congregata \
              /beegfs/project/horizon/PDV_Segments/Results/Known_cases/run_mmseqs2/Genes_Cotesia_congregata_Genome_Microplitis_demolitor/_result_mmseqs2  \
	     /beegfs/project/horizon/PDV_Segments/Results/Known_cases/run_mmseqs2/Genes_Cotesia_congregata_Genome_Microplitis_demolitor/_result_mmseqs2.m8  \
              --search-type 3

/beegfs/project/horizon/bin/miniconda3/bin/mmseqs convertalis \
	--format-output 'query,qlen,tlen,target,pident,alnlen,mismatch,gapopen,qstart,qend,tstart,tend,evalue,bits,qaln,tcov,qcov' \
              /beegfs/project/horizon/PDV_Segments/Results/Known_cases/Genome_db/Cotesia_congregata_mmseqs2_db  \
              /beegfs/project/horizon/PDV_Segments/Results/Known_cases/Genes_db/Genes_Microplitis_demolitor \
              /beegfs/project/horizon/PDV_Segments/Results/Known_cases/run_mmseqs2/Genes_Microplitis_demolitor_Genome_Cotesia_congregata/_result_mmseqs2  \
	     /beegfs/project/horizon/PDV_Segments/Results/Known_cases/run_mmseqs2/Genes_Microplitis_demolitor_Genome_Cotesia_congregata/_result_mmseqs2.m8  \
                --search-type 3





/beegfs/project/horizon/bin/miniconda3/bin/mmseqs convertalis \
	--format-output 'query,qlen,tlen,target,pident,alnlen,mismatch,gapopen,qstart,qend,tstart,tend,evalue,bits,qaln,tcov,qcov' \
              /beegfs/project/horizon/PDV_Segments/Results/Known_cases/Genome_db/Hyposoter_didymator_mmseqs2_db  \
              /beegfs/project/horizon/PDV_Segments/Results/Known_cases/Genes_db/Genes_Cotesia_congregata \
              /beegfs/project/horizon/PDV_Segments/Results/Known_cases/run_mmseqs2/Genes_Cotesia_congregata_Genome_Hyposoter_didymator/_result_mmseqs2  \
	     /beegfs/project/horizon/PDV_Segments/Results/Known_cases/run_mmseqs2/Genes_Cotesia_congregata_Genome_Hyposoter_didymator/_result_mmseqs2.m8  \
                --search-type 3

/beegfs/project/horizon/bin/miniconda3/bin/mmseqs convertalis \
	--format-output 'query,qlen,tlen,target,pident,alnlen,mismatch,gapopen,qstart,qend,tstart,tend,evalue,bits,qaln,tcov,qcov' \
              /beegfs/project/horizon/PDV_Segments/Results/Known_cases/Genome_db/Hyposoter_didymator_mmseqs2_db  \
              /beegfs/project/horizon/PDV_Segments/Results/Known_cases/Genes_db/Genes_Microplitis_demolitor \
              /beegfs/project/horizon/PDV_Segments/Results/Known_cases/run_mmseqs2/Genes_Microplitis_demolitor_Genome_Hyposoter_didymator/_result_mmseqs2  \
	     /beegfs/project/horizon/PDV_Segments/Results/Known_cases/run_mmseqs2/Genes_Microplitis_demolitor_Genome_Hyposoter_didymator/_result_mmseqs2.m8  \
               --search-type 3



//beegfs/project/horizon/bin/miniconda3/bin/mmseqs convertalis\
	--format-output 'query,qlen,tlen,target,pident,alnlen,mismatch,gapopen,qstart,qend,tstart,tend,evalue,bits,qaln,tcov,qcov' \
              /beegfs/project/horizon/PDV_Segments/Results/Known_cases/Genome_db/Microplitis_demolitor_mmseqs2_db  \
              /beegfs/project/horizon/PDV_Segments/Results/Known_cases/Genes_db/Genes_Hyposoter_didymator \
              /beegfs/project/horizon/PDV_Segments/Results/Known_cases/run_mmseqs2/Genes_Hyposoter_didymator_Genome_Microplitis_demolitor/_result_mmseqs2  \
	     /beegfs/project/horizon/PDV_Segments/Results/Known_cases/run_mmseqs2/Genes_Hyposoter_didymator_Genome_Microplitis_demolitor/_result_mmseqs2.m8  \
               --search-type 3


/beegfs/project/horizon/bin/miniconda3/bin/mmseqs convertalis \
	--format-output 'query,qlen,tlen,target,pident,alnlen,mismatch,gapopen,qstart,qend,tstart,tend,evalue,bits,qaln,tcov,qcov' \
              /beegfs/project/horizon/PDV_Segments/Results/Known_cases/Genome_db/Cotesia_congregata_mmseqs2_db  \
              /beegfs/project/horizon/PDV_Segments/Results/Known_cases/Genes_db/Genes_Hyposoter_didymator \
              /beegfs/project/horizon/PDV_Segments/Results/Known_cases/run_mmseqs2/Genes_Hyposoter_didymator_Genome_Cotesia_congregata/_result_mmseqs2  \
	     /beegfs/project/horizon/PDV_Segments/Results/Known_cases/run_mmseqs2/Genes_Hyposoter_didymator_Genome_Cotesia_congregata/_result_mmseqs2.m8  \
              --search-type 3



