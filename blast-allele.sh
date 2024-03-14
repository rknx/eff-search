__ts() { echo -e [`date "+%m/%d %H:%M:%S"`]: $1; }
export -f __ts

ncpu=8

# Effector AA file
# ls eff_prot_xg.faa

# File with inks for downlaoding PGAP genbank files
# ls ncbi_accn.txt

## Get PGAP files --------------------------------------------------------

mkdir annots
__get_genome(){
    strain=$(basename ${1%%_genome*} | tr '[:lower:]' '[:upper:]')
    wget -O annots/$strain.gb "$1"
    cat annots/$strain.gb | ./gb2faa.py > annots/"$strain"_cds.ffn
}

export -f __get_genome

parallel --will-cite __get_genome :::: ncbi_accn.txt

## Creating separate query list ------------------------------------------------
# cp eff_prot_xg.faa 
awk -vRS=">" -vORS="" -vFS="\n" -vOFS="" '
    NR>1{
        split($1, a, "_");
        $1="";
        b=">"a[1]"\n"$0;
        system("mkdir -p query/"a[1]);
        system("echo -e \"" b "\" > query/"a[1]"/prot.fa")
    }
' eff_prot_xg.faa

## Creating BLAST databases ----------------------------------------------------

__ts "Making BLAST databases."

# Load modules
module -q reset; module load parallel ncbi_blast

# Get list of genomes
genlist=( $(
	find annots/ -name "*_cds.ffn" | \
		sed 's/^annots\///; s/_cds\.ffn$//' | sort
) )
__ts "${#genlist[@]} genomes found."

__make_blast(){

	# Creating cds BLAST database
    [[ ! -s db/$1/cds.pto ]] && \
		__ts "Creating BLAST database for $1 genes." && \
        mkdir -p db/$1 && \
        ln -s $(pwd)/annots/"$1"_cds.ffn db/$1/cds && \
        makeblastdb -dbtype nucl \
            -in db/$1/cds \
            -title $1 \
            -out db/$1/cds || \
		__ts "BLAST database is present for $1 genes. Skipping."

}

# Make function global
export -f __make_blast

# Run the polishing function
parallel --will-cite -j $(( $ncpu > 8 ? 8 : $ncpu )) __make_blast ::: ${genlist[@]}

__ts "Finished creating BLAST database."

## Doing BLAST -----------------------------------------------------------------

__ts "Blasting sequences."

# Load modules
module -q reset; module load parallel ncbi_blast samtools

# Get list of queries
qrylist=( $(
	find query/ -name "prot.fa" | \
		sed 's/^query\///; s/\/prot\.fa$//' | sort
) )
__ts "${#qrylist[@]} queries found."

# Get list of blast databases
dblist=( $(
	find db/ -name "*.nto" | \
		sed 's/^db\///; s/\/cds\.nto$//' | \
        sort
) )
__ts "${#dblist[@]} databases found."

__blast_query(){

	# Blasting queries to databases
    [[ ! -s out/$2/$1.seq ]] && \
		__ts "BLASTing $2 to $1." && \
        mkdir -p out/$2 && \
        tblastn -outfmt 6 -qcov_hsp_perc 60 \
            -query query/$2/prot.fa \
            -db db/$1/cds | \
            awk -vdb="$1" '$3>60 {
                "samtools faidx -n0 db/" db "/cds " $2 " | grep -v \">\"" | getline seq
                split($2,stat,":")
                print stat[2] seq
            }' > \
            out/$2/$1.seq || \
		__ts "$2 BLAST output already exists for $1. Skipping."

}

# Make function global
export -f __blast_query

# Run the BLAST function
parallel --will-cite -j $(( $ncpu > 8 ? 8 : $ncpu )) __blast_query \
    ::: ${dblist[@]} ::: ${qrylist[@]}

__ts "Finished BLASTing."

## Getting allele list ---------------------------------------------------------

__ts "Gathering alleles."

# Load modules
module -q reset; module load parallel

# Get list of queries
qrylist=( $(
	find query/ -name "prot.fa" | \
		sed 's/^query\///; s/\/prot\.fa$//' | sort
) )
__ts "${#qrylist[@]} queries found."

# Function to gather alleles
__allele_list(){

	# 
    [[ ! -s out/$1/alleles.tab ]] && \
		__ts "Gathering outputs." && \
        cat out/$1/*.seq | \
            sort -u | \
            awk -vOFS="\t" '{print NR, $0}' > \
            out/$1/alleles.tab || \
		__ts "BLAST outputs are already gathered. Skipping."
}

# Make function global
export -f __allele_list

# Run the polishing function
parallel --will-cite -k -j $(( $ncpu > 8 ? 8 : $ncpu )) __allele_list ::: ${qrylist[@]}

__ts "Finished gathering alleles."

## Assigning allele number -------------------------------------------------------

__ts "Assigning number to alleles."

# Load modules
module -q reset; module load parallel

# Get list of queries
qrylist=( $(
	find out/ -name "alleles.tab" | \
		sed 's/^out\///; s/\/alleles\.tab$//' | sort
) )
__ts "${#qrylist[@]} queries found."

# Get list of blast databases
dblist=( $(
	find  -name "*.nto" | \
		sed 's/^db\///; s/\/cds\.nto$//' | \
        sort
) )
__ts "${#dblist[@]} databases found."

__allele_number(){

	# Assigning allele number
    ! grep -qs $1 output/$2.allele && \
		__ts "Begin assigning allele number." && \
		mkdir -p output && \
        awk -vstr="$1" -vqry="$2" -vOFS="\t" '
            BEGIN{allele="0"}
            NR==FNR{a[$2]=$1; next}
            $1 in a{allele=allele"/"a[$1]}
            END{gsub("0/","",allele); print qry, str, allele}
        ' out/$2/alleles.tab out/$2/$1.seq >> \
        output/$2.allele || \
		__ts "Allele number is already assigned. Skipping."
}

# Make function global
export -f __allele_number

# Run the polishing function
parallel -k -j $(( $ncpu > 8 ? 8 : $ncpu )) __allele_number ::: ${dblist[@]} ::: ${qrylist[@]}

__ts "Allele number assignment complete."

## Combining results -------------------------------------------------------

__ts "Combining allele tables."

# Load modules
module -q reset

# Get list of queries
tablist=( $(
	find output/ -name "*.allele" | sort
) )
__ts "${#tablist[@]} queries found." 

__allele_number(){

	# 
    [[ ! -s allele.table ]] && \
		__ts "Generating combined allele table." && \
        cat "$@" | \
        awk -vOFS="\t" '
            {row[$1]=1; col[$2]=1; val[$1" "$2]=$3}
            END {

                m=asorti(col, icol)
                n=asorti(row, irow)

                printf "%s", "Effector"
                for (i=1; i<=m; i++) { split(icol[i], a, "_"); printf "\t%s",a[1] }; printf "\n"

                for (j=1; j<=n; j++) {
                    printf "%s", irow[j]
                    for (i=1; i<=m; i++) { printf "\t%s",val[irow[j]" "icol[i]] }
                    printf "\n"
                }
            }
        ' > allele.table || \
		__ts "Combined allele table exists. Skipping."
}

# Make function global
export -f __allele_number

# Run the polishing function
__allele_number ${tablist[*]}

__ts "All done."
