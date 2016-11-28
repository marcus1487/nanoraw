g1Dir="signif_testing/group1/"
g2Dir="signif_testing/group2/"
genomeFn="genomes.e_coli.s_aureus.m_smegmatis.fa"
genomeLocs='"S_aureus:2064835" "S_aureus:2064935"'

printf "********* Testing help commands **********\n"
nanoraw -h

nanoraw genome_resquiggle -h
nanoraw plot_max_coverage -h
nanoraw plot_genome_location -h
nanoraw plot_kmer_centered -h

nanoraw plot_max_difference -h
nanoraw plot_most_significant -h
nanoraw plot_kmer_with_stats -h

nanoraw plot_correction -h
nanoraw plot_multi_correction -h

nanoraw cluster_most_significant -h
nanoraw plot_kmer -h

nanoraw write_most_significant -h
nanoraw write_wiggle -h

printf "\n\n********* Testing genome correction scripts **********\n"
nanoraw genome_resquiggle \
        $g1Dir ./graphmap $genomeFn \
        --timeout 60 --cpts-limit 100 --normalization-type median --overwrite \
        --failed-reads-filename testing.signif_group1.failed_read.txt \
        --2d --processes 4
nanoraw genome_resquiggle \
        $g2Dir ./graphmap $genomeFn \
        --timeout 60 --cpts-limit 100 --normalization-type median --overwrite \
        --failed-reads-filename testing.signif_group2.failed_read.txt \
        --2d --processes 4

printf "\n\n********* Testing pA normalization **********\n"
nanoraw genome_resquiggle \
        $g1Dir ./graphmap $genomeFn \
        --timeout 60 --cpts-limit 100 --normalization-type ont \
        --corrected-group RawGenomeCorrected_pA_000 --overwrite \
        --failed-reads-filename testing.signif_group1.pA.failed_read.txt \
        --2d --processes 4

printf "\n\n********* Testing single sample genome-anchored plotting functions **********\n"
nanoraw plot_max_coverage --fast5-basedirs $g1Dir --2d \
        --num-bases 21 --overplot-threshold 1000
nanoraw plot_max_coverage --fast5-basedirs $g1Dir --2d \
        --num-bases 21 --overplot-threshold 1000 \
        --obs-per-base-filter 99:200 100:5000
nanoraw plot_genome_location --fast5-basedirs $g1Dir \
        --genome-locations $genomeLocs \
        --2d --num-bases 21 --overplot-threshold 1000
nanoraw plot_kmer_centered --fast5-basedirs $g1Dir --kmer ATC \
        --genome-fasta $genomeFn --2d \
        --num-bases 21 --overplot-threshold 1000
nanoraw plot_kmer_centered --fast5-basedirs $g1Dir --kmer ATC \
        --genome-fasta $genomeFn --2d \
        --num-bases 21 --overplot-threshold 1000 --deepest-coverage

printf "\n\n********* Testing mutliple sample genome-anchored plotting functions **********\n"
nanoraw plot_max_coverage --fast5-basedirs $g1Dir \
        --fast5-basedirs2 $g2Dir --2d \
        --num-bases 21 --overplot-threshold 1000
nanoraw plot_max_coverage --fast5-basedirs $g1Dir \
        --fast5-basedirs2 $g2Dir --2d \
        --num-bases 21 --overplot-threshold 1000 --obs-per-base-filter 99:200 100:5000
nanoraw plot_genome_location --fast5-basedirs $g1Dir \
        --fast5-basedirs2 $g2Dir \
        --genome-locations $genomeLocs \
        --2d --num-bases 21 --overplot-threshold 1000
nanoraw plot_kmer_centered --fast5-basedirs $g1Dir --kmer ATC \
        --genome-fasta $genomeFn \
        --fast5-basedirs2 $g2Dir --2d \
        --num-bases 21 --overplot-threshold 1000 --deepest-coverage

printf "\n\n********* Testing mutliple sample statistical testing genome-anchored plotting functions **********\n"
nanoraw plot_max_difference --fast5-basedirs $g1Dir \
        --fast5-basedirs2 $g2Dir --2d \
        --num-bases 21 --overplot-threshold 1000
nanoraw plot_most_significant --fast5-basedirs $g1Dir \
        --fast5-basedirs2 $g2Dir --2d \
        --num-bases 21 --overplot-threshold 1000
nanoraw plot_kmer_with_stats --fast5-basedirs $g1Dir \
        --fast5-basedirs2 $g2Dir --motif ATC --2d \
        --overplot-threshold 1000 --test-type mw_utest \
        --genome-fasta $genomeFn

printf "\n\n********* Testing overplotting options **********\n"
nanoraw plot_max_coverage --fast5-basedirs $g1Dir \
        --fast5-basedirs2 $g2Dir --2d \
        --num-bases 21 --overplot-threshold 20 --overplot-type Downsample \
        --pdf-filename Nanopore_read_coverage.max_coverage.Downsample.pdf
nanoraw plot_max_coverage --fast5-basedirs $g1Dir \
        --fast5-basedirs2 $g2Dir --2d \
        --num-bases 21 --overplot-threshold 20 --overplot-type Boxplot \
        --pdf-filename Nanopore_read_coverage.max_coverage.Boxplot.pdf
nanoraw plot_max_coverage --fast5-basedirs $g1Dir \
        --fast5-basedirs2 $g2Dir --2d \
        --num-bases 21 --overplot-threshold 20 --overplot-type Quantile \
        --pdf-filename Nanopore_read_coverage.max_coverage.Quantile.pdf
nanoraw plot_max_coverage --fast5-basedirs $g1Dir \
        --fast5-basedirs2 $g2Dir --2d \
        --num-bases 21 --overplot-threshold 20 --overplot-type Violin \
        --pdf-filename Nanopore_read_coverage.max_coverage.Violin.pdf

printf "\n\n********* Testing correction plotting commands **********\n"
nanoraw plot_correction --fast5-basedirs $g1Dir --region-type random
nanoraw plot_multi_correction --fast5-basedirs $g1Dir

printf "\n\n********* Testing other plotting commands **********\n"
nanoraw cluster_most_significant --fast5-basedirs $g1Dir \
        --fast5-basedirs2 $g2Dir --2d \
        --genome-fasta $genomeFn --num-regions 100
nanoraw cluster_most_significant --fast5-basedirs $g1Dir \
        --fast5-basedirs2 $g2Dir --2d \
        --genome-fasta $genomeFn --num-regions 100 \
        --r-data-filename testing.cluster_data.RData \
        --statistics-filename testing.significance_values.txt
nanoraw plot_kmer --fast5-basedirs $g1Dir
nanoraw plot_kmer --fast5-basedirs $g1Dir --read-mean

printf "\n\n********* Testing auxilliary commands **********\n"
nanoraw write_most_significant --fast5-basedirs $g1Dir \
        --fast5-basedirs2 $g2Dir --2d \
        --genome-fasta $genomeFn
nanoraw write_wiggle --fast5-basedirs $g1Dir
