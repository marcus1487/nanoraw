g1Dir="signif_testing/group1/"
g2Dir="signif_testing/group2/"
poreModel="template_median68pA.5mers.model"
genomeFn="genomes.e_coli.s_aureus.m_smegmatis.fa"
genomeLocs='"S_aureus:2064835" "S_aureus:2064935"'
strandGenomeLocs='"S_aureus:2064835:-" "S_aureus:2064935"'

printf "********* Testing help commands **********\n"
nanoraw -h

nanoraw genome_resquiggle -h
nanoraw plot_max_coverage -h
nanoraw plot_genome_location -h
nanoraw plot_motif_centered -h

nanoraw plot_max_difference -h
nanoraw plot_most_significant -h
nanoraw plot_motif_with_stats -h

nanoraw plot_correction -h
nanoraw plot_multi_correction -h

nanoraw cluster_most_significant -h
nanoraw plot_kmer -h

nanoraw write_most_significant_fasta -h
nanoraw write_wiggles -h

printf "\n\n********* Testing genome correction scripts **********\n"
nanoraw genome_resquiggle \
        $g1Dir $genomeFn --graphmap-executable ./graphmap \
        --timeout 60 --cpts-limit 100 --normalization-type median \
        --failed-reads-filename testing.signif_group1.failed_read.txt \
        --2d --processes 4 --overwrite
nanoraw genome_resquiggle \
        $g2Dir $genomeFn --graphmap-executable ./graphmap \
        --timeout 60 --cpts-limit 100 --normalization-type median \
        --failed-reads-filename testing.signif_group2.failed_read.txt \
         --overwrite --2d --processes 4

printf "\n\n********* Testing Alternative BWA-MEM Mapper **********\n"
nanoraw genome_resquiggle \
        $g1Dir $genomeFn --bwa-mem-executable ./bwa \
        --timeout 60 --cpts-limit 100 --normalization-type median \
        --corrected-group RawGenomeCorrected_bwamem_000 --overwrite \
        --failed-reads-filename testing.group1.bwamem.failed_read.txt \
        --2d --processes 4

printf "\n\n********* Testing pA normalization **********\n"
nanoraw genome_resquiggle \
        $g1Dir $genomeFn --graphmap-executable ./graphmap \
        --timeout 60 --cpts-limit 100 --normalization-type pA_raw \
        --corrected-group RawGenomeCorrected_pA_raw_000 --overwrite \
        --failed-reads-filename testing.signif_group1.pA.failed_read.txt \
        --2d --processes 4
nanoraw genome_resquiggle \
        $g1Dir $genomeFn --graphmap-executable ./graphmap \
        --timeout 60 --cpts-limit 100 --normalization-type pA \
        --pore-model-filename $poreModel \
        --corrected-group RawGenomeCorrected_pA_000 --overwrite \
        --failed-reads-filename testing.signif_group1.pA.failed_read.txt \
        --2d --processes 4

printf "\n\n********* Testing single sample genome-anchored plotting functions **********\n"
nanoraw plot_max_coverage --fast5-basedirs $g1Dir --2d \
        --num-bases 21 --overplot-threshold 1000 \
        --pdf-filename testing.max_cov.1_samp.pdf
nanoraw plot_max_coverage --fast5-basedirs $g1Dir --2d \
        --num-bases 21 --overplot-threshold 1000 \
        --obs-per-base-filter 99:200 100:5000 \
        --pdf-filename testing.max_cov.1_samp.filt.pdf
nanoraw plot_genome_location --fast5-basedirs $g1Dir \
        --genome-locations $genomeLocs \
        --2d --num-bases 21 --overplot-threshold 1000 \
        --pdf-filename testing.genome_loc.1_samp.pdf
nanoraw plot_motif_centered --fast5-basedirs $g1Dir --motif ATC \
        --genome-fasta $genomeFn --2d \
        --num-bases 21 --overplot-threshold 1000 \
        --pdf-filename testing.motif_centered.1_samp.pdf
nanoraw plot_motif_centered --fast5-basedirs $g1Dir --motif ATC \
        --genome-fasta $genomeFn --2d \
        --num-bases 21 --overplot-threshold 1000 --deepest-coverage \
        --pdf-filename testing.motif_centered.deepest.1_samp.pdf

printf "\n\n********* Testing mutliple sample genome-anchored plotting functions **********\n"
nanoraw plot_max_coverage --fast5-basedirs $g1Dir \
        --fast5-basedirs2 $g2Dir --2d \
        --num-bases 21 --overplot-threshold 1000 \
        --pdf-filename testing.max_cov.2_samp.pdf
nanoraw plot_max_coverage --fast5-basedirs $g1Dir \
        --fast5-basedirs2 $g2Dir --2d \
        --num-bases 21 --overplot-threshold 1000 \
        --obs-per-base-filter 99:200 100:5000 \
        --pdf-filename testing.max_cov.2_samp.filt.pdf
nanoraw plot_genome_location --fast5-basedirs $g1Dir \
        --fast5-basedirs2 $g2Dir \
        --genome-locations $genomeLocs \
        --2d --num-bases 21 --overplot-threshold 1000 \
        --pdf-filename testing.genome_loc.2_samp.pdf
nanoraw plot_motif_centered --fast5-basedirs $g1Dir --motif ATC \
        --genome-fasta $genomeFn \
        --fast5-basedirs2 $g2Dir --2d \
        --num-bases 21 --overplot-threshold 1000 --deepest-coverage \
        --pdf-filename testing.motif_centered.2_samp.pdf

printf "\n\n********* Testing mutliple sample statistical testing genome-anchored plotting functions **********\n"
nanoraw plot_max_difference --fast5-basedirs $g1Dir \
        --fast5-basedirs2 $g2Dir --2d \
        --num-bases 21 --overplot-threshold 1000 \
        --pdf-filename testing.max_diff.pdf
rm testing.significance_values.txt
nanoraw plot_most_significant --fast5-basedirs $g1Dir \
        --fast5-basedirs2 $g2Dir --2d \
        --num-bases 21 --overplot-threshold 1000 \
        --statistics-filename testing.significance_values.txt \
        --pdf-filename testing.most_signif.pdf
nanoraw plot_most_significant --fast5-basedirs $g1Dir \
        --fast5-basedirs2 $g2Dir --2d \
        --num-bases 21 --overplot-threshold 1000 \
        --statistics-filename testing.significance_values.txt \
        --pdf-filename testing.most_signif.re_calc.pdf
nanoraw plot_motif_with_stats --fast5-basedirs $g1Dir \
        --fast5-basedirs2 $g2Dir --motif ATC --2d \
        --overplot-threshold 1000 --test-type mw_utest \
        --pdf-filename testing.motif_w_stats.pdf

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
nanoraw plot_multi_correction --fast5-basedirs $g1Dir \
        --genome-locations $strandGenomeLocs

printf "\n\n********* Testing other plotting commands **********\n"
nanoraw cluster_most_significant --fast5-basedirs $g1Dir \
        --fast5-basedirs2 $g2Dir --2d \
        --genome-fasta $genomeFn --num-regions 100
nanoraw cluster_most_significant --fast5-basedirs $g1Dir \
        --fast5-basedirs2 $g2Dir --2d \
        --genome-fasta $genomeFn --num-regions 100 \
        --r-data-filename testing.cluster_data.RData \
        --statistics-filename testing.significance_values.txt
nanoraw plot_kmer --fast5-basedirs $g1Dir \
        --pdf-filename testing.kmer_dist.median.all_events.pdf
nanoraw plot_kmer --fast5-basedirs $g1Dir --read-mean \
        --num-kmer-threshold 2 \
        --pdf-filename testing.kmer_dist.median.pdf
nanoraw plot_kmer --fast5-basedirs $g1Dir --read-mean \
        --corrected-group RawGenomeCorrected_pA_raw_000 \
        --num-kmer-threshold 2 \
        --pdf-filename testing.kmer_dist.pA_raw.pdf
nanoraw plot_kmer --fast5-basedirs $g1Dir --read-mean \
        --corrected-group RawGenomeCorrected_pA_000 \
        --num-kmer-threshold 2 \
        --pdf-filename testing.kmer_dist.pA.pdf

printf "\n\n********* Testing auxilliary commands **********\n"
nanoraw write_most_significant_fasta --fast5-basedirs $g1Dir \
        --fast5-basedirs2 $g2Dir --2d \
        --genome-fasta $genomeFn
nanoraw write_wiggles --fast5-basedirs $g1Dir \
        --fast5-basedirs2 $g2Dir --2d \
        --wiggle-types coverage signal signal_sd length pvals \
        qvals difference \
        --statistics-filename testing.significance_values.txt
