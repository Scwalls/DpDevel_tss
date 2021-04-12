#setwd("/scratch/tswenty/GoSTRIPES_sing/STRIPES/tsr_Fs")
setwd("/scratch/scwalls/STRIPEseq_pipelines/GoSTRIPES_sing/STRIPES/figures/4_7_21")

library(TSRchitect)
library(TSRexploreR)
library(GenomicRanges)
library(DESeq2)
library(ChIPseeker)
library(ggseqlogo)
library(ggplot2)

#load("/scratch/tswenty/GoSTRIPES_sing/STRIPES/tsr_Es/PdSTRIPE_complete.RData")
load("/scratch/scwalls/STRIPEseq_pipelines/GoSTRIPES_sing/STRIPES/tsr_Fs/PdSTRIPE_complete.RData")

#creating the annotation and assembly files
#update both paths below
Dp.annot <- "/data/LynchLabCME/Daphnia/DaphniaTissues/GoSTRIPES_sing_DpTissues/STRIPES/DpGENOME/PA42.4.0.gff"
Dp.assembly <- "/data/LynchLabCME/Daphnia/DaphniaTissues/GoSTRIPES_sing_DpTissues/STRIPES/DpGENOME/PA42.4.1.fasta"

# Generate sample sheet
sample_sheet <- data.frame(
  sample_name=stringr::str_glue("DpDevel_{seq_len(3)}"),
  file_1=NA, file_2=NA,
  condition=rep("TimePoint_C", 3)
)

#writing the tss files to the workspace
tss.1 <- PdSTRIPE@tssCountData[[1]]
tss.2 <- PdSTRIPE@tssCountData[[2]]
tss.3 <- PdSTRIPE@tssCountData[[3]]
#tss.4 <- PdSTRIPE@tssCountData[[4]]

colnames(tss.1) <- c("seq","TSS", "strand", "score")
colnames(tss.2) <- c("seq","TSS", "strand", "score")
colnames(tss.3) <- c("seq","TSS", "strand", "score")
#colnames(tss.4) <- c("seq","TSS", "strand", "score")

#making granges files from tss data frames
tss.1.gr <- makeGRangesFromDataFrame(tss.1,
                                     keep.extra.columns = TRUE,
                                     seqnames.field="seq",
                                     start.field="TSS",
                                     end.field="TSS",
                                     strand.field = "strand"
)

tss.2.gr <- makeGRangesFromDataFrame(tss.2,
                                     keep.extra.columns = TRUE,
                                     seqnames.field="seq",
                                     start.field="TSS",
                                     end.field="TSS",
                                     strand.field = "strand"
)

tss.3.gr <- makeGRangesFromDataFrame(tss.3,
                                     keep.extra.columns = TRUE,
                                     seqnames.field="seq",
                                     start.field="TSS",
                                     end.field="TSS",
                                     strand.field = "strand"
)

#tss.4.gr <- makeGRangesFromDataFrame(tss.4,
#                                     keep.extra.columns = TRUE,
#                                     seqnames.field="seq",
#                                     start.field="TSS",
#                                     end.field="TSS",
#                                     strand.field = "strand"
#)

Dp.tss <- list(tss.1.gr, tss.2.gr, tss.3.gr)
names(Dp.tss) <- c("F1", "F2", "F3")

#Creating the TSR explorer object
exp <- tsr_explorer(TSSs=Dp.tss, 
                    genome_annotation=Dp.annot, genome_assembly=Dp.assembly,
                    sample_sheet = sample_sheet
)

#Initial TSS processing

exp <- format_counts(exp, data_type="tss")

#Normalize TSSs
exp <- normalize_counts(exp, data_type = "tss", method = "DESeq2")

## TSS annotation
exp <- annotate_features(exp, data_type = "tss", feature_type="transcript")

##### Sequence analysis
## creating a truncated object for sequence analysis
## some intervals are too close to the edges of short scaffolds, so this was my workaround

plot_threshold_exploration(exp, samples="F1", point_size=1) + scale_color_viridis_c()
ggsave(file="threshold_F1.png")

plot_threshold_exploration(exp, samples="F2", point_size=1) + scale_color_viridis_c()
ggsave(file="threshold_F2.png")

plot_threshold_exploration(exp, samples="F3", point_size=1) + scale_color_viridis_c()
ggsave(file="threshold_F3.png")

#plot_threshold_exploration(exp, samples="As4", point_size=1) + scale_color_viridis_c()
#ggsave(file="threshold_Es4.png")

## Correlation plot: replicates
plot_correlation(
  exp, data_type="tss",
  use_normalized=TRUE, font_size=12,
  heatmap_colors=viridis::viridis(100)
)

ggsave(file="Fs_correlation_matrix.png")

### Genomic distribution analysis:
#### To update with new UTR-less annotation

plot_genomic_distribution(exp, data_type="tss", samples=c("F1", "F2","F3")) +
  scale_fill_viridis_d(direction=-1, name="Annotation")

ggsave(file="genomic_distribution_Fs.png")

### promoter fraction plot
#### To update with new UTR-less annotation
plot_detected_features(exp, data_type="tss", samples=c("F1", "F2","F3")) +
  scale_fill_viridis_d(direction=-1)

ggsave(file="promoter_fraction_Fs.png")

### Density plot
plot_density(exp, data_type="tss", samples="F1")
ggsave(file="TSS_density_CDS_F1.png")

plot_density(exp, data_type="tss", samples="F2")
ggsave(file="TSS_density_CDS_F2.png")

plot_density(exp, data_type="tss", samples="F3")
ggsave(file="TSS_density_CDS_F3.png")

#plot_density(exp, data_type="tss", samples="As4")
#ggsave(file="TSS_density_CDS_Es4.png")

### TSS pileup heatmap
## need to update annotation/colour
plot_heatmap(
  exp, data_type="tss", samples="F1",
  upstream=250, downstream=250,
  use_normalized=TRUE,
  rasterize=TRUE, raster_dpi=150
)

ggsave(file="TSS_pileup_F1.png")

plot_heatmap(
  exp, data_type="tss", samples="F2",
  upstream=250, downstream=250,
  use_normalized=TRUE,
  rasterize=TRUE, raster_dpi=150
)

ggsave(file="TSS_pileup_Fs2.png")

plot_heatmap(
  exp, data_type="tss", samples="F3",
  upstream=250, downstream=250,
  use_normalized=TRUE,
  rasterize=TRUE, raster_dpi=150
)

ggsave(file="TSS_pileup_F3.png")

#plot_heatmap(
#  exp, data_type="tss", samples="As4",
#  upstream=250, downstream=250,
#  use_normalized=TRUE,
#  rasterize=TRUE, raster_dpi=150
#)

#ggsave(file="TSS_pileup_Es4.png")

### Sequence logo analysis- all three replicates
plot_sequence_logo(exp, samples="F1")
ggsave(file="F1_seq_logo.png")

plot_sequence_logo(exp, samples="F2")
ggsave(file="F2_seq_logo.png")

plot_sequence_logo(exp, samples="F3")
ggsave(file="F3_seq_logo.png")

#plot_sequence_logo(exp, samples="As4")
#ggsave(file="As4_seq_logo.png")

### Dinucleotide frequency- all three replicates
#### TODO: repeat As1 with truncated tss data (get code from Taylor)
plot_dinucleotide_frequencies(exp, samples="F1") +
  scale_fill_viridis_c()
ggsave(file="dinucleotide_frequency_F1.png")

plot_dinucleotide_frequencies(exp, samples="F2") +
  scale_fill_viridis_c()
ggsave(file="dinucleotide_frequency_F2.png")

plot_dinucleotide_frequencies(exp, samples="F3") +
  scale_fill_viridis_c()
ggsave(file="dinucleotide_frequency_F3.png")

#plot_dinucleotide_frequencies(exp, samples="As4") +
#  scale_fill_viridis_c()
#ggsave(file="dinucleotide_frequency_Es4.png")

### TSS Sequence colour map

plot_sequence_colormap(exp, samples="F1", rasterize=TRUE)
ggsave(file="sequence_colormap_F1.png")
plot_sequence_colormap(exp, samples="F2", rasterize=TRUE)
ggsave(file="sequence_colormap_F2.png")
plot_sequence_colormap(exp, samples="F3", rasterize=TRUE)
ggsave(file="sequence_colormap_F3.png")
#plot_sequence_colormap(exp, samples="As4", rasterize=TRUE)
#ggsave(file="sequence_colormap_Es4.png")

### identify TSRs using clustering
exp <- tss_clustering(exp, threshold=3, n_samples=3, max_distance = 25)

# Associate TSSs with TSRs
exp <- associate_with_tsr(exp)

# Annotate TSRs
exp <- annotate_features(exp, data_type="tsr", upstream=250, downstream=100,
                         feature_type="transcript")

# Mark dominant TSS per TSR
exp <- mark_dominant(exp, data_type="tss")

# Calculate TSR metrics
exp <- tsr_metrics(exp)

###### Add tsr analysis from Standard Analysis documentation

### plot selected tsr metrics
#### TODO: to re-generate using custom script
plot_tsr_metric(exp, tsr_metrics=c("score", "width"), log2_transform=TRUE, samples="all")
ggsave(file="plot_tsr_metrics_Fs.png")

# Dinucleotide motifs by TSR shape
#### TODO: examine As2 and As3, which look odd
plot_sequence_logo(exp, dominant=TRUE, samples="all",
                   data_conditions=conditionals(data_grouping=shape_class)
)
ggsave(file="dinucl_motif_plot_shape_Fs.png")

plot_sequence_logo(exp, dominant=TRUE, samples="F1",
                   data_conditions=conditionals(data_grouping=shape_class)
)

ggsave(file="dinucl_motif_plot_shape_F1.png")

plot_sequence_logo(exp, dominant=TRUE, samples="F2",
                   data_conditions=conditionals(data_grouping=shape_class)
)

ggsave(file="dinucl_motif_plot_shape_F2.png")

plot_sequence_logo(exp, dominant=TRUE, samples="F3",
                   data_conditions=conditionals(data_grouping=shape_class)
)

ggsave(file="dinucl_motif_plot_shape_F3.png")

#plot_sequence_logo(exp, dominant=TRUE, samples="As4",
#                   data_conditions=conditionals(data_grouping=shape_class)
#)

#ggsave(file="dinucl_motif_plot_shape_Es4.png")

### Plot selected gene track
#### TODO: Taylor check with Bob about rendering more clearly
gene_tracks(
  exp, feature_name="dp_gene7125", use_normalized = TRUE, promoter_only=TRUE,
  samples=c(TSS="F1", TSR="F1")
)

ggsave(file="gene_track_F1_gene7125.png")

