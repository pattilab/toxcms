# install package from local file path
file_path <- "/Users/Lingjue/Documents/Research/PattiLab/Bioinformatics/Rcode/toxcms_1.0.5.tar.gz" #Mac OS
file_path <- "E:/Box Sync/Box Sync/PattiLab/Bioinformatics/Rcode/toxcms_1.0.5.tar.gz" #Windows
install.packages(pkgs = file_path, repos = NULL, type="source")

# OR install package from github
devtools::install_github("pattilab/toxcms",ref = "master",dependencies = TRUE)

# load required packages
library(toxcms)
library(data.table)

# upload example dataset from toxcms
data_path <- system.file("extdata", "example.csv", package="toxcms")
feature <- data.table(read.csv(data_path))

# Step 1 statistical analysis
etom_dosestat <- calcdosestat(Feature = feature, Dose_Levels = c("_0uM","_10uM","_50uM","_200uM"), multicomp = "none",
                              p.adjust.method = "none",projectName = "etomoxir_dataset")

# Step 2.1 monotonic trend filtering
etom_drreport_mono <- trendfilter(etom_dosestat, pval_cutoff = 0.05, pval_thres = 1, anova_cutoff = NULL, trend = "mono",
                                           relChange_cutoff = 0.05, export = FALSE)
# Step 3.1 ED50 modeling and estimation
etom_drreport_fit = fitdrc(DoseResponse_report=etom_drreport_mono, Dose_values=c(0, 10, 50, 200), ED=0.5, export = TRUE,
                           mz_tag = "mzmed", rt_tag = "rtmed", plot=FALSE)

# Step 4.1 clustering
etom_drreport_clust = clusttrend(etom_drreport_mono, reference_index = NULL, sort.method =c("clust","layer"), sort.thres = 20, dist.method = "euclidean", hclust.method = "average",
                     mztag = "mzmed", rttag = "rtmed", heatmap.on = TRUE, plot.all = TRUE, filename = "testdataset_hclust_20clusters.pdf")

# Step 5.1 PCA of all features with color-coded ED50 values
plotpca(DoseResponse_report = etom_drreport_fit, DoseStat = etom_dosestat, EDrange = c(0,200))

# Step 2.2 inflection trend filtering
etom_drreport_reverse <- trendfilter(etom_dosestat, pval_cutoff = 0.05, pval_thres = 1, anova_cutoff = 0.05, trend = "reverse",
                                     relChange_cutoff = 0.05, export = TRUE)

# Step 3.2 inflection trend plotting
plottrend(etom_drreport_reverse, Dose_conditions = c("0uM","10uM","50uM","200uM"), y_transform = T,
          mz_tag = "mzmed", rt_tag = "rtmed")
