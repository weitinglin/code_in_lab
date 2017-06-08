#library in the source R code
# manage the time series 
source("/Users/Weitinglin/Documents/R_scripts/Lab/microarray/analysis/microarry_function.R")


#=====================input the data================================
data.path  <- file.path("/Users/Weitinglin/Documents/2016 實驗室資料處理/201510 microarray/raw data 20100114/set 4 time series")
experiment.set <- c("set4_p14" , "set4_p3", rep("set4_p6" ,3), "set4_p8")


celfile.set <- do_phenodata(data.path = data.path, experiment.set = experiment.set)

#=====================preprocess the data ======
#use the mas5,in order not to log twice
Exprs.data <- mas5 (celfile.set , normalize = FALSE, analysis = "absolute", sc = 500)

##expression profile
ecelfile.set <- exprs ( Exprs.data )

#========================quantile normalization===========================
#Method Quantile Normalization+log2

norm <- quantile_normalization(ecelfile.set = Exprs.data, method = "median")

norm.time.series <- norm
save(norm.time.series, file = "norm_qn_log_time.Rdata")



#



