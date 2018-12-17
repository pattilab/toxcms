plottrend<- function(DoseResponse_report, Dose_conditions=NULL, y_tranform=TRUE, mz_tag = "mzmed", rt_tag = "rtmed",file="Reverse Reponse Curve.pdf"){
  Dose_Replicates <- DoseResponse_report$Dose_Replicates
  if(is.null(Dose_conditions)){
    cat("Dose conditions are not specified. Previous dose names are applied.\n")
    Dose_conditions=names(Dose_Replicates)
  } else if(length(Dose_conditions)!=length(Dose_Replicates)){
    cat("Dose conditions",Dose_conditions,"do not match",names(Dose_Replicates),".\n")
    stop("Abort")
  }
  # print Dose conditions and corresponding Dose names
  cat("Dose conditions",Dose_conditions,"are assigned to",names(Dose_Replicates),".\n")
  doses <- rep(Dose_conditions,times = Dose_Replicates)
  # assign doses as factor and assign dose levels for plotting
  doses <- factor(doses,levels=as.factor(Dose_conditions))
  num_feature <- nrow(DoseResponse_report$Feature)
  Feature <- DoseResponse_report$Feature
  
  mz <- Feature %>% dplyr::select(contains(mz_tag)) # get the m/z  from Feature
  rt <- Feature %>% dplyr::select(contains(rt_tag)) # get the rt from Feature
  mz <- round(mz,4)
  rt <- round(rt)
  # checking mz_tag and rt_tag
  if(length(mz)==1){
    cat("Succefully found",colnames(mz)," in feature table.\n")
  } else if(length(mz)==0) {
    stop("The specified '",mz_tag, "' does not match any columns in feature table.")
  } else if (length(mz)>1){
    stop("Multiple columns including '",mz_tag, "'. Please check the column names of feature table." )
  } else{
    stop("Error occurs when matching ", mz_tag)
  }
  if(length(rt)==1){
    cat("Succefully found",colnames(rt),"in feature table.\n")
  } else if(length(rt)==0){
    stop("The specified '",rt_tag, "'does not match any columns in feature table.")
  } else if (length(rt)>1){
    stop("Multiple columns including '",rt_tag, "'. Please check the column names of feature table." )
  } else{
    stop("Error occurs when matching ", rt_tag)
  }
  new_index <- sort(mz[[1]],index.return = TRUE)$ix
  last_index <- tail(new_index,1)
  # creat pdf for printing
  pdf(file=file, width = 8.5, height = 11)
  temp_grid <- list()
  j=1
  # looping to plot
  for (i in new_index){
    Normalized_Response <- as.numeric(DoseResponse_report$Normalized_Response[i,-1L])
    data <- data.frame(Doses=doses, Normalized_Response=Normalized_Response)
    a <- ggplot(data,aes(x=Doses,y=Normalized_Response,group=1)) + geom_point(size=5,alpha=I(0.8),color="red",shape=19) + geom_line(stat = 'summary',fun.y=mean,size=1.2)
    #y tranform
    if (y_tranform){
      a <- a + ggplot2::scale_y_sqrt()
    }
    # Set parameters for background
    a <- a + ggplot2::theme(strip.text.x = ggplot2::element_text(size = 16))
    a <- a + ggplot2::theme(panel.grid.minor = ggplot2::element_blank())
    a <- a + ggplot2::theme(panel.grid.major = ggplot2::element_blank())
    a <- a + ggplot2::theme_bw()
    a <- a + ggplot2::theme(axis.title.x = ggplot2::element_text(size = 10, margin = ggplot2::margin(15, 0, 0, 0)))
    a <- a + ggplot2::theme(axis.title.y = ggplot2::element_text(size = 10, margin = ggplot2::margin(0, 15, 0, 0)))
    a <- a + ggplot2::theme(axis.text.x = ggplot2::element_text(size = 10))
    a <- a + ggplot2::theme(axis.text.y = ggplot2::element_text(size = 10))
    # add labeles
    a <- a + ggplot2::labs(title=paste("index:", i, "  mz:", mz[i], "  rt:", rt[i]),x="Dose",y="Normalized Response") + ggplot2::theme(plot.title=element_text(size=10, hjust = 0.5, face="bold", color="black", lineheight=1))
    # print plot
    temp_grid[[j]] <- a
    j=j+1
    if(j==7 || i==last_index){
      do.call(grid.arrange,c(temp_grid,ncol=2,nrow=3,newpage=TRUE))
      temp_grid <- list()
      j=1
    }
  }
  dev.off()
  cat("Plotting finished. a pdf file is generated under:\n", getwd())
}