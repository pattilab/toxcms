plotDr4pl <- function(x,
                      type.curve = "all",
                      text.title = "Dose-response plot",
                      text.x = "Dose",
                      text.y = "Normalized Response",
                      indices.outlier = NULL,
                      breaks.x = NULL,
                      breaks.y = NULL,
                      dose_transform = TRUE,
                      ...) {
  
  ### Check whether function arguments are appropriate
  if(!is.character(text.title)) {
    stop("Title text should be characters.")
  }
  if(!is.character(text.x)) {
    stop("The x-axis label text should be characters.")
  }
  if(!is.character(text.y)) {
    stop("The y-axis label text should be characters.")
  }
  ### Draw a plot
  n <- x$sample.size
  color.vec <- rep("blue", n)
  shape.vec <- rep(19, n)
  if(!is.null(indices.outlier)) {
    color.vec[indices.outlier] <- "red"
    shape.vec[indices.outlier] <- 19
  }
  a <- ggplot2::ggplot(aes(x = x$data$Dose, y = x$data$Response), data = x$data)
  if(type.curve == "all") {
    a <- a + ggplot2::stat_function(fun = MeanResponse,
                                    args = list(theta = x$parameters),
                                    size = 1.2)
  }
  a <- a + ggplot2::geom_point(size = 5, alpha = I(0.8), color = color.vec,
                               shape = shape.vec)
  a <- a + ggplot2::labs(title = text.title,
                         x = text.x,
                         y = text.y)
  # Set parameters for the grids
  a <- a + ggplot2::theme(strip.text.x = ggplot2::element_text(size = 16))
  a <- a + ggplot2::theme(panel.grid.minor = ggplot2::element_blank())
  a <- a + ggplot2::theme(panel.grid.major = ggplot2::element_blank())
  if(!is.null(breaks.x)) { 
    a <- a + ggplot2::scale_x_sqrt(breaks = breaks.x)
  } else if (dose_transform){ 
    a <- a + ggplot2::scale_x_sqrt()
  } else {}
  if(!is.null(breaks.y)) {
    a <- a + ggplot2::scale_y_continuous(breaks = breaks.y)
  } else { 
    a <- a + ggplot2:: scale_y_continuous()
  }
  a <- a + ggplot2::theme_bw()
  # Set parameters for the titles and text / margin(top, right, bottom, left)
  a <- a + ggplot2::theme(axis.title.x = ggplot2::element_text(size = 10, margin = ggplot2::margin(15, 0, 0, 0)))
  a <- a + ggplot2::theme(axis.title.y = ggplot2::element_text(size = 10, margin = ggplot2::margin(0, 15, 0, 0)))
  a <- a + ggplot2::theme(axis.text.x = ggplot2::element_text(size = 10))
  a <- a + ggplot2::theme(axis.text.y = ggplot2::element_text(size = 10))
  return(a)
}
