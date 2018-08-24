#' Performs mini batch kmeans clustering on an hdf5 file
#' @author Shenzhi Shi
#' @import BiocSklearn
#' @param file name of hdf5 file.
#' @param name name of object in hdf5 file.
#' @param mbk MiniBatchKMeans object. Use to specify n_clusters, etc.
#'
#' @return An array of integers indicating the cluster of each point
#' @examples
#' mbk <- SklearnEls()$skcl$MiniBatchKMeans(n_clusters = 3L)
#' testmbk <- h5mbk(file = "py-mouse-data.h5", name = "mousedata", mbk = mbk)
h5mbk <- function(file = "", name = "", mbk = SklearnEls()$skcl$MiniBatchKMeans(), prediction = array(), first = TRUE, i = 1){
  #get data
  f <- SklearnEls()$h5py$File(file, "r")
  data <- f$get(name)
  #if it's the first time running, set return value to the same length as the number of samples
  if(first){
    prediction <- array(dim = data$shape[2])
  }
  #check for the end of the samples
  if ((i+2) > data$shape[2]){
    if(i <= data$shape[2]){
      mbk$partial_fit(t(data$value[,i:data$shape[2]]))
      prediction[i:(i+2)] <- mbk$predict(t(data$value[,i:(i+2)]))
      f$close()
      return(prediction)
    }
    else{
      f$close()
      return(prediction)
    }
  }
  #read in a slice of 3 samples and update prediction with minibatchkmeans partial_fit
  else {
    mbk$partial_fit(t(data$value[,i:(i+2)]))
    prediction[i:(i+2)] <- mbk$predict(t(data$value[,i:(i+2)]))
    f$close()
    #continue loop
    return(h5mbk(file, name, mbk, prediction, first = FALSE,(i+3)))
  }
}
