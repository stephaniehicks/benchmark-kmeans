#' Performs mini batch kmeans clustering on an hdf5 file
#' 
#' @param file name of hdf5 file.
#' @param name name of object in hdf5 file.
#' @param mbk MiniBatchKMeans object. Use to specify n_clusters, etc.
#' @return An array of integers indicating the cluster of each point
#' @examples
#' testmbk <- h5mbk("py-mouse-data.h5", "mousedata", SklearnEls()$skcl$MiniBatchKMeans(n_clusters = 3L))
h5mbk <- function(file = "", name = "", mbk = SklearnEls()$skcl$MiniBatchKMeans(), prediction = array(), first = TRUE, i = 1){
  f <- SklearnEls()$h5py$File(file, "r")
  data <- f$get(name)
  if(first){
    prediction <- array(dim = data$shape[2])
  }
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
  else {
    mbk$partial_fit(t(data$value[,i:(i+2)]))
    prediction[i:(i+2)] <- mbk$predict(t(data$value[,i:(i+2)]))
    f$close()
    return(h5mbk(file, name, mbk, prediction, first = FALSE,(i+3)))
  }
}