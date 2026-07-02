############################################################
## io.R
############################################################

GenerateSeuratObject <- function(sample){
  
  datadir <- file.path("../data", sample)
  
  message("Reading ", datadir)
  
  mydata <- Read10X(data.dir = datadir)
  
  CreateSeuratObject(
    mydata,
    project = sample
  )
  
}

ReadSeuratObjects <- function(samples, directory = "RDSfiles") {
  
  objs <- lapply(samples, function(sample) {
    
    readRDS(
      file.path(directory, paste0(sample, ".rds"))
    )
    
  })
  
  names(objs) <- samples
  
  return(objs)
}

CreateSampleDirectories <- function(samples, base.dir = "results") {
  
  # Create the base directory if needed
  if (!dir.exists(base.dir)) {
    dir.create(base.dir, recursive = TRUE)
  }
  
  # Create one directory per sample
  sample.dirs <- file.path(base.dir, samples)
  
  invisible(
    lapply(sample.dirs, function(dir) {
      if (!dir.exists(dir)) {
        dir.create(dir, recursive = TRUE)
      }
    })
  )
  
  sample.dirs <- setNames(
    file.path(base.dir, samples),
    samples
  )
  
  return(sample.dirs)
}