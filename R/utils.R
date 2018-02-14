# copied from plyr to avoid the dependency on plyr
compact <- function (l)
  Filter(Negate(is.null), l)

LdFlags <- function(){
 
  lib <- "/libflowWorkspace.a"
    
  libpaths <- paste0("lib", Sys.getenv("R_ARCH"), lib)
  libpaths <- lapply(libpaths, function(libpath)tools::file_path_as_absolute(base::system.file(libpath, package = "flowWorkspace" )))
  cat(paste(libpaths, collapse = " "))
}

#Taken from limma (don't want to import and create a dependency)
trimWhiteSpace<-function (x)
{
  sub("[ \t\n\r]*$", "", sub("^[ \t\n\r]*", "", x))
}

#' make a formula from a character vector
#'
#' construct a valid formula to be used by flowViz::xyplot
#'
#' @param dims a \code{character} vector that contains y , x axis, if it is unnamed, then treated as the order of c(y,x)
#' @param isChar \code{logical} flag indicating whehter to return a formula or a pasted string
#' @return when \code{isChar} is TRUE, return a character, otherwise coerce it as a \code{formula}
#' @examples
#' all.equal(mkformula(c("SSC-A", "FSC-A")),`SSC-A` ~ `FSC-A`)#unamed vecotr
#' all.equal(mkformula(c(x = "SSC-A", y = "FSC-A")),`FSC-A` ~ `SSC-A`)#named vector
#' @export
mkformula<-function(dims,isChar=FALSE){
  if(length(dims)==1){
    form<-paste(c("",sapply((dims), function(x) paste("`",x, "`", sep = ""))), collapse = "~")
  }else{
    
    dnames <- names(dims)
    if(!is.null(dnames)){
      if(isTRUE(all.equal(sort(dnames), c("x", "y"))))
        dims <-  dims[rev(order(names(dims)))]
      else
        warning("invalid axis names: ", paste(dnames, collapse = ","), "(expect 'x' or 'y')")
    }
    form <- paste(sapply((dims),function(x)paste("`",x,"`",sep="")),collapse="~")
  }
  if(!isChar)
    form<-as.formula(form)
  return(form)
}
#' Get Cell Population Statistics and Sample Metadata
#'
#' @param object a \code{GatingSet} or \code{GatingSetList}
#' @param ... additional arguments passed to \code{getPopStats}
#'
#' @return a \code{data.table} of merged population statistics with sample metadata.
#' @export
#' @importFrom dplyr inner_join
#' @examples
#'  \dontrun{
#'     #G is a GatingSetList
#'     stats = getMergedStats(G)
#'   }
getMergedStats = function(object,...){
	if(!inherits(object,"GatingSet")&!inherits(object,"GatingSetList")){
		stop("object must be a GatingSet or GatingSetList")
	}
	stats = getPopStats(object,...)
	#process name column so that it contains XXX.fcs
	message("Processing sample names..")
	stats[, sampleName:=name]
	stats[,name := gsub("(fcs).*","\\1",name)]
	pd = pData(object)
	message("merging..")
	ret = inner_join(stats,pd,by="name")
	message("Done")
	return(ret)
}

#' save the event counts parsed from xml into c++ tree structure
#'
#' It is for internal use by the diva parser
#'
#' @param gh GatingHierarchy
#' @param node the unique gating path that uniquely identifies a population node
#' @param count integer number that is events count for the respective gating node directly parsed from xml file
#' @export
#' @examples
#' \dontrun{
#' set.count.xml(gh, "CD3", 10000)
#' }
set.count.xml <- function(gh, node, count){
  .set.count.xml(gh@pointer, sampleNames(gh), node, count)
}