#' Parse a flowJo Workspace
#' 
#' Deprecated by CytoML::flowJo_to_GatingSet function.
#' @param ...
#' @export 
parseWorkspace <-function(...){
     .Deprecated("CytoML::flowjo_to_gatingset")
    }
    
#' Open/Close a flowJo workspace
#' @param ... other arguments passed to \code{\link{xmlTreeParse}}
#' @export 
openWorkspace <- function(...)
{
  .Deprecated("CytoML::open_flowjo_workspace")
  
}