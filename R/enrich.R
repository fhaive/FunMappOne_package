#library(gProfilerR)

#' @keywords internal
enrich = function(x, type, org, pval, adjust_method,sig = TRUE, mis = 0, only_annotated = TRUE){
  if(only_annotated){
    domain_size = "annotated"
  }else{
    domain_size = "known"
  }

  #print("before gprofiler")
  #out = gProfileR::gprofiler(query = as.character(x[,1]),src_filter=type,organism=org,domain_size = domain_size,
  #                           max_p_value = pval, correction_method = adjust_method,significant = sig,min_isect_size = mis)

  ##out = gProfileR::gprofiler(query = x[,1],src_filter=type,organism=org,max_p_value = pval,domain_size = "known")
  ##print(out)

  ###print("after gprofiler")
  #out=out[,c("term.id","intersection","p.value","p.value","term.name")]
  #colnames(out) = c( "annID","gID","pValue","pValueAdj","Description")
  out = gprofiler2::gost(query = as.character(x[,1]),sources = type,organism = org,domain_scope = domain_size,
                        user_threshold = pval,correction_method = adjust_method,significant = sig, evcodes = TRUE)
  out = out$result

  if(is.null(out)) out = data.frame(query= character(), significant= integer(), p_value= integer(), term_size = integer(),            
                                  query_size =  integer(), intersection_size =  integer(),precision =  integer(),recall =  integer(),               
                                  term_id = character(),source = character(),term_name = character(),effective_domain_size = integer(),
                                  source_order = integer(),parents = character(),evidence_codes = character(),intersection = character())

  ### UPDATE gprofiler version
  # out = gProfileR::gprofiler(query = as.character(x[,1]),src_filter=type,organism=org,domain_size = domain_size,
  #                            max_p_value = pval, correction_method = adjust_method,significant = sig,min_isect_size = mis)

  # out=out[,c("term.id","intersection","p.value","p.value","term.name")]
  out=out[,c("term_id","intersection","p_value","p_value","term_name")]
  colnames(out) = c( "annID","gID","pValue","pValueAdj","Description")

  #print(out)
  # if(nrow(out)>0){
  #   out$annID= gsub("KEGG:","",out$annID)
  #   if(org=="hsapiens")
  #     out$annID= gsub("REAC:","R-HSA-",out$annID)
  #   else{
  #     out$annID= gsub("REAC:","R-MMU-",out$annID)
  #   }

  if(nrow(out)>0)
    out$annID= gsub("KEGG:|REAC:","",out$annID)


  return(out)
}
