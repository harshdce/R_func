library(lavaan)
CDnV<- function(model,data)
{
  extract_factors <- function(model_string) {
    # Use regular expressions to find factor names
    factors <- gregexpr("([a-zA-Z0-9_]+)\\s*=~", model_string, perl = TRUE)
    factor_names <- regmatches(model_string, factors)
    
    # Clean the extracted names
    factor_names <- unlist(factor_names)
    factor_names <- gsub("\\s*=~", "", factor_names)
    return(factor_names)
  }
  
  AVE <- function(std_solution, latent_variable)   {
    loadings <- std_solution[std_solution$lhs == latent_variable & std_solution$op == "=~", "est.std"]
    AVE_value <- sum(loadings^2) / length(loadings)
    return(AVE_value)
  }
  CFA<-cfa(model,data)
  std_solution <- standardizedSolution(CFA)
  n_fac=length(extract_factors(model))
  factors_name=extract_factors(model)
  val_matrix <- matrix(NA, nrow = n_fac, ncol = n_fac) 
  rownames(val_matrix) <- colnames(val_matrix) <- factors_name
  
  latent_cors <- lavInspect(CFA, "cor.lv")
  # Assign off-diagonal value (squared correlation)
  for (i in 1:n_fac) {
    fn <- factors_name[i]
    AVE_fn <- AVE(std_solution, fn)
    val_matrix[factors_name[i], factors_name[i]] <- sqrt(AVE_fn)
    
    for (j in 1:n_fac) {
      if (j< i) {
        fm <- factors_name[j]
        val_matrix[factors_name[i], factors_name[j]] <- latent_cors[fn, fm]
      }
    }
  }
  return(val_matrix)
}
