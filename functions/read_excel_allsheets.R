read_excel_allsheets <- function(filename, tibble = FALSE) {
  sheets <- readxl::excel_sheets(filename)
  x <- lapply(sheets, function(X){ 
    
    y = readxl::read_excel(filename, sheet = X)
    y = as.data.frame(y)
    #y[,1] = as.character(as.vector(y[[1]]))
    rownames(y) = y[,1]
    y = y[,-1]
    y
  })
  
  if(!tibble) x <- lapply(x, as.data.frame)
  names(x) <- sheets
  return(x)
  
}

# DF = read_excel_allsheets(filename = "../pheno_list.xlsx",tibble = FALSE)
# Exp = read_excel_allsheets(filename = "../exp_mat_file.xlsx",tibble = FALSE)
