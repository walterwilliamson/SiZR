#' Wrangle data for fitting in SiZR model
#'

#'
#' @param data
#'  A data frame containing the outcome variable and covariates.
#' @param shape
#'  A character string specifying the shape of the data. Options are "long" or "wide".
#' @param predictor
#'  A character string specifying the type of predictor. Options are "functional" or "non-functional".
#' @param Lmat
#'  An optional matrix of weights for the functional predictor. If NULL, equal weights will be used.
#'
#' @return A data frame formatted for SiZR model fitting.
#' @export
#'

make_data <- function(data,
                      shape = c("long", "wide"),
                      predictor = c("functional", "non-functional"),
                      Lmat = NULL){

  shape <- match.arg(shape)
  predictor <- match.arg(predictor)

  if(anyNA(data)){
    warning("Input data contains missing values.")
  }


  if (shape == "wide"){
    return(data)
  }

  if (shape == "long") {
    # Assumes 'Y' is the response and remaining columns are predictors
    if (predictor == "non-functional") {
      # In non-functional case, no need to expand anything
      return(data)

    } else if (predictor == "functional") {
      # 'X' should be a list-column of length-J numeric vectors
      if (!"X" %in% names(data)) {
        stop("Data must have a column named 'X' containing functional predictors.")
      }

      J <- length(data$X[[1]])
      n <- nrow(data)

      if (is.null(Lmat)) {
        Lmat <- matrix(1/J, nrow = n, ncol = J)
      } else {
        # Validate user-provided Lmat
        if (!is.matrix(Lmat) || nrow(Lmat) != n || ncol(Lmat) != J) {
          stop("Lmat must be a matrix with dimensions n x J (rows = observations, cols = time points).")
        }
      }

      colnames(Lmat) <- paste0("L", 1:J)
      Lmat_df <- as.data.frame(Lmat)

      # Convert list-column X into a wide format data frame
      Xmat <- do.call(rbind, data$X)
      colnames(Xmat) <- paste0("X", 1:J)
      X_df <- as.data.frame(Xmat)

      new_data <- cbind(data["Y"], X_df, Lmat_df)
      return(new_data)
    }
  }
}
