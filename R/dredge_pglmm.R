#' Model selection for pglmm
#'
#' This is an adaptation of the \code{dredge} function from the \code{MuMIn} package to
#' phylogenetic generalized linear mixed models (\code{pglmm} from the \code{phyr}). It
#' generates a model selection table of models as would do the \code{dredge}
#' function from \code{MuMIn} with combinations (subsets) of fixed effect terms in
#' the global model, with optional model inclusion rules.
#'
#' @param formulaRE A character string specifying the random effects structure (e.g., "~ 1 + (1 | Species)"). Must include an intercept (`1`) as fixed effect.
#' @param fixed A character vector of fixed effects to test. Interaction terms are allowed (e.g., "Sepal.Width:Petal.Length"), but some post-processing of tested models might be required.
#' @param data A data.frame containing the variables named in formula.
#' @param rank Criterion used to rank the models. Choose either "AIC" (default) or "AICc".
#' @param family Distribution family to use in model fitting. Options are "gaussian", "binomial", or "poisson".
#' @param cov_ranef A named list of covariance matrices of random terms. The names should be the group variables that are used as random terms with specified covariance matrices (without the two underscores, e.g. list(sp = tree1, site = tree2)). The actual object can be either a phylogeny with class "phylo" or a prepared covariance matrix. If it is a phylogeny, pglmm will prune it and then convert it to a covariance matrix assuming Brownian motion evolution. pglmm will also standardize all covariance matrices to have determinant of one. Group variables will be converted to factors and all covariance matrices will be rearranged so that rows and columns are in the same order as the levels of their corresponding group variables.
#' @param estimate Logical. If `TRUE` (default), includes fixed effect estimates in the output.
#' @param std.err Logical. If `TRUE`, includes standard errors for the fixed effects. Default is `FALSE`.
#' @param round Integer. Number of decimal places to round numeric results to. Default is `3`.
#'
#' @returns A data frame summarizing model selection results, ranked by the chosen criterion.
#' @export
#' @importFrom stats as.formula
#' @importFrom utils combn
#' @examples
#'
#' mod1 <- phyr::pglmm(log_onset ~ 1 + (1|Species__),
#'                   data = senescence, family = "gaussian", cov_ranef = list(Species = tree_ultra),
#'                   REML = TRUE, verbose = FALSE, s2.init = .1)
#'
#' dredge_pglmm(
#'        formulaRE = "log_onset ~ 1 + (1 | Species)",
#'        fixed = c("log_afr", "log_mass"),
#'        data = senescence,
#'        rank = "AICc",
#'        family = "gaussian",
#'        cov_ranef = list(Species = tree_ultra)
#'  )
dredge_pglmm <- function(formulaRE, fixed, data,
                         rank=c("AICc", "AIC"),
                         family=c("gaussian", "binomial", "poisson"),
                         cov_ranef,
                         estimate=T, std.err=F, round=3) {
  # global options
  options(na.print="")

  # NULL MODEL
  formula.RE <- as.formula(formulaRE) # converts the character random effect formula into a formula
  mod.RE <- phyr::pglmm(formula.RE, data, family,
                  cov_ranef = cov_ranef,
                  REML=FALSE, verbose=FALSE, s2.init=.1) # runs RE model

  df.RE <- nrow(lme4::fixef(mod.RE))+nrow(lme4::ranef(mod.RE)) # extracts df of null model
  logLik.RE <- mod.RE$logLik # extracts log-likelihood of the null model
  AICc.RE <- mod.RE$AIC + ((2*df.RE^2)+2*df.RE)/(nrow(data)-df.RE-1) # extracts AICc of the null model
  AIC.RE <- mod.RE$AIC # extracts AIC of the null model

  # FIXED-EFFECT MODELS
  df.fix <- c() # creates an empty vector to store df of non-null models
  logLik.fix <- c() # creates an empty vector to store log-likelihood of non-null models
  AICc.fix <- c() # creates an empty vector to store AICc of non-null models
  AIC.fix <- c() # creates an empty vector to store AIC of non-null models
  form.fix <- c() # creates an empty vector to store formula of non-null models
  model_list.fix <- c() # creats an empty vector to store the list of all non-null models tested

  comb.fix <- Map(combn, list(fixed), seq_along(fixed), simplify = FALSE) # creates all possible combination of fixed effects

  # runs all fixed-effect models
  for (i in 1:length(fixed)) {
    for (j in 1:length(comb.fix[[i]])) {
      ifelse(i==1, l <- comb.fix[[i]][[j]], l <- paste(comb.fix[[i]][[j]], collapse="+")) # puts all combination of fixed effects in the proper format
      form <- paste0("formula_", i,j) # creates variable names for each combination of fixed effects
      form.fix <- c(form.fix, assign(form, stringr::str_replace(formulaRE, "1", l) |> as.formula())) # includes fixed effects within the null formula provided in the input

      mod.fix <- phyr::pglmm(eval(parse(text = form.fix)), data, family,
                            cov_ranef = cov_ranef,
                            REML=FALSE, verbose=FALSE, s2.init=.1) # runs fixed-effect models

      df.fix <- c(df.fix, nrow(lme4::fixef(mod.fix))+nrow(lme4::ranef(mod.fix))) # extracts df of each model
      logLik.fix <- c(logLik.fix, mod.fix$logLik) # extracts log-likelihood of each model
      AICc.fix <- c(AICc.fix, mod.fix$AIC + ((2*(nrow(lme4::fixef(mod.fix))+nrow(lme4::ranef(mod.fix)))^2)+2*(nrow(lme4::fixef(mod.fix))+nrow(lme4::ranef(mod.fix))))/(nrow(data)-(nrow(lme4::fixef(mod.fix))+nrow(lme4::ranef(mod.fix)))-1)) # extracts AICc of the null model
      AIC.fix <- c(AIC.fix, mod.fix$AIC) # extracts AIC values for fixed-effects models

      model_list.fix <- c(model_list.fix, list(mod.fix))
    }
  }

  # DREDGE TABLE
  # Adds the null model to the list of models and to the df/AICc/logLik vectors
  #model_list <<- c(model_list.fix, list(mod.RE))
  model_list <- c(model_list.fix, list(mod.RE))
  df <- c(df.fix, df.RE)
  AICc <- c(AICc.fix, AICc.RE)
  AIC <- c(AIC.fix, AIC.RE)
  logLik <- c(logLik.fix, logLik.RE)

  all_fixed_effects <- unique(unlist(lapply(model_list, function(mod) rownames(lme4::fixef(mod)))))   # gets the unique fixed effects from all models

  # Creates an empty dataframe with columns named after the unique fixed effects
  dredge.table <- data.frame(matrix(ncol = length(all_fixed_effects), nrow = 0))
  colnames(dredge.table) <- all_fixed_effects

  # Iterates through each model and add a row to the dataframe with fixed effects values
  for (i in seq_along(model_list)) {
    model_fixed_effects <- rownames(lme4::fixef(model_list[[i]]))
    model_values <- as.data.frame(t(lme4::fixef(model_list[[i]])[,1]))
    model_se <- as.data.frame(t(lme4::fixef(model_list[[i]])[,2]))

    # Creates a template row with empty values
    template_row <- data.frame(matrix(NA, ncol = length(all_fixed_effects), nrow = 1))
    colnames(template_row) <- all_fixed_effects

    if (estimate==TRUE & std.err==TRUE){
      template_row[, model_fixed_effects] <- paste(round(model_values,round), paste0("(", round(model_se,round), ")"), sep = " ") # fills in the values and se for fixed effects present in the model
    }else{
      if (estimate==TRUE & std.err==FALSE){
        template_row[, model_fixed_effects] <- round(model_values,round) # fills in the values fixed effects present in the model
      }else{
        template_row[, model_fixed_effects] <- "+"
      }
    }
    dredge.table <- rbind(dredge.table, template_row) # adds the row to the dataframe
  }

  if (rank=="AICc"){
    dredge.table <- cbind(dredge.table, df, logLik, AICc) # adds corresponding dfs, logLiks and AICcs to the dataframe

    dredge.table <- dredge.table[order(dredge.table$AICc),] # orders the dataframe according to AICc value
    dredge.table$delta_AICc <- dredge.table$AICc - min(dredge.table$AICc) # calculates the delta_AICc column
    dredge.table$weight <- exp(-0.5*dredge.table$delta_AICc)/sum(exp(-0.5*dredge.table$delta_AICc)) # calculates the AICc weight columns
  }else{
    dredge.table <- cbind(dredge.table, df, logLik, AIC) # adds corresponding dfs, logLiks and AICcs to the dataframe

    dredge.table <- dredge.table[order(dredge.table$AIC),] # orders the dataframe according to AIC value
    dredge.table$delta_AIC <- dredge.table$AIC - min(dredge.table$AIC) # calculates the delta_AIC column
    dredge.table$weight <- exp(-0.5*dredge.table$delta_AIC)/sum(exp(-0.5*dredge.table$delta_AIC)) # calculates the AIC weight columns
  }

  # dredge.table <<- dredge.table

  # DISPLAYING
  full_model <- as.character(as.expression(model_list.fix[[length(model_list.fix)]]$formula_original))
  matches <- gregexpr("\\(1 \\|[^)]+\\)", formulaRE)
  selected_parts <- regmatches(formulaRE, matches)[[1]]
  random_part <- paste(selected_parts, collapse=" + ")

  print(paste("Global model: pglmm(", full_model,")", sep=""))
  cat("Data:\n")
  print(substitute(data))

  cat("---\nModel selection table\n")
  print(dredge.table)
  if (rank=="AICc"){
    cat("Models ranked by AICc\nRandom terms (all models):\n", random_part, sep="     ")
  }else{
    cat("Models ranked by AIC\nRandom terms (all models):\n", random_part, sep="     ")
  }
  cat("\nCovariance matrix of random effects:\n")
  print(substitute(cov_ranef))
}
