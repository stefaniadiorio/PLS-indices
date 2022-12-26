library(randcorr)
library(MASS)
#library(OrdNor)

pois.func <- function(y, lambda) { return(lambda^y * exp(-lambda) / factorial(y)) }

exploreOutcomes <- function(outcomes, probs, name, ...) {
  maxvalue <- if (is.character(name) && name == "poisson") 50 else 1000
  minvalue <- if (is.character(name) == "poisson") -50 else -1000
  
  outcomes[outcomes == Inf] <- maxvalue
  outcomes[outcomes == -Inf] <- minvalue
  
  myprobs <- probs(outcomes[1]:outcomes[2], ...)
  myouts <- (outcomes[1]:outcomes[2])[!is.nan(myprobs)]
  myprobs <- myprobs[!is.nan(myprobs)]
  
  return(myouts[myprobs > .Machine$double.eps^0.5])
}

RV <- function(outcomes, probs = NULL, odds = NULL, fractions = (class(probs) != "function"), range = any(is.infinite(outcomes)), verifyprobs = TRUE, id = rnorm(1), ...) {
  outtry <- outcomes
  if (length(outcomes) == 1 && outcomes == "bernoulli") {
    probs <- function(x, ...) { p <- ifelse(length(list(...)) == 0, .5, list(...)[[1]]); p^x * (1 - p)^(1 - x) }
    outtry <- 0:1
  } else if (length(outcomes) == 1 && outcomes == "poisson") {
    probs <- function(x, ...) { lambda <- ifelse(length(list(...)) == 0, 5, list(...)[[1]]); lambda^x * exp(-lambda) / factorial(x) }
    range <- TRUE
    fractions <- FALSE
    outtry <- c(0, Inf)
  } else if (length(outcomes) == 1 && outcomes == "geometric") {
    probs <- function(x, ...) { p <- ifelse(length(list(...)) == 0, .5, list(...)[[1]]); (1 - p)^(x - 1) * p }
    range <- TRUE
    fractions <- FALSE
    outtry <- c(1, Inf)
  }
  
  test <- fractions # TODO: Fix
  old.out <- outtry
  
  if (range) outcomes <- suppressWarnings(exploreOutcomes(outtry, probs, outcomes, ...))     else outcomes <- outtry
  
  if (class(probs) == "function") probs <- suppressWarnings(probs(outcomes, ...))
  
  pr <- probs
  if (is.null(pr)) pr <- odds
  
  probsSum <- sum(pr)
  
  if ((probsSum > (1 + .Machine$double.eps^0.5)) && is.null(odds) && verifyprobs) stop("Probabilities sum to over 1")
  if (any(pr < 0)) stop("Probabilities cannot be negative")
  
  isOdds <- !is.null(odds)
  
  if (length(outcomes) < length(pr)) {
    stop("More probabilities/odds than outcomes provided")
  } else if (length(outcomes) > length(pr)) {
    pr <- c(pr, rep(ifelse(isOdds, 1, (1 - probsSum) / (length(outcomes) - length(pr))), length(outcomes) - length(pr)))
  }
  
  ## Convert to probs
  if (verifyprobs) probs <- pr / sum(pr)
  names(probs) <- outcomes
  
  ## Remove zero prob events
  ind <- (probs > .Machine$double.eps^0.5)
  outcomes <- outcomes[ind]
  probs <- probs[ind]
  
  if (any(duplicated(outcomes))) {
    out <- NULL
    my.df <- data.frame(out = as.vector(outcomes), pr = as.numeric(probs))
    my.sum <- ddply(my.df, .(out), summarise, pr = sum(pr))
    outcomes <- type.convert(as.character(my.sum$out))
    probs <- as.numeric(my.sum$pr)
  }
  return(list(categories = outcomes, probabilities = probs))
}

create_plist <- function(p,max_num_categories, min_num_categories=2){
  plist <- vector(mode="list",p)
  categories = max_num_categories
  for (i in 1:p){
    p_outcomes <- round(runif(1,min_num_categories,categories),0)
    print(c(i,"-th variable has ",p_outcomes," categories"))
    lmda <- round(runif(1,1,10),0)
    p_dist = RV(outcomes = 1:p_outcomes, probs = pois.func,lambda=lmda)
    p_probs = p_dist$probabilities
    print(p_probs)
    p_cdf = c(cumsum(p_probs)[-p_outcomes])
    print(p_cdf)
    plist[[i]] <- p_cdf
  }
  return(plist)
}
