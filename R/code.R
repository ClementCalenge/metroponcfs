#' metroponcfs: A package for implementing the Metropolis algorithm
#'
#' The metroponcfs provides many interesting functions for
#' implementing the Metropolis algorithm, but the two main help pages
#' of interest are \code{GeneralSingleMetropolis}, which implements
#' and demonstrate how to use the Metropolis algorithm, and
#' \code{c.GSMetrop}, which demonstrate how to combine different
#' chains in an object.
#'
#'
#' @docType package
#' @name metroponcfs
NULL

#' @title Common Teal captured and Ringed in Abberton Reservoir, UK,
#' and Camargue, France.
#'
#' @details
#' Dataset describing the recapture locations of common teals
#' initially captured and ringed in two Western European places: (i)
#' in Abberton Reservoir, UK, and (ii) in Camargue, sourthern
#' France. See Guillemain et al. for a complete description of this
#' dataset.  To circumvent copyright issues, we provide here an
#' altered version of the dataset used by these authors: we selected a
#' random sample of 75\% of the recaptures of the original data,
#' keeping only the location information (i.e. only the x and y
#' coordinates of the recaptures), and we added a random noise to
#' these locations (we moved every bird recapture location randomly by
#' a distance comprised between 0 and 100 km).
#'
#'
#' @format A list with the following elements:
#'
#' \describe{
#'
#' \item{recaptures}{this element is a data.frame containing the
#' following information for each recapture of common teal: (i) a
#' variable named \code{date}, storing the number of years elapsed
#' between the initial capture and recapture, (ii) a variable named
#' \code{Abberton}, indicating whether the recaptured animal was
#' initially captured at Abberton Reservoir, UK (=1) or in Camargue,
#' southern France (=0), and (iii) two variables named \code{x} and
#' \code{y} containing the coordinates of the recapture (coordinate
#' system: Lambert azimuthal equal-area).}
#'
#' \item{rotationMatrix}{this 2x2 matrix M contains the two vectors
#' \code{cbind(m1, m2)} used by Guillemain et al. to rotate the
#' geographical coordinates in a new coordinate system (see Guillemain
#' et al.). Thus, if \code{x} is a vector of length two containing the
#' x and y coordinates of a location in the Lambert azimuthal
#' equal-area system, \code{z <- x\%*\%M} contains the coordinates of
#' this point in this new system.}
#'
#' \item{inverseRotationMatrix}{this 2x2 matrix \code{R} allows to
#' transform the coordinates of a point from the new coordinate system
#' to the old one. Thus, if \code{z} is a vector of length two
#' containing the coordinates of a point in the new coordinate system,
#' \code{x <- z\%*\%R} contains the coordinates of this point in the
#' original Lambert azimuthal equal-area system.}
#'
#' \item{knots}{this vector contains the coordinates of
#' the 26 knots in the new coordinate system used to define the
#' B-spline basis in the paper of Guillemain et al.}
#'
#' \item{lipum}{a list containing the parameters of the updating
#' mechanisms used in the Metropolis algorithm.}
#' }
#'
#' @references Guillemain M., Calenge C., Champagnon J. and Hearn
#' R. in prep. Determining the boundaries and plasticity of migratory
#' bird flyways: a Bayesian model for Common Teal Anas crecca
#' in Western Europe.
#'
"recteal"



#' @title Bayesian Model Fitted to Estimate the Boundaries of Common
#' Teal Flyways
#'
#' @details This model is fully described in the vignette "flyways".
#' Please type \code{vignette("flyways")} for more details.
#'
#' @format An object of class \code{CMM} (see
#' \code{help("c.GSMetrop")} for additional details on objects of
#' class \code{CMM}.
#'
#' @references Guillemain M., Calenge C., Champagnon J. and Hearn
#' R. in prep. Determining the boundaries and plasticity of migratory
#' bird flyways: a Bayesian model for Common Teal Anas crecca
#' in Western Europe.
"rectealmodel"


#' @title Bayesian Model Fitted to Estimate the Boundaries of Common
#' Teal Flyways (time effect)
#'
#' @details This model is fully described in the vignette "flyways".
#' Please type \code{vignette("flyways")} for more details.
#'
#' @format An object of class \code{CMM} (see
#' \code{help("c.GSMetrop")} for additional details on objects of
#' class \code{CMM}.
#'
#' @references Guillemain M., Calenge C., Champagnon J. and Hearn
#' R. in prep. Determining the boundaries and plasticity of migratory
#' bird flyways: a Bayesian model for Common Teal Anas crecca
#' in Western Europe.
"rectealmodeltime"




#' @title Find Starting Values for the Metropolis Algorithm
#' @export
#'
#' @details \code{findStartingValues} allows to find a list of starting
#' values for the Metropolis algorithm for the MCMC fitting of Bayesian models.
#'
#' @param support named list with one element per parameter (the names
#' of the list correspond to the name of the parameters).  Each
#' element should be a vector of length 2 indicating the limits of the
#' support of the parameters (i.e. the min and max values in which the
#' starting values will be searched).
#' @param logposterior a function for the calculation of the log
#' posterior for the current model. This function must rely only on
#' the data provided in \code{lidat}. This function must have three
#' parameters: (i) \code{par} is the list of parameters, (ii)
#' \code{lidat} is the list containing the data, and (iii) \code{ctrl}
#' is a named list (see the help page of
#' \code{GeneralSingleMetropolis}). Note that the argument \code{ctrl}
#' is not used by the function here.
#' @param lidat a named list containing the data (the names of the
#' list correspond to the name of the variables used in logposterior,
#' combined with the parameters), required for the calculation of the
#' posterior
#' @param multiple a named list describing which parameters are
#' vectors of parameters. For example, if a model relies on a single
#' parameter \code{"a"} and two vectors of parameters \code{"theta1"}
#' and \code{"theta2"} of length 5 and 6 respectively, the argument
#' \code{multiple} should take the value: \code{list(theta1=5,
#' theta2=6)}.  Note that single parameters are not indicated in this
#' list.
#' @param nrand the number of starting values to generate before a
#' result is returned (see details).
#' @param method the method to be used by the function (see details).
#' @param ndispersed integer. When \code{method="dispersed"},
#' indicates how many dispersed starting values should be
#' returned. Should be greater than \code{nrand}.
#' @param info whether information on the generation process should be
#' displayed to the user.
#'
#' @details Three methods are proposed to generate one or several
#' vectors of starting values. All of them rely on the generation of
#' \code{nrand} plausible starting values (i.e. leading to finite
#' log-posterior) randomly sampled in the support of the
#' parameters. Then: (i) the (default) method \code{"onebest"} chooses
#' the starting value leading to the highest log-posterior, (ii) the
#' method \code{"dispersed"} chooses \code{ndispersed} values the most
#' different among each other (to start several chains), (iii) the
#' method \code{"all"} return all sampled values.  The user can play
#' with the "support" of the parameters to find more or less dispersed
#' values.
#'
#' @return If \code{method=="onebest"}, a named list of starting
#' values for the parameters, that can be used directly as an argument
#' to \code{GeneralSingleMetropolis}.  For other methods, a list of
#' lists of starting values.
#' @seealso \code{\link{GeneralSingleMetropolis}}
#' @author Clement Calenge, \email{clement.calenge@@oncfs.gouv.fr}
#' @import graphics
#' @import stats
#' @examples
#'
#' ## We generate data
#' x1 <- rnorm(100)
#' x2 <- rnorm(100)
#' x3 <- rnorm(100)
#' y <- 1+x1+0.5*x2+1.5*x3+rnorm(100, sd=0.5)
#'
#' ## Our data
#' lid <- list(x=cbind(x1,x2,x3), y=y)
#'
#' ## 5 parameters: intercept + 3 slopes + residual variance
#' ## For the sake of demonstration, we put the three slopes in
#' ## a vector theta for the estimation.
#' ## The list of starting values should therefore contain: (i) theta,
#' ## (ii) intercept, (iii) residprec
#' ##
#' ## The function for the posterior is:
#' lpost <- function(par, lidat, ctrl)
#' {
#'      ## prior:
#'      parb <- par
#'      parb$residprec <- NULL
#'      lprior <- sum(dnorm(unlist(parb), mean=0, sd=100)) ## vague prior
#'      ## (uniform prior for the residual precision, we do not add this constant)
#'
#'      mu <- lidat$x%*%par$theta+par$intercept
#'      ## log-posterior
#'      return(sum(lprior + dnorm(lidat$y, mean=mu, sd = 1/sqrt(par$residprec), log=TRUE)))
#' }
#'
#' ##  We define large support for all parameters
#' support <- list(theta=c(-10,10), intercept=c(-10,10), residprec=c(0.0001, 10))
#' ## only theta is multiple (3 elements), so:
#' multiple <- list(theta=3)
#'
#' ## Default method: finds the best starting value among 1000 sampled values
#' ## The plot shows the distribution of the sampled posterior (up to a constant)
#' ## the best value is returned in sv
#' (sv <- findStartingValues(support, lpost, lid, multiple, nrand=1000))
#'
#' ## Same, but keeps 4 most dispersed starting values among 100 generated values:
#' (sv2 <- findStartingValues(support, lpost, lid, multiple, nrand=100,
#'                            method="dispersed", ndispersed=4))
#'
findStartingValues <- function(support, logposterior, lidat,
                              multiple=list(), nrand=10,
                              method=c("onebest","dispersed", "all"),
                              ndispersed=3, info=TRUE)
{
    if (length(names(support))<1)
        stop("support should be a named list, with one element per parameter")
    if (!all(sapply(support, length)==2))
        stop("Each element of support should be a vector of length 2\ncontaining the limits (min, max) of the support for the corresponding parameter")
    if (!all(names(multiple)%in%names(support)))
        stop("all the names of the parameter multiple should correspond to names in support")
    method <- match.arg(method)
    if (method=="dispersed") {
        if (ndispersed>nrand)
            stop("ndispersed should be <= nrand")
    }
    lpar <- lapply(1:nrand, function(r) {
                       if (info)
                           cat("Iteration", r, "out of", nrand,"\r")
                       ok <- FALSE
                       while (!ok) {
                           parinit <- lapply(1:length(support), function(i) {
                                                 n <- ifelse(names(support)[i]%in%names(multiple),
                                                             multiple[[names(support)[i]]],1)
                                                 return(runif(n, support[[i]][1], support[[i]][2]))
                                             })
                           names(parinit) <- names(support)
                           post <- logposterior(parinit, lidat, list(namePar="No221cc9328Parameter23",iter=1))
                           if (!is.infinite(post)) {
                               ok <- TRUE
                           }
                       }
                       attr(parinit,"lposterior") <- logposterior(parinit,lidat, list(namePar="No221cc9328Parameter23",iter=1))
                       return(parinit)
                   })
    sp <- sapply(lpar, function(x) attr(x, "lposterior"))
    if (info)
        plot(sort(sp, decreasing=TRUE), ylab = "log-posterior", xlab = "Trials")
    if (method=="all")
        return(lpar)
    lp <- lpar[[which.max(sp)]]
    if (method=="onebest")
        return(lp)
    if (method=="dispersed") {
        if (info)
            cat("Identifying dispersed values...\r")
        la <- do.call(rbind,lapply(lpar, function(x) unlist(x)))
        la2 <- scale(la)
        ## Solution to find proposed here:
        ## https://stackoverflow.com/questions/22152482/choose-n-most-distant-points-in-r
        di <- as.matrix(dist(la2))
        pa <- 1:nrow(la2)
        while (length(pa) > ndispersed) {
            cdists = rowSums(di)
            closest <- which(cdists == min(cdists))[1]
            pa <- pa[-closest]
            di <- di[-closest,-closest]
        }
        return(lpar[pa])
    }
}


.iterMetrop <- function(par, lidat, logposterior, lipum, NamesP, UpdatingMechanism,
                        liopt)
{
    ## Posterior courante (inutilisée si optimisation de la mise à jour d'un paramètre)
    postcur <- logposterior(par, lidat, list(namePar="No221cc9328Parameter23",iter=1))
    opti <- FALSE

    ## "working parameter"
    parw <- par

    ## liste d'acceptation
    accept <- list()

    ## Mise à jour de a:
    for (j in 1:length(parw)) {

        ## Les paramètres nécessaires
        nam <- NamesP[j]
        fun <- UpdatingMechanism[j]

        ## Si optimisation pour le paramètre précédent et pas pour le paramètre courant,
        ## recalculer la posterior courante
        if (opti&(!(nam%in%liopt))) {
            postcur <- logposterior(par, lidat, list(namePar="NoBPo221ccParameter23",iter=1))
        }

        ## Mise à jour de opti pour le paramètre courant
        opti <- nam%in%liopt

        ## Mise en œuvre de la fonction
        liparam <- list(par=parw, lidat=lidat, nam=nam, pum=lipum[[nam]],
                        ctrl=list(namePar=nam,iter=1),
                        logposterior=logposterior,
                        postcur=postcur, optimiser=opti)
        resm <- do.call(fun, liparam)

        ## Résultats
        parw[[nam]] <- resm$parm
        accept[[nam]] <- resm$accept
        postcur <- resm$postcur
    }

    return(list(par=parw, accept=accept, postcur=postcur))
}



#' @title Updating mechanisms used in the metropolis algorithm
#' @export
#' @details \code{singleStep} implements an updating mechanism for a
#' unique parameter based on a normal distribution.
#' \code{MultipleIndependentSteps} implements the same updating
#' mechanism for a vector of parameters. \code{MultiNormalStep}
#' implements an updating mechanism relying on the sampling of a
#' multinormal distribution. These functions are used by the package,
#' but are not to be used directly by the user.
#'
#' @param par a named list containing the parameters of the model (the
#' names of the list correspond to the names of the parameters)
#' @param lidat a named list containing the variables and constants
#' required by the model (especially the function \code{logposterior}
#' @param nam character string. The name of the parameter to be
#' updated.
#' @param pum for \code{singleStep}, the standard deviation of the
#' Gaussian distribution used to generate the proposal. For
#' \code{MultipleIndependentSteps}, a vector containing the standard
#' deviations for the Gaussian distribution used to generate the
#' proposal of every element of the vector. For \code{MultiNormalStep},
#' the covariance matrix of the multinormal distribution.
#' @param ctrl a named list with two parameters: (i) \code{namePar} is
#' the name of the parameter in \code{par} that is updated, (ii) if
#' this parameter is a vector where each component is updated in turn
#' (e.g. for \code{MultipleIndependentSteps}), \code{iter} indicates
#' the element of this vector that is updated (this element is ignored
#' in other cases).
#' @param logposterior a function used to calculate the log-posterior
#' probability of \code{par} (up to a constant).  This function must
#' have two arguments: \code{par} is the list of parameters and
#' \code{lidat} is the list containing the data used in the model.
#' @param postcur the current value of the posterior, before updating.
#' @param optimiser logical. If \code{TRUE}, the log-posterior
#' calculation is optimized for the current parameter (i.e. the
#' posterior is proportional to the log-posterior but the
#' proportionality constant is not necessarily the same for other
#' optimized parameters, so that the current posterior should be
#' calculated before every step).  If \code{FALSE} the function
#' \code{logposterior} returns the exact log-posterior.  In other
#' words, if \code{optimiser} is TRUE, the function \code{singleStep},
#' etc. calculates the current posterior (prior to updating) and does
#' not take into account the value of \code{postcur} passed to the
#' function. If FALSE, the function relies on the value of
#' \code{postcur} passed as argument for the updating.
#'
#' @return A list with three elements: (i) parm is the updated vector
#' of parameters, (ii) postcur is the value of the (optionally
#' optimized) posterior after updating, and (iii) accept is a vector
#' containing TRUE if an element proposed has been accepted and FALSE
#' otherwise
#'
#' @author Clement Calenge, \email{clement.calenge@@oncfs.gouv.fr}
#' @seealso \code{\link{GeneralSingleMetropolis}}
#' @importFrom MASS mvrnorm
#' @import graphics
#' @import stats
singleStep <- function(par, lidat, nam, pum, ctrl, logposterior,
                       postcur, optimiser)
{
    parw <- par
    prev <- parw[[nam]]
    ## si optimisation, on doit calculer la posterior à priori
    ## sinon, on prend l'argument
    if (optimiser) {
        postcur <- logposterior(parw, lidat, ctrl)
    }
    nex <- prev+rnorm(1, mean=0, sd=pum)
    parw[[nam]] <- nex
    postnext <- logposterior(parw, lidat, ctrl)
    if (is.infinite(postnext)&is.infinite(postcur)) {
        tf <- tempfile("dumpmetrop", tmpdir=".", fileext = ".Rdata")
        li <- list(par=par, lidat=lidat, nam=nam, ctrl=ctrl, logposterior=logposterior,
                   postcur=postcur, postnext=postnext, optimiser=optimiser, i=NA,
                   nex=nex, prev=prev, parw=parw, functionError="SingleStep")
        save(li, file=tf)
        msg <- paste0("Error with update of parameter ", nam,
                      ":\ninfinite current and next log-posterior.\n",
                      "All necessary information has been dumped in file", tf,"\n",
                      "Reload this file to explore the state of the chain when the error occured")
        stop(msg)
    }
    alpha <- min(c(1,exp(postnext-postcur)))
    un <- runif(1)<alpha
    accept <- as.numeric(un)
    parm <- ifelse(un, nex, prev)
    pc <- ifelse(un, postnext, postcur)
    return(list(parm=parm, postcur=postcur, accept=accept))
}


#' @rdname singleStep
#' @export
MultipleIndependentSteps <- function(par, lidat, nam, pum, ctrl, logposterior,
                                     postcur, optimiser)
{
    parw <- par
    parm <- parw[[nam]]
    accept <- rep(0,length(parm))

    for (i in 1:length(parm)) {
        ctrl$iter <- i
        if (optimiser) {
            postcur <- logposterior(parw, lidat, ctrl)
        }
        prev <- parm[i]
        nex <- prev+rnorm(1, mean=0, sd=pum[i])
        parm[i] <- nex
        parw[[nam]] <- parm
        postnext <- logposterior(parw, lidat, ctrl)
        if (is.infinite(postnext)&is.infinite(postcur)) {
            tf <- tempfile("dumpmetrop", tmpdir=".", fileext = ".Rdata")
            li <- list(par=par, lidat=lidat, nam=nam, ctrl=ctrl, logposterior=logposterior,
                       postcur=postcur, postnext=postnext, optimiser=optimiser, i=i,
                       nex=nex, prev=prev, parw=parw, functionError="MultipleIndependentSteps")
            save(li, file=tf)
            msg <- paste0("Error with update of parameter ", nam,
                          ":\ninfinite current and next log-posterior.\n",
                          "All necessary information has been dumped in file", tf)
            stop(msg)
        }
        alpha <- min(c(1,exp(postnext-postcur)))
        un <- runif(1)<alpha
        accept[i] <- un
        parm[i] <- ifelse(un, nex, prev)
        postcur <- ifelse(un, postnext, postcur)
    }

    return(list(parm=parm, postcur=postcur, accept=accept))
}


#' @rdname singleStep
#' @export
MultiNormalStep <- function(par, lidat, nam, pum, ctrl, logposterior,
                            postcur, optimiser)
{
    parw <- par
    prev <- parw[[nam]]
    if (optimiser) {
        postcur <- logposterior(parw, lidat, ctrl)
    }
    nex <- prev+MASS::mvrnorm(1, mu=rep(0,ncol(pum)), Sigma=pum)
    parw[[nam]] <- nex
    postnext <- logposterior(parw, lidat, ctrl)
    if (is.infinite(postnext)&is.infinite(postcur)) {
        tf <- tempfile("dumpmetrop", tmpdir=".", fileext = ".Rdata")
        li <- list(par=par, lidat=lidat, nam=nam, ctrl=ctrl, logposterior=logposterior,
                   postcur=postcur, postnext=postnext, optimiser=optimiser, i=NA,
                   nex=nex, prev=prev, parw=parw, functionError="MultiNormalStep")
        save(li, file=tf)
        msg <- paste0("Error with update of parameter ", nam,
                      ":\ninfinite current and next log-posterior.\n",
                      "All necessary information has been dumped in file", tf)
        stop(msg)
    }

    alpha <- min(c(1,exp(postnext-postcur)))
    un <- runif(1)<alpha
    accept <- as.numeric(un)
    if (un) {
        parm <- nex
        postcur <- postnext
    } else {
        parm <- prev
    }
    return(list(parm=parm, postcur=postcur, accept=accept))
}



#' @title Implements one Markov chain with the Metropolis algorithm
#' @export
#'
#' @details \code{GeneralSingleMetropolis} implements the Metropolis
#' algorithm for the MCMC fitting of Bayesian models.
#' \code{restartGSM} can be used to restart a model for \code{another}
#' iterations. \code{reloadGSM} can be used to return the GSMetro
#' object saved in a file by the backup mechanism of
#' \code{GeneralSingleMetropolis}. \code{defaultListUGSM} can be used
#' to generate default values for \code{listUpdating} (defaulting to
#' \code{"sis"} for parameters of length 1 and to \code{"mis"} for
#' vectors of parameters).  This list can then be manually altered
#' (e.g. changing the updating mechanism to "mns" for some
#' parameters).
#'
#' @param parInit named list with one element per parameter (the
#' names of the list correspond to the name of the parameters),
#' containing the initial values for the parameters
#' @param lidat a named list containing the data (the names of the
#' list correspond to the name of the variables used in logposterior,
#' combined with the parameters), required for the calculation of the
#' posterior.
#' @param logposterior a function for the calculation of the log
#' posterior for the current model. This function must rely only on
#' the data provided in \code{lidat}. This function must have three
#' parameters: (i) \code{par} is the list of parameters, (ii)
#' \code{lidat} is the list containing the data, and (iii)
#' \code{ctrl} is a named list with two elements: (i) an element
#' \code{namePar} is a character string containing the name of the
#' parameter in \code{par} that is updated, (ii) if this parameter is
#' a vector where each component is updated in turn (e.g. when
#' \code{MultipleIndependentSteps} is the updating mechanism),
#' \code{iter} indicates the element of this vector that is updated
#' (this element is ignored in other cases). It is possible to define
#' a log-posterior that does not make use of \code{ctrl} (i.e., which
#' returns the same value whatever the parameter that is updated),
#' but even in this case, the function should have an argument named
#' \code{ctrl}. Note that the argument \code{ctrl} is essentially
#' useful when the model contains many parameters
#' (e.g. overdispersion residuals).
#' @param lipum a named list containing the parameters of the
#' updating mechanisms (see the help page of the updating mechanisms
#' for a description of the required parameters).  The names of the
#' list correspond to the name of the parameters of the model
#' (i.e. \code{names(parInit)}).
#' @param listUpdating named list with one element parameter (the
#' names of the list correspond to the name of the parameters),
#' each one containing one of the following character strings:
#' (i) "mis", which implements the function \code{MultipleIndependentSteps}
#' as updating function for the corresponding parameter, (ii) "sis",
#' which implements the function \code{singleStep} as updating function
#' for the corresponding parameter, (iii) "mns", which implements
#' the function \code{MultiNormalStep}. If \code{NULL}, the function
#' \code{defaultListUGSM} is used internally to generate a default
#' value ("mis" or "sis").
#' @param liopt optionally, a vector containing the name of the
#' parameters for which the log-posterior is optimized (i.e. the
#' function logposterior returns a value that is proportional to the
#' log-posterior, and not the log-poserior itself, see the help page
#' of \code{singleStep}).
#' @param nrepet an integer giving the number of iterations required
#' for the MCMC algorithm.
#' @param nburnin an integer giving the number of burnin iterations
#' required
#' @param thinPar an integer giving the number of steps separating
#' the storage of two parameter vectors generated by the algorithm
#' @param thinAcc an integer giving the number of steps separating
#' the storage of the vector describing whether the proposal has been
#' accepted or not
#' @param verbose logical. Whether information should be displayed to
#' the user.
#' @param saveResults logical. Whether the results should be saved in a file.
#' @param saveEvery integer. When should the program save the list of
#' generated values (for backup)?
#' @param fileSave character string. The name of the Rdata file used to
#' save the results.
#' @param x an object of class GSMetrop.
#' @param resuParms an object of class GSMetrop.
#' @param another an integer giving the number of iterations required
#' for the MCMC algorithm.
#' @param \dots additional arguments to be passed to or from other
#' functions
#'
#'
#' @details This functions backs up regularly the state of the chain
#' (every \code{saveEvery} iterations) in a list of class
#' \code{GSMetrop} named \code{libackup}, saved in a file named
#' \code{fileSave}. This object can be retrieved with \code{reloadGSM}.
#'
#' @return A list of class "GSMetrop", where each element corresponds
#' to one parameter, and where each element is a matrix containing the
#' generated values as rows.
#' @author Clement Calenge, \email{clement.calenge@@oncfs.gouv.fr}
#' @examples
#' ## Generates a dataset:
#' x <- rnorm(100)
#' y <- 1+ 2*x + rnorm(100, sd = 0.1)
#'
#' ## Linear regression: estimate a, b, sigma in
#' ## y = a + b*x + epsilon
#' ## epsilon ~ dnorm(0, sigma)
#' ## We will model lsigma=log(sigma) to ensure it is >0
#'
#' ## Log-posterior
#' logp <- function(par, lidat, ctrl)
#' {
#'    ## Current values of the parameters
#'    a <- par$a
#'    b <- par$b
#'    sigma <- exp(par$lsigma) ## sigma=exp(log(sigma))
#'
#'    ## log-Prior:
#'    lprior <- dnorm(a, mean=0, sd=100, log=TRUE)
#'    lprior <- lprior + dnorm(b, mean=0, sd=100, log=TRUE)
#'    ## uniform prior on log-sigma
#'
#'    ## log-likelihood
#'    llike <- sum(dnorm(lidat$y, mean=lidat$x*b+a, sd=sigma, log=TRUE))
#'
#'    ## log-posterior
#'    return(llike+lprior)
#' }
#'
#' ## Elements required for the fit
#' parinit <- list(a=0, b=1, lsigma=0) ## initial values
#'
#' lipum <- list(a=0.1, b=0.1, lsigma=1) ## test standard deviation of the proposals
#'
#' lidat <- list(x=x, y=y) ## the data
#'
#' \dontrun{
#' gsm <- GeneralSingleMetropolis(parinit, lidat, logp, lipum,
#'                                nrepet=10000, saveResults = FALSE)
#' gsm
#'
#' ## Plot the chain
#' plot(gsm)
#'
#' ## We increase the burnin
#' gsm2 <- increaseBurnin(gsm, newBurnin=200)
#'
#' ## Ok
#' plot(gsm2)
#'
#' ## show a summary of the results
#' summary(gsm2, parameters=TRUE, accept=TRUE)
#'
#' ## Try to use a multinormal updating for a and b
#' ## we change the logposterior, and define theta
#' ## as c(a,b)
#' logp2 <- function(par, lidat, ctrl)
#' {
#'    ## Current values of the parameters
#'    theta <- par$theta
#'    a <- theta[1]
#'    b <- theta[2]
#'    sigma <- exp(par$lsigma)
#'
#'    ## log-Prior:
#'    lprior <- dnorm(a, mean=0, sd=100, log=TRUE)
#'    lprior <- lprior + dnorm(b, mean=0, sd=100, log=TRUE)
#'    ## uniform prior on sigma
#'
#'    ## log-likelihood
#'    llike <- sum(dnorm(lidat$y, mean=lidat$x*b+a, sd=sigma, log=TRUE))
#'
#'    ## log-posterior
#'    return(llike+lprior)
#' }
#'
#' ## We get the last value of the parameters generated in gsm2
#' pr <- getParameterVector(gsm2, nrow(gsm2[[1]]))
#'
#' ## We change the variables
#' parinit2 <- list(theta=c(pr$a, pr$b), lsigma=pr$lsigma) ## initial values
#' listUpdating2 <- list(theta="mns", lsigma="sis") ## now, theta is updated with mns
#' lipum2 <- list(theta=cov(cbind(gsm2$a, gsm2$b)), lsigma=1)
#'
#' ## And we start again, saving the results in ficxample231.Rdata
#' gsm3 <- GeneralSingleMetropolis(parinit2, lidat, logp2, lipum2, listUpdating2,
#'                                 nrepet=10000, fileSave="ficxample231.Rdata")
#'
#' gsm3
#'
#' plot(gsm3)
#'
#' summary(gsm3, parameters=TRUE, accept=TRUE)
#'
#' ## Demonstrate how to continue the iterations:
#' gsm4 <- restartGSM(gsm3, another=10000)
#' gsm4
#'
#' ## demonstrates how to reload a file:
#' gsm5 <- reloadGSM("ficxample231.Rdata")
#'
#' ## housekeeping
#' file.remove("ficxample231.Rdata")
#' }
GeneralSingleMetropolis <- function(parInit, lidat, logposterior, lipum, listUpdating=NULL,
                                    liopt = NULL, nrepet=999, nburnin=10, thinPar = 1,
                                    thinAcc = thinPar,
                                    verbose=TRUE, saveResults=TRUE, saveEvery=1000,
                                    fileSave="saveMetrop.Rdata")
{
    ## Starts duration
    startTime <- proc.time()

    ## Generates listUpdating
    if (is.null(listUpdating))
        listUpdating <- defaultListUGSM(parInit)

    ## Checks length
    if (length(parInit)!=length(listUpdating))
        stop("listUpdating does not have the same length as parInit")
    ## Checks names
    if (!all(names(parInit)==names(listUpdating)))
        stop("listUpdating has names different from the parInit")

    ## checks updating mechanisms
    if (!all(unlist(listUpdating%in%c("mis", "mns", "sis"))))
        stop("Not allowed updating mechanisms")

    ## checks arguments of logposterior
    ar <- names(formals(logposterior))
    if (!all((ar==c("par", "lidat", "ctrl"))))
        stop("The log-posterior must have three arguments named par, lidat and ctrl")

    ## Checks that the logposterior works for all parameters (at least for the
    ## first element)
    restmp <- lapply(1:length(parInit), function(i) {
                         ii <- try(logposterior(parInit, lidat,
                                                list(namePar=names(parInit)[i], iter=1)))
                         if (inherits(ii, "try-error"))
                             stop(paste("The evaluation of the log-posterior fails for parameter",
                                        names(parInit)[i]))
                     })

    lp <- logposterior(parInit, lidat,
                       list(namePar=names(parInit)[1], iter=1))
    if (is.infinite(lp))
        stop("bad starting values for the parameters (infinite log-posterior)")

    ## We pass the list to par
    par <- parInit

    ## The names of the parameters
    NamesP <- names(par)

    ## we transform the code for updating parameters into function names
    upd <- as.numeric(factor(sapply(names(par), function(x) listUpdating[[x]]),
                             levels=c("mis","mns","sis")))
    UpdatingMechanism <- c("MultipleIndependentSteps","MultiNormalStep","singleStep")[upd]

    ## List of generated parameters
    liparms <- list()

    ## List of accepted mechanisms
    liacc <- list()


    ## function used for saving
    saveMetrop <- function(r)
    {
        resuParms <- lapply(names(par), function(na) {
                                do.call(rbind,lapply(liparms, function(x) {
                                                         x[[na]]
                                                     }))
                            })
        resuAcc <- lapply(names(par), function(na) {
                              do.call(rbind,lapply(liacc, function(x) {
                                                       x[[na]]
                                                   }))
                          })
        names(resuParms) <- names(par)
        names(resuAcc) <- names(par)
        attr(resuParms, "AcceptationRate") <- resuAcc
        attr(resuParms, "thinPar") <- thinPar
        attr(resuParms, "thinAcc") <- thinAcc
        attr(resuParms, "nburnin") <- nburnin
        prestime <- proc.time()
        attr(resuParms, "duration") <- (prestime-startTime)
        attr(resuParms, "listUpdating") <- listUpdating
        attr(resuParms, "lidat") <- lidat
        attr(resuParms, "lipum") <- lipum
        attr(resuParms, "logposterior") <- logposterior
        attr(resuParms, "verbose") <- verbose
        attr(resuParms, "saveResults") <- saveResults
        attr(resuParms, "saveEvery") <- saveEvery
        attr(resuParms, "fileSave") <- fileSave
        if (!is.null(liopt))
            attr(resuParms,"liopt") <- liopt
        attr(resuParms, "fileSave") <- fileSave
        class(resuParms) <- "GSMetrop"

        libackup <- resuParms
        if (saveResults)
            save(libackup, file=fileSave)
        return(resuParms)
    }



    ## Main loop
    for (r in 1:nrepet) {
        if (verbose)
            cat("Iteration", r, "\r")

        im <- .iterMetrop(par, lidat, logposterior, lipum, NamesP, UpdatingMechanism, liopt)
        par <- im$par
        accept <- im$accept

        ## Storage
        if ((r%%thinPar==0)&(r>nburnin)) {
            liparms[[length(liparms)+1]] <- par
        }

        ## Acceptation
        if (r%%thinAcc==0&(r>nburnin)) {
            liacc[[length(liacc)+1]] <- accept
        }

        ## save?
        if ((r%%saveEvery==0)&(r>nburnin)&(!is.na(fileSave))) {
            resuParms <- saveMetrop(r)
        }
    }

    ## Results
    if (!is.na(fileSave)) {
        resuParms <- saveMetrop(r)
    }

    class(resuParms) <- "GSMetrop"
    return(resuParms)
}


#' @rdname GeneralSingleMetropolis
#' @export
restartGSM <- function(resuParms, another=1000)
{
    ## Checks
    if (!inherits(resuParms, "GSMetrop"))
        stop("incorrect class")

    ## No need to check the other arguments, the class
    ## "backupMetropolis" says it all

    ## Reassigns arguments
    liparms <- lapply(1:nrow(resuParms[[1]]), function(i) {
                          lapply(resuParms, function(x) x[i,])
                      })
    resuAcc <- attr(resuParms, "AcceptationRate")
    liacc <- lapply(1:nrow(resuAcc[[1]]), function(i) {
                        lapply(resuAcc, function(x) x[i,])
                    })

    resuAcc <- attr(resuParms, "AcceptationRate")
    thinPar <- attr(resuParms, "thinPar")
    thinAcc <- attr(resuParms, "thinAcc")
    nburnin <- attr(resuParms, "nburnin")
    liopt <- attr(resuParms, "liopt")
    nburnin <- attr(resuParms, "nburnin")
                                        #     (prestime-startTime) <- attr(resuParms, "duration")
    listUpdating <- attr(resuParms, "listUpdating")
    lidat <- attr(resuParms, "lidat")
    lipum <- attr(resuParms, "lipum")
    logposterior <- attr(resuParms, "logposterior")
    verbose <- attr(resuParms, "verbose")
    saveResults <- attr(resuParms, "saveResults")
    saveEvery <- attr(resuParms, "saveEvery")
    fileSave <- attr(resuParms, "fileSave")

    ## current value of the parameters
    par <- getParameterVector(resuParms, nrow(resuParms[[1]]))

    ## The names of the parameters
    NamesP <- names(par)

    ## we transform the code for updating parameters into function names
    upd <- as.numeric(factor(sapply(names(par), function(x) listUpdating[[x]]),
                             levels=c("mis","mns","sis")))
    UpdatingMechanism <- c("MultipleIndependentSteps","MultiNormalStep","singleStep")[upd]

    ## measure Time
    duration <- attr(resuParms, "duration")
    startTime <- proc.time()

    ## function used for saving
    saveMetrop <- function(r)
    {
        resuParms <- lapply(names(par), function(na) {
                                do.call(rbind,lapply(liparms, function(x) {
                                                         x[[na]]
                                                     }))
                            })
        resuAcc <- lapply(names(par), function(na) {
                              do.call(rbind,lapply(liacc, function(x) {
                                                       x[[na]]
                                                   }))
                          })
        names(resuParms) <- names(par)
        names(resuAcc) <- names(par)
        attr(resuParms, "AcceptationRate") <- resuAcc
        attr(resuParms, "thinPar") <- thinPar
        attr(resuParms, "thinAcc") <- thinAcc
        attr(resuParms, "nburnin") <- nburnin
        prestime <- proc.time()
        attr(resuParms, "duration") <- duration+(prestime-startTime)
        attr(resuParms, "listUpdating") <- listUpdating
        attr(resuParms, "lidat") <- lidat
        attr(resuParms, "lipum") <- lipum
        attr(resuParms, "logposterior") <- logposterior
        attr(resuParms, "verbose") <- verbose
        attr(resuParms, "saveEvery") <- saveEvery
        attr(resuParms, "saveResults") <- saveResults
        attr(resuParms, "fileSave") <- fileSave
        if (!is.null(liopt))
            attr(resuParms,"liopt") <- liopt
        class(resuParms) <- "GSMetrop"

        libackup <- resuParms
        if (saveResults)
            save(libackup, file=fileSave)
        return(resuParms)
    }


    ## Main loop
    for (r in 1:another) {
        if (verbose)
            cat("Iteration", r, "\r")

        im <- .iterMetrop(par, lidat, logposterior, lipum, NamesP, UpdatingMechanism,
                          liopt)
        par <- im$par
        accept <- im$accept

        ## Storage
        if (r%%thinPar==0) {
            liparms[[length(liparms)+1]] <- par
        }

        ## Acceptation
        if (r%%thinAcc==0) {
            liacc[[length(liacc)+1]] <- accept
        }

        ## save?
        if ((r%%saveEvery==0)&(!is.na(fileSave))) {
            resuParms <- saveMetrop(r)
        }
    }

    ## Results
    if (!is.na(fileSave)) {
        resuParms <- saveMetrop(r)
    }

    return(resuParms)
}


#' @rdname GeneralSingleMetropolis
#' @export
defaultListUGSM <- function(parInit)
{
    if (length(names(parInit))<1)
        stop("parInit should be a named list, with one element per parameter")
    listUpdating <- lapply(1:length(parInit), function(i) {
                               if (length(parInit[[i]])==1)
                                   return("sis")
                               return("mis")
                           })
    names(listUpdating) <- names(parInit)
    return(listUpdating)
}


#' @rdname GeneralSingleMetropolis
#' @export
print.GSMetrop <- function(x, ...)
{
    if (!inherits(x, "GSMetrop"))
        stop("x should inherit the class \"GSMetrop\".")
    cat("Object of class \"GSMetrop\"\n\n")
    su <- summary(x)
    cat(length(x), "parameters or vectors of parameters in the model:\n")
    print(names(x))
    cat("\n")
    if ((attr(x,"nburnin"))<attr(x,"thinPar")) {
        nite <- nrow(x[[1]])*attr(x,"thinPar")
    } else {
        nite <- nrow(x[[1]])*attr(x,"thinPar")+(attr(x,"nburnin"))%/%(attr(x,"thinPar"))
    }
    cat(nite, "iterations (including",
        attr(x,"nburnin"), " burn-in samples)\n")
    cat("Thin parameters every:", attr(x,"thinPar"), "iterations\n")
    cat("Number of iterations stored:", nrow(x[[1]]),"\n")
    cat("Calculation took",  round(attr(x,"duration")[3]),"seconds (total time)\n\n")

    if (nrow(su[[1]])<10) {
        cat("Distribution of the acceptation rates:\n")
        print(su[[1]][,2])
    } else {
        cat("Summary of the acceptation rates: \n")
        print(summary(su[[1]]$AcceptationRate))
    }
}



#' @title Summary information on a 'GSMetrop' object
#' @export
#'
#' @details Calculate acceptation rates, mean and SD from a
#' 'GSMetrop' object
#'
#' @param object an object of class \code{GSMetrop}
#' @param which the name of the parameters to summarize
#' @param accept logical. whether to calculate acceptation rates
#' @param parameters logical. whether to calculate mean and sd
#' @param \dots additional parameters to be passed
#'
#' @seealso \code{\link{GeneralSingleMetropolis}} for examples of use.
#' @author Clement Calenge, \email{clement.calenge@@oncfs.gouv.fr}
summary.GSMetrop <- function(object, which=names(object),
                             accept=TRUE, parameters=FALSE, ...)
{
    if (!inherits(object, "GSMetrop"))
        stop("object should inherit the class \"GSMetrop")


    resu <- list()
    if (accept) {
        ac <- attr(object, "AcceptationRate")[which]
        moya <- unlist(lapply(ac, function(y) apply(y,2,mean)))
        naac <- unlist(lapply(1:length(ac), function(i) {
                                  if (ncol(ac[[i]])>1)
                                      return(paste0(names(ac)[i], "[", 1:ncol(ac[[i]]),"]"))
                                  return(names(ac[i]))
                              }))
        resu$AcceptationRate <- data.frame(Parameter=naac, AcceptationRate=moya)
        row.names(resu$AcceptationRate) <- 1:nrow(resu$AcceptationRate)
    }
    if (parameters) {
        object <- object[which]
        napa <- unlist(lapply(1:length(object), function(i) {
                                  if (ncol(object[[i]])>1)
                                      return(paste0(names(object)[i], "[",
                                                    1:ncol(object[[i]]),"]"))
                                  return(names(object[i]))
                              }))
        moyp <- unlist(lapply(object, function(y) apply(y,2,mean)))
        sdp <- unlist(lapply(object, function(y) apply(y,2,sd)))
        resu$Parameters <- data.frame(Parameter=napa, Mean=moyp, SD=sdp)
        row.names(resu$Parameters) <- 1:nrow(resu$Parameters)
    }
    return(resu)
}


#' @title graphical display of the result of a Metropolis run
#' @export
#' @details plot the mixing of a metropolis run
#'
#' @param x an object of class \code{GSMetrop}
#' @param which vector of name of parameters to be displayed
#' @param show the number of plot to show on the graphical device
#' @param \dots additional arguments to be passed to and from other
#' functions
#'
#' @seealso \code{\link{GeneralSingleMetropolis}} for examples of use.
#' @author Clement Calenge, \email{clement.calenge@@oncfs.gouv.fr}
#' @importFrom grDevices n2mfrow
plot.GSMetrop <- function(x, which=names(x), show=100, ...)
{
    if (!inherits(x, "GSMetrop"))
        stop("object should inherit the class \"GSMetrop")

    x <- x[which]
    npl <- sum(sapply(x,ncol))
    ngr <- min(c(show,npl))
    opar <- par(mfrow = n2mfrow(ngr),
                ask=show<npl, mar=c(0,0,2,0))
    on.exit(par(opar))
    na <- lapply(1:length(x), function(i) {
                     if (ncol(x[[i]])>1)
                         return(paste0(names(x)[i], "[", 1:ncol(x[[i]]),"]"))
                     return(names(x[i]))
                 })
    for (i in 1:length(na)) {
        for (j in 1:length(na[[i]])) {
            plot(x[[i]][,j], ty="l", main=na[[i]][j], axes=FALSE)
            box()
        }
    }
}




#' @title transformations used in the package
#' @export
#'
#' @details \code{logit} calculates the logit of the argument, and
#' \code{invlogit} calculates the inverse logit transform.
#'
#' @param x parameter of a model.
#' @seealso \code{\link{GeneralSingleMetropolis}} for examples of use.
#' @author Clement Calenge, \email{clement.calenge@@oncfs.gouv.fr}
logit <- function(x) return(log(x/(1-x)))


#' @rdname logit
#' @export
invlogit <- function(x) return(exp(x)/(1+exp(x)))


#' @rdname GeneralSingleMetropolis
#' @export
reloadGSM <- function(fileSave)
{
    .myDataEnv <- new.env(parent=emptyenv())
    load(fileSave, envir=.myDataEnv)
    if (!exists("libackup", envir=.myDataEnv))
        stop("incorrect file")
    libackup <- get("libackup", envir=.myDataEnv)
    if (!inherits(libackup, "GSMetrop"))
        stop("incorrect file")
    return(libackup)
}



#' @title Various functions used to manipulate results of the Metropolis Algorithm
#' @export
#' @details \code{getParameterVector} extracts the parameter vector
#' corresponding to a given iteration.  \code{increaseBurnin} can be
#' used to increase a posteriori the size of the burnin
#' sample. \code{nit} returns the number of stored iterations in a
#' GSMetrop object.
#'
#' @param resuParms an object of class \code{GSMetrop}
#' @param r iteration to be extracted
#' @param newBurnin new size of the burnin sample (in number of
#' *performed* iterations: in other words, if the objects results from
#' a calculation which carried out 5000 iterations thinned every 5
#' iterations, with a 10 iterations burnin period, then the original
#' object will contain 990 stored iterations. Setting
#' \code{newBurnin=500} will remove the 500 first performed
#' iterations, i.e. the first (500 - 10 old burnin)/(thin=5) = 98
#' iterations stored in the object).
#' @seealso \code{\link{GeneralSingleMetropolis}} for examples of use.
#' @author Clement Calenge, \email{clement.calenge@@oncfs.gouv.fr}
getParameterVector <- function(resuParms, r)
{
    if (!inherits(resuParms, "GSMetrop"))
        stop("incorrect class")
    par <- lapply(resuParms, function(x) x[r,])
    names(par) <- names(resuParms)
    return(par)
}

#' @rdname getParameterVector
#' @export
nit <- function(resuParms)
{
    if (!inherits(resuParms, "GSMetrop"))
        stop("incorrect class")
    return(nrow(resuParms[[1]]))
}

#' @rdname getParameterVector
#' @export
increaseBurnin <- function(resuParms, newBurnin=attr(resuParms,"nburnin"))
{
    if (!inherits(resuParms, "GSMetrop"))
        stop("incorrect class")
    if (attr(resuParms,"nburnin")>newBurnin)
        stop("newBurnin cannot be lower than the present value of nburnin")
    if (((nrow(resuParms[[1]])*attr(resuParms, "thinPar"))+attr(resuParms,"nburnin"))<=(newBurnin))
        stop("newBurnin too large")
    nf <- 1+(newBurnin-attr(resuParms,"nburnin"))%/%attr(resuParms,"thinPar")
    attr(resuParms,"nburnin") <- newBurnin
    nr <- nrow(resuParms[[1]])
    for (i in 1:length(resuParms))
        resuParms[[i]] <- resuParms[[i]][nf:nr,, drop=FALSE]
    return(resuParms)
}


#' @title Combine Several Single Metropolis Objects Into One CMM objects
#' @export
#' @aliases CMM
#' @details \code{c.GSMetrop} transforms several objects of class
#' \code{GSMetrop} into an object of class \code{CMM} (combined
#' multiple Metropolis).
#'
#' @param \dots a list of objects of class \code{GSMetrop}
#' @return an object of class "CMM".  If only one object is passed,
#' the function returns the input object without transformation.
#' @seealso \code{\link{GeneralSingleMetropolis}}
#' @author Clement Calenge, \email{clement.calenge@@oncfs.gouv.fr}
#' @examples
#' #############################################
#' ##
#' ## First start with the same examples as in
#' ## the help page of GeneralSingleMetropolis
#' ##
#' ## Generates a dataset:
#' x <- rnorm(100)
#' y <- 1+ 2*x + rnorm(100, sd = 0.1)
#'
#' ## Linear regression: estimate a, b, sigma in
#' ## y = a + b*x + epsilon
#' ## epsilon ~ dnorm(0, sigma)
#'
#' ## Log-posterior
#' logp <- function(par, lidat, ctrl)
#' {
#'    ## Current values of the parameters
#'    a <- par$a
#'    b <- par$b
#'    sigma <- exp(par$lsigma)
#'
#'    ## log-Prior:
#'    lprior <- dnorm(a, mean=0, sd=100, log=TRUE)
#'    lprior <- lprior + dnorm(b, mean=0, sd=100, log=TRUE)
#'    ## uniform prior on sigma
#'
#'    ## log-likelihood
#'    llike <- sum(dnorm(lidat$y, mean=lidat$x*b+a, sd=sigma, log=TRUE))
#'
#'    ## log-posterior
#'    return(llike+lprior)
#' }
#'
#' ## Elements required for the fit
#' parinit <- list(a=0, b=1, lsigma=1) ## initial values
#' listUpdating <- list(a="sis", b="sis", lsigma="sis") ## all elements are updated with sis
#' lidat <- list(x=x, y=y)
#' lipum <- list(a=0.1, b=0.1, lsigma=1) ## test standard deviation of the proposals
#'
#' \dontrun{
#' #############################################
#' ##
#' ## Then demonstrate the class CMM
#' ##
#' ## Three runs
#' gsm1 <- GeneralSingleMetropolis(parinit, lidat, logp, lipum,
#'                                 nrepet=5000, saveResults=FALSE, nburnin=500)
#' gsm2 <- GeneralSingleMetropolis(parinit, lidat, logp, lipum,
#'                                 nrepet=5000, saveResults=FALSE, nburnin=500)
#' gsm3 <- GeneralSingleMetropolis(parinit, lidat, logp, lipum,
#'                                 nrepet=5000, saveResults=FALSE, nburnin=500)
#'
#' ## Combine these three runs
#' ae <- c(gsm1, gsm2, gsm3)
#'
#' ## Show various information about these models
#' ae
#' summary(ae)
#' plot(ae)
#'
#' ## Extraction features:
#' ae[1:2]
#' ae[2]
#'
#' ## work with the package coda:
#' library(coda)
#' raftery.diag(tocoda(gsm1))
#'
#' gelman.diag(tocoda(ae))
#' }
#'
c.GSMetrop <- function(...)
{
    li <- list(...)
    ## correct class?
    if (any(sapply(li,class)!="GSMetrop"))
        stop("objects should be of class GSMetrop")
    if (length(li)==1)
        return(li[[1]])
    ## Same number of attributes?
    liat <- lapply(li, function(x) {
        attr <- attributes(x)
        attr$AcceptationRate <- NULL
        attr$duration <- NULL
        return(attr)
    })
    if (any(sapply(li, length)!=length(li[[1]])))
        stop("all objects do not have the same number of attributes")
    ## Same names of attributes?
    if (!all(apply(sapply(liat, names),1, function(x) all(x==x[1]))))
        stop("not the same attribute names")
    ## same values for the attributes
    liatb <- lapply(liat, function(x) {x$logposterior <- NULL;return(x)})
    liatb <- lapply(liatb, function(x) {x$fileSave <- NULL;return(x)})
    tmp <- lapply(1:length(liatb[[1]]), function(i) {
        la <- lapply(1:length(liatb), function(j) liatb[[j]][[i]])
        if (!all(sapply(1:length(la), function(j) identical(la[[j]],la[[1]]))))
            stop(paste0("Non consistent value for attribute ", names(liatb[[1]])[i]))
    })
    ## same values for the posterior
    tmp <- lapply(liat, function(x) x$logposterior)
    if (!all(sapply(tmp, function(y) body(y)==body(tmp[[1]]))))
        stop("different log-posterior")


    ## ok, combines
    cmb <- lapply(li, function(x) {
        acc <- attr(x,"AcceptationRate")
        dur <- attr(x,"duration")
        fs <- attr(x,"fileSave")
        nam <- names(x)
        attributes(x) <- NULL
        names(x) <- nam
        attr(x,"duration") <- dur
        attr(x,"AcceptationRate") <- acc
        attr(x,"fileSave") <- fs
        return(x)
    })
    liat <- lapply(li, function(x) {
        attr(x,"AcceptationRate") <- NULL
        attr(x,"duration") <- NULL
        attr(x,"fileSave") <- NULL
        names(x) <- NULL
        return(attributes(x))
    })
    attributes(cmb) <- liat[[1]]
    class(cmb) <- "CMM"
    return(cmb)
}



#' @title Extract elements from a CMM object
#' @export
#' @details \code{Extract} function can be used to extract one or
#' several elements from an object of class \code{CMM} (combined
#' multiple Metropolis).
#'
#' @param x an object of class \code{CMM}
#' @param i a vector of positive integer values used as indices of the
#' elements to extract.
#' @return an object of class "CMM".  If only one object is passed,
#' the function returns the input object without transformation.
#' @note negative indices (for element removal) are not allowed.
#' @seealso \code{\link{c.GSMetrop}} for examples
#' @author Clement Calenge, \email{clement.calenge@@oncfs.gouv.fr}
"[.CMM" <- function(x, i)
{
    if (!inherits(x, "CMM"))
        stop("x should be of class \"CMM\"")
    if (max(i)>length(x)) {
        stop("too large indices")
    }
    if (any(i<=0))
            stop("negative or zero indices not allowed")

    if (length(i)>1) {
        atr <- attributes(x)
        x <- unclass(x)
        y <- x[i]
        attributes(y) <- atr
        return(y)
    }
    aa <- attributes(x)
    li <- x[[i]]
    ab <- attributes(li)
    attributes(li) <- aa
    for (j in 1:length(ab))
        attr(li, names(ab)[j]) <- ab[[j]]
    class(li) <- "GSMetrop"
    return(li)
}


#' @title Prints a CMM object
#' @export
#' @details \code{print.CMM} function displays a short summary of an
#' object of class \code{CMM} (combined multiple Metropolis).
#'
#' @param x an object of class \code{CMM}
#' @param \dots additional arguments to be passed to other functions.
#' @seealso \code{\link{c.GSMetrop}} for examples
#' @author Clement Calenge, \email{clement.calenge@@oncfs.gouv.fr}
print.CMM <- function(x, ...)
{
    if (!inherits(x, "CMM"))
        stop("x should inherit the class \"CMM\".")
    cat("Object of class \"CMM\"\n\n")
    cat("Number of chains:", length(x),"\n")
    cat(length(x), "parameters or vectors of parameters in the model:\n")
    print(names(x[[1]]))
    cat("\n")
    cat(nrow(x[[1]][[1]])*attr(x,"thinPar")+attr(x,"nburnin"), "iterations per chain (including",
        attr(x,"nburnin"), " burn-in samples)\n")
    cat("Thin parameters every:", attr(x,"thinPar"), "iterations\n")
    cat("Number of iterations stored for each chain:", nrow(x[[1]][[1]]),"\n")
    cat("Calculation took",  round(sum(sapply(1:length(x), function(i) attr(x[[i]],"duration")[3]))),
                                   "seconds (total time over all chains)\n\n")
}

#' @title graphical display of the result of several Metropolis runs
#' @export
#' @details plot the mixing of several metropolis runs
#'
#' @param x an object of class \code{CMM}
#' @param which vector of name of parameters to be displayed
#' @param show the number of plot to show on the graphical device
#' @param col a character vector giving the colors of the different chains
#' @param \dots additional arguments to be passed to and from other
#' functions
#'
#' @seealso \code{\link{c.GSMetrop}} for examples of use.
#' @author Clement Calenge, \email{clement.calenge@@oncfs.gouv.fr}
#' @importFrom grDevices n2mfrow rainbow
plot.CMM <- function(x, which=names(x[1]), show=100, col=rainbow(length(x)), ...)
{
    if (!inherits(x, "CMM"))
        stop("object should inherit the class \"CMM")

    xt <- x[1][which]

    npl <- sum(sapply(xt,ncol))
    ngr <- min(c(show,npl))
    opar <- par(mfrow = n2mfrow(ngr),
                ask=show<npl, mar=c(0,0,2,0))
    on.exit(par(opar))
    na <- lapply(1:length(xt), function(i) {
                     if (ncol(xt[[i]])>1)
                         return(paste0(names(xt)[i], "[", 1:ncol(xt[[i]]),"]"))
                     return(names(xt[i]))
                 })
    for (i in 1:length(na)) {
        for (j in 1:length(na[[i]])) {
            gg <- range(unlist(lapply(x, function(y) y[[which[i]]][,j])))
            plot(x[[1]][[i]][,j], ty="n", main=na[[i]][j], axes=FALSE, ylim = gg)
            tmp <- lapply(1:length(x), function(k) lines(x[[k]][[which[i]]][,j], col=col[k]))
            box()
        }
    }
}


#' @title Summary information on a 'CMM' object
#' @export
#'
#' @details Calculate mean and SD of parameters from a
#' 'CMM' object
#'
#' @param object an object of class \code{CMM}
#' @param which the name of the parameters to summarize
#' @param perchain logical. whether the statistics should be calculated for each chain separately or globally
#' @param \dots additional parameters to be passed
#'
#' @return a data.frame with the required statistics.
#'
#' @seealso \code{\link{c.GSMetrop}} for examples of use.
#' @author Clement Calenge, \email{clement.calenge@@oncfs.gouv.fr}
summary.CMM <- function(object, which=names(object[1]), perchain=TRUE, ...)
{
    if (!inherits(object, "CMM"))
        stop("object should inherit the class \"CMM\"")
    if (perchain) {
        ta <- lapply(1:length(object), function(i) {
            su <- summary(object[i], which=which, accept = FALSE,
                          parameters = TRUE)$Parameters
            names(su) <- paste0("Chain.",i,".", names(su))
            return(su)
        })
        ta2 <- cbind.data.frame(Parameter=ta[[1]][,1], lapply(ta,function(x) x[,-1]))
        return(ta2)
    } else {
        objectb <- lapply(1:length(object[[1]]), function(i)
            do.call(rbind, lapply(object, function(x) x[[i]])))
        names(objectb) <- names(object[[1]])
        object <- objectb[which]
        napa <- unlist(lapply(1:length(object), function(i) {
                                  if (ncol(object[[i]])>1)
                                      return(paste0(names(object)[i], "[",
                                                    1:ncol(object[[i]]),"]"))
                                  return(names(object[i]))
                              }))
        moyp <- unlist(lapply(object, function(y) apply(y,2,mean)))
        sdp <- unlist(lapply(object, function(y) apply(y,2,sd)))
        resu <- data.frame(Parameter=napa, Mean=moyp, SD=sdp)
        row.names(resu) <- 1:nrow(resu)
        return(resu)
    }

}



#' @title Extracts the name of the parameters from a model
#' @export
#'
#' @details Extracts the name of the parameters from a \code{CMM} or
#' \code{GSMetrop} object
#'
#' @param x an object of class \code{CMM} or \code{GSMetrop}
#'
#' @return a character vector giving the names of the parameters.
#'
#' @seealso \code{\link{c.GSMetrop}} for examples of use.
#' @author Clement Calenge, \email{clement.calenge@@oncfs.gouv.fr}
param <- function(x)
{
    if ((!inherits(x, "GSMetrop"))&(!inherits(x,"CMM")))
        stop("x should inherits GSMetrop or CMM")
    if (inherits(x, "GSMetrop"))
        return(names(x))
    return(names(x[[2]]))
}



#' @title Exports object toward the coda package
#' @export
#'
#' @details \code{tocoda} exports objects created with the package
#' metroponcfs toward classes suitable for further analysis with the
#' coda package.
#'
#' @param x an object of class \code{CMM} or \code{GSMetrop}
#' @param which a character vector indicating which parameters should be exported
#'
#' @return if \code{x} is an object of class GSMetrop, an object of
#' class \code{mcmc}. If \code{x} is an object of class CMM, an object
#' of class {mcmc.list}
#'
#' @seealso \code{\link{c.GSMetrop}} for examples of use.
#' @author Clement Calenge, \email{clement.calenge@@oncfs.gouv.fr}
#' @importFrom coda mcmc mcmc.list
tocoda <- function(x, which=param(x))
{
    if ((!inherits(x, "GSMetrop"))&(!inherits(x,"CMM")))
        stop("x should inherits GSMetrop or CMM")

    if (inherits(x,"GSMetrop")) {
        xt <- x[which]
        aa <- do.call(cbind,xt)
        napa <- unlist(lapply(1:length(xt), function(i) {
                                  if (ncol(xt[[i]])>1)
                                      return(paste0(names(xt)[i], "[",
                                                    1:ncol(xt[[i]]),"]"))
                                  return(names(xt)[i])
                              }))
        colnames(aa) <- napa
        ## calcul du start
        if ((attr(x,"nburnin")+1)<attr(x,"thinPar")) {
            start <- attr(x,"thinPar")
        } else {
            start <- (((attr(x, "nburnin"))%/%attr(x,"thinPar"))+1)*attr(x,"thinPar")
        }
        return(coda::mcmc(aa, start=start,
                          thin=attr(x,"thinPar")))

    } else {
        return(do.call(coda::mcmc.list,lapply(1:length(x), function(i) {
            tocoda(x[i], which=which)
        })))
    }
}


#' @title Deviance Information Criterion
#' @export
#'
#' @details Calculates the deviance information criterion of a model.
#'
#' @param x an object of class \code{CMM} or \code{GSMetrop}
#' @param loglikelihood a function for the calculation of the log
#' likelihood for the current model. This function must rely only on
#' the data provided in \code{lidat}. This function must have three
#' parameters: (i) \code{par} is the list of parameters, (ii)
#' \code{lidat} is the list containing the data, and (iii) \code{ctrl}
#' is a named list (see the help page of
#' \code{GeneralSingleMetropolis}). Note that the argument \code{ctrl}
#' is not used by the function here.
#'
#' @return the value of the DIC.
#'
#' @seealso \code{\link{GeneralSingleMetropolis}}.
#' @author Clement Calenge, \email{clement.calenge@@oncfs.gouv.fr}
DIC <- function(x, loglikelihood)
{
    if ((!inherits(x, "GSMetrop"))&(!inherits(x,"CMM")))
        stop("x should inherits GSMetrop or CMM")
    if (inherits(x, "CMM")) {
        dev <- unlist(lapply(1:length(x), function(i) {
                                 mo <- x[i]
                                 dev <- -2*sapply(1:nit(mo), function(i) {
                                                      loglikelihood(getParameterVector(mo,i),
                                                                    attr(mo, "lidat"))
                                                  })
                             }))
    } else {
        dev <- -2*sapply(1:nit(x), function(i) {
                             loglikelihood(getParameterVector(x,i),
                                           attr(x, "lidat"))
                         })
    }
    return(mean(dev)+var(dev)/2)
}
