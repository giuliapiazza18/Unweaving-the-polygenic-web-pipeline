getArgs <- function(args)
{
  if (length(args)>0)
  {
    isqgraph <- sapply(args,function(x)"qgraph"%in%class(x))
    argLists <- c(lapply(args[isqgraph],'[[','Arguments'),lapply(args[isqgraph],'[','layout'))
    args <- args[!isqgraph]
    newArgs <- lapply(argLists,getArgs)
    for (l in newArgs) args <- c(args,l[!names(l)%in%names(args)])
  }
  return(args)
}

isColor <- function(x) {
  sapply(x, function(X) {
    if (!is.logical(X)) tryCatch(is.matrix(col2rgb(X)), 
                                 error = function(e) FALSE) else FALSE
  })
}


Fade <- function(col,alpha,bg)
{
  # col = color to fade
  # bg = color to fade to
  # alpha = inverse transparency, 1 = fully visable, 0 = fully transparant.
  
  if (missing(bg)) bg <- par("bg")
  if (length(bg)!=1) stop("'bg' must be of length 1")
  if (length(alpha)==1) alpha <- rep(alpha,length(col))
  if (length(col)==1) col <- rep(col,length(alpha))
  if (length(col)!=length(alpha)) stop("Length of 'col' not equal to length of 'alpha'")
  
  n <- length(col)
  
  rgbCols <- col2rgb(col)/255
  rgbBG <-  col2rgb(bg)/255
  
  colAlpha <- col2rgb(col,alpha=TRUE)[4,]/255
  
  Mix <- rgbCols*rep(alpha,each=3) + rgbBG%*%t(1-alpha)
  
  return(rgb(Mix[1,],Mix[2,],Mix[3,],colAlpha))
}


# Map user space to inches space:
usr2inX <- function(x)
{
  usr <- par("usr")
  pin <- par("pin")
  (x-usr[1])/(usr[2]-usr[1]) * pin[1]
}

usr2inY <- function(x)
{
  usr <- par("usr")
  pin <- par("pin")
  (x-usr[3])/(usr[4]-usr[3]) * pin[2]
}

# Same but about origin (for atan2):
usr2inX2 <- function(x)
{
  usr <- par("usr")
  pin <- par("pin")
  x/(usr[2]-usr[1]) * pin[1]
}

usr2inY2 <- function(x)
{
  usr <- par("usr")
  pin <- par("pin")
  x/(usr[4]-usr[3]) * pin[2]
}
atan2usr2in <- function(x,y) atan2(usr2inX2(x),usr2inY2(y))%%(2*pi)

# Map inches space to user space:
in2usrX <- function(x)
{
  usr <- par("usr")
  pin <- par("pin")
  usr[1] + x/pin[1] * (usr[2] - usr[1])
}

in2usrY <- function(x)
{
  usr <- par("usr")
  pin <- par("pin")
  usr[3] + x/pin[2] * (usr[4] - usr[3])
}

## Find perpundicular poin to quantile of line:
PerpMid <- function(xy0,xy1,ang=1,cex=1,q=0.5)
{
  # Change xy0 to quantile:
  xy0 <- xy0 + q * (xy1 - xy0)
  
  # Fixed inches size:
  cexIn <-  cex * 0.025 * sqrt(sum(par("pin")^2))
  
  # Rotate about origin:
  xyr <- matrix(c(0,ang,-ang,0),2,2) %*% (c(usr2inX(xy1[1]),usr2inY(xy1[2])) - c(usr2inX(xy0[1]),usr2inY(xy0[2])))  
  
  # Rescale:
  xyr <- xyr * cexIn/sqrt(sum(xyr^2))
  
  # Add origin:
  xyr <- c(usr2inX(xy0[1]),usr2inY(xy0[2])) + xyr
  
  # Map to usr and return:
  return(c(in2usrX(xyr[1]),in2usrY(xyr[2])))
}

