## Works Cited:
## Construction of a genetic toggle switch in Escherichia coli Timothy S. Gardner*2, Charles R. Cantor* & James J. Collins*2 * Department of Biomedical Engineering, 2 Center for BioDynamics and 3 Center for Advanced Biotechnology, Boston University, 44 Cummington Street, Boston, Massachusetts 02215, USA
library(FME)

pr <- c(Ku=2, Kv=2, Cv=2, Cu=2)
st <- c(U=1, V=1)
state <- st
params <- pr
sts <- lapply(seq(0,1,.1), function(x) return(c(U=as.double(x), V=as.double(x+1))))

solveBistable <- function(state, params) {
  derivs <- function(t, state, params) {
    with(as.list(c(state, params)), {
      dU = ( Ku / (1+V^Cv) ) - U
      dV = ( Kv / (1+U^Cu) ) - V   
      return(list(c(dU, dV)))
    })
  }
  times <- seq(0,10,0.1)
  out <- ode(y = state, parms = params, times = times, func = derivs)
  as.data.frame(out)
}
solve.states <- function(states=c(), parms) {
   out <- lapply(states, function(state){
    out <- solveBistable(state, parms)
    return(out)
  })
  return(c(out))
}

solve.inducer <- function(course=function(t){return(0)}, times=seq(0,10,0.1)) {
  derivs <- function(t, state, params) {
    with(as.list(c(state, params)), {
      dU = ( Ku / (1+V^(Cv)) ) - U
      dV = ( Kv / (1+U^(Cu - course(t))) ) - V   
      return(list(c(dU, dV)))
    })
  }
  out <- ode(y = state, parms = params, times = times, func = derivs)
  return(as.data.frame(out))
}
solve.iptg <- function(state, params, course, times) {
    derivs <- function(t, state, params) {
      with(as.list(c(state, params)), {
        dU = ( Ku / (1 + course(t)) ) - U
        dV = ( Kv / (1+U^Cu) ) - V
        return(list(c(dU, dV)))
      })
    }
    out <- ode(y = state, parms = params, times = times, func = derivs)
    return(as.data.frame(out))
}
solve.states.iptg <- function(states, parms, course, times) {
  out <- lapply(states, function(state){
    out <- solve.iptg(state, parms, course, times)
    return(out)
  })
  return(c(out))
}
solve.g <- function(course, times) {
  derivs <- function(t, state, params) {
    with(as.list(c(state, params)), {
      dU = ( Ku / (1+(V+course(t))) ) - U
      dV = ( Kv / (1+(U+course(t))) ) - V
      return(list(c(dU, dV)))
    })
  }
  out <- ode(y = state, parms = params, times = times, func = derivs)
  return(as.data.frame(out))
}
course.inducer <- function(t) {
  if(t>1 & t<5){
    return(t^2 * 0.5)
  }else{
    return(0)
  }
}
course.iptg <- function(t) {
  if(t>2 & t<5){
    return(-t/5)
  }else{
    return(0)
  }
}
c.iptg <- function(t) {
  if(t>2 & t<5){
    return(-t*0.1)
  }else{
    return(0)
  }
}

par(mfrow=c(1,1))
out_inducer <- solve.inducer(course=course.inducer)
plot(out_inducer$time, out_inducer$V, col='red', type='l')  

par(mfrow=c(2,2))
out <- solve.iptg(state=state, params=params, course=c.iptg, times=seq(1,10,0.1))
plot(out$time, out$V, type='l')
plot(out$time, out$U, type='l')

out <- solve.states(sts, pr)

#U plots
plot(out[[1]]$time, out[[1]]$U, type='l')
for(o in c(2:length(out))){
  lines(out[[o]]$time, out[[o]]$U, col='red')
}
              
#V plots
plot(out[[1]]$time, out[[1]]$V, type='l')

for(o in c(2:length(out))){
  lines(out[[o]]$time, out[[o]]$V, col='blue')
}
              
#nullc lines
plot(out[[1]]$V, out[[1]]$U, type='l', col ='blue')
for(o in c(2:length(out))){
  lines(out[[o]]$V, out[[o]]$U, type='l', col='blue')
  
}
              
##--
out <- solveBistable(st, pr)
plot(out$time, out$U, type='l', col='red')
plot(out$time, out$V, type='l', col='blue')

##--
par(mfrow=c(2,2))
out <- solve.states.iptg(sts, pr, c.iptg, seq(1,10,0.1))
plot(out[[1]]$time, out[[1]]$U, type='l')
for(o in c(2:length(out))){
  lines(out[[o]]$time, out[[o]]$U, col='red')
}
              
#V plots
plot(out[[1]]$time, out[[1]]$V, type='l')

for(o in c(2:length(out))){
  lines(out[[o]]$time, out[[o]]$V, col='blue')
}
              
#nullc lines
plot(out[[1]]$V, out[[1]]$U, type='l', col ='blue')
for(o in c(2:length(out))){
  lines(out[[o]]$V, out[[o]]$U, type='l', col='blue')
  
}

