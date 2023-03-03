addlog <- function(a,b)
{
  if(b > a + 200)
    return(b)
  else if(a > b + 200)
    return(a)
  else
    return(a + log1p(exp(b-a)));
}


u<-runif(10, 2, 13)
u1<-u/sum(u)
sum(u1)
v <- log(u)
w <- log(0)
for(i in 1:10)
  w <- addlog(w, v[i])
u2<-exp(v-w)
cbind(u1,u2)
