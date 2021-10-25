functions{
    vector ode_base ( real t,
                      vector y,
                      real lamdba,
                      real d,
                      real beta,
                      real a,
                      real k,
                      real u)
                      
                      
                      vector[y] dydt;
                      
                      dydt[1] = lambda  - d * y[1] - beta * y[1] * y[3]  ;
                      dydt[2] = - beta * y[1] * y[3] - a * y[2];
                      dydt[3] =  k * [2] - u * y[3];
                      
                      return(dydt);
                      
}

data {
  int<lower=0> N;
  real t[N];
  vector[3] y0;

}
