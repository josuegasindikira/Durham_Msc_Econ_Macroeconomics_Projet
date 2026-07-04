// Summarative assignment

var c k y inv w r repla a zx zc; 
//With R the replacement price of capital
// 3 exogenous shocks: tfp, investment-specific technology, preference
varexo eps_a eps_zx eps_zc;

parameters 

beta alpha sigma delta rho_a rho_zx rho_zc abar zxbar zcbar;

abar   = 1;
zxbar  = 1;
zcbar  = 1;
beta   = 0.99;
delta  = 0.02;
rho_a  = 0.9;
rho_zx = 0.9;
rho_zc = 0.9;
alpha  = 0.35;
sigma  = 1.5;

// Exact equation system
model;

  // Euler equation
  zc*c^(-sigma)*(1/zx)
    = beta*(zc(+1)*c(+1)^(-sigma)*
           (a(+1)*alpha*k^(alpha-1) + (1-delta)/zx(+1)));

  // Resource constraint
  c + (k - (1-delta)*k(-1))/zx = a*k(-1)^alpha;

  // GDP
  y = a*k(-1)^alpha;

  // Investment
  inv = (k - (1-delta)*k(-1))/zx;

  // Wage
  w = (1-alpha)*y;

  //replacement price of capital
  repla = 1- ((1-delta)/zx);

  // Rental price (net)
  r = alpha*a*k(-1)^(alpha-1) - repla;

  // Stochastic processes
  zc - zcbar = rho_zc*(zc(-1) - zcbar) + eps_zc;
  zx - zxbar = rho_zx*(zx(-1) - zxbar) + eps_zx;
  a  - abar  = rho_a *(a(-1)  - abar)  + eps_a;

end;

steady_state_model;

  a  = abar;
  zc = zcbar;
  zx = zxbar;

  k   = ((1 - beta*(1-delta)) / (alpha*beta*a*zx))^(1/(alpha-1));
  c   = a*k^alpha - (delta/zx)*k; 
  y   = a*k^alpha;
  inv = (delta/zx)*k;
  w   = (1-alpha)*y;
repla = 1- ((1-delta)/zx);
  r   = alpha*a*k^(alpha-1) - repla;

end;

steady;
check;

shocks;
  var eps_a;  stderr 0.01;
  var eps_zx; stderr 0.01;
  var eps_zc; stderr 0.01;
end;

stoch_simul(irf = 100, order = 1) c k y inv w r repla a zx zc ;
