n = 100
p = 3
t = 10
nz = 300

x = buildStData(X = array(rnorm(n = n*p*t), dim = c(n,p,t)), 
                Y = matrix(rnorm(n = n*t), nrow = n), 
                Z = matrix(rnorm(n = nz*t), nrow = nz), 
                tLabs = 1:t, 
                coords.r = matrix(rnorm(n = nz*2), ncol = 2),
                coords.s = matrix(rnorm(n = n*2), ncol = 2), 
                Q = NULL, 
                X.lab = 'xlab', Y.lab = 'ylab', Z.lab = 'zlab')