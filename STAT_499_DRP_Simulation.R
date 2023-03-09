install.packages("igraph")
install.packages("expm")
library(expm)
library(igraph)

setwd("~/STAT_390_Files")

# First set up a symmetric random walk

sym_walk <- function(x0, n)
{
  values <- c(0)
  current <- x0
  while (length(values) < n)
  {
    change = 0
    if (runif(1) > 0.5)
    {
      change = 1
    }
    else
    {
      change = -1
    }
    current <- current + change
    values <- c(values, current)
  }
  values
}

plot(sym_walk(0, 100), xlab ="Number of Steps (n)",
     ylab = "Sn", main = "Simple Symmetric Random
     Walk with Probability p = 0.5", type = "l")

# Second set up a random walk without symmetry

prob_walk <- function(x0, p, n)
{
  values <- c(0)
  current <- x0
  while (length(values) < n)
  {
    change = 0
    if (runif(1) > p)
    {
      change = 1
    }
    else
    {
      change = -1
    }
    current <- current + change
    values <- c(values, current)
    }
  values
}

plot(prob_walk(0, 0.25, 100), xlab =
    "Number of Steps (n)", ylab = "Sn", 
    main = "Simple Symmetric Random
     Walk with Probability p = 0.25", type = "l")

# Third graph random walk on graph w/ equal prob

G <- make_graph(~ 1 - 2, 3 - 4, 4 - 1:2, 5 - 1:2:3)
plot(G, vertex.color = "green")

# Degree Matrix of Graph G
deg <- degree(G)
G_deg <- diag(deg, 5, 5)

# Adjacency Matrix of Graph G
G_adj <- matrix(as_adj(G), byrow = FALSE, nrow = 5)

# Manually Computed Laplacian
G_man_lap <- G_deg - G_adj

# Computer Computed Laplacian
G_com_lap <- matrix(laplacian_matrix(G), 
                    byrow = FALSE, nrow = 5)

# Manually Computed Normalized Laplacian
deg_neg_half <- degree(G)^(-0.5)
G_deg_neg_half <- diag(deg_neg_half, 5, 5)
G_man_norm_lap <- G_deg_neg_half %*% G_man_lap %*%
  G_deg_neg_half

# Computer Computer Normalized Laplacian
G_com_norm_lap <- matrix(laplacian_matrix(G, 
                        normalized = TRUE), 
                         byrow = FALSE, nrow = 5)

# That proves the Laplacian matrix as well as our
# equation to normalize it.

# Now to do some stuff with the original M matrix
H <- make_graph(~ A - B:C, B - C:D)
plot(H, vertex.color = "light blue")

trans_prob_step <- function(graph, t)
{
  nodes = vcount(graph)
  d <- degree(graph)
  M <- matrix(0, nodes, nodes)
  for (i in 1:nodes)
  {
    for (j in 1:nodes)
    {
      if (are_adjacent(graph, i, j) == TRUE)
      {
        M[i, j] = 1/d[i]
      }
      else
      {
        M[i, j] = 0
      }
    }
  }
  if (t == 0)
  {
    print("Invalid Operation.")
  }
  else if (t == 1)
  {
    M <- M
  }
  else
  {
    for (k in 1:(t-1))
    {
      M <- M %*% M
    }
  }
  N <- M
}

N <- trans_prob_step(H, 2)

prob_cloud <- function(graph, t, i)
{
  A <- trans_prob_step(graph, t)
  v <- A[i,]
  v
}

row <- prob_cloud(H, 2, 1)

lazy_matrix_step <- function(graph, t)
{
  nodes = vcount(graph)
  d <- degree(graph)
  M <- matrix(0, nodes, nodes)
  if (runif(1) > 0.5)
  {
    for (i in 1:nodes)
    {
      for (j in 1:nodes)
      {
        if (are_adjacent(graph, i, j) == TRUE)
        {
          M[i, j] = 1/d[i]
        }
        else
        {
          M[i, j] = 0
        }
      }
    }
  }
  else
  {
    M <- diag(nodes)
  }
  if (t >= 2)
  {
    k = 1
    while (k < t)
    {
      S <- matrix(0, nodes, nodes)
      if (runif(1) > 0.5)
      {
        for (i in 1:nodes)
        {
          for (j in 1:nodes)
          {
            if (are_adjacent(graph, i, j) == TRUE)
            {
              S[i, j] = 1/d[i]
            }
            else
            {
              S[i, j] = 0
            }
          }
        }
      }
      else
      {
        S <- diag(nodes)
      }
      M = M %*% S
      k = k + 1
    }
  }
  A <- M
}

L <- lazy_matrix_step(H, 2)

I <- make_graph(~ a - b, c - a:b, d - a:b, e - a:c)
plot(I, vertex.color = "light green")

convergence <- function(graph)
{
  nodes = vcount(graph)
  val = 1
  R <- trans_prob_step(graph, val)
  while ((R[1, 1] + R[2, 1]) != (2 * R[1, 1]))
  {
    val = val + 1
    R <- trans_prob_step(graph, val)
  }
  val
}

num <- convergence(I)

converg_matrix <- function(graph)
{
  t <- convergence(graph)
  R <- trans_prob_step(graph, t)
  R
}

decomp_trans_prob <- function(graph, p)
{
  M <- trans_prob_step(graph, 1)
  nodes = vcount(graph)
  deg <- degree(graph)
  D <- diag(deg, nodes, nodes)
  half_D <- diag(degree(graph)^(0.5), nodes, nodes)
  neg_half_D <- diag(degree(graph)^(-0.5), 
                     nodes, nodes)
  S <- half_D %*% M %*% neg_half_D
  evs <- eigen(S)
  evals <- evs$values
  evects <- evs$vectors
  lambda <- diag(evals, nodes, nodes)
  phi <- neg_half_D %*% evects
  psi <- half_D %*% evects
  E <- matrix(0, nodes, nodes)
  for (a in 1:nodes)
  {
    r <- t(t(phi[,a]))
    y <- t(psi[,a])
    E <- E + ((evals[a]^p) * r %*% y)
  }
  E <- round(E, digits = 7)
}

E <- decomp_trans_prob(H, 2)
M <- trans_prob_step(H, 2)

J <- make_graph(~ 1 - 2:3, 2 - 3:4, 3 - 4:5, 4 - 5)
plot(J, vertex.color = "light green")

F <- converg_matrix(J)
power <- convergence(J)
P <- trans_prob_step(J, 9)