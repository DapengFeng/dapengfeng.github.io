---
title: Band Storage Gaxpy
categories:
- Matrix Computation
---

## Introduction
Suppose \\(\color{black}{A\in\mathcal{R}^{n\times n}}\\) has lower bandwith \\(\color{black}{p}\\) and upper banwidth \\(\color{black}{q}\\) and assume that \\(\color{black}{p}\\) and \\(\color{black}{q}\\) are much smaller than \\(\color{black}{n}\\). Such a matrix can be stored in a \\(\color{black}{(p+q+1)-\text{ by }-n}\\) array \\(\color{black}{A.band}\\) with the convention that

\\[\color{black}{a_{ij}=A.band(i-j+q+1,j)}\\]

for all \\(\color{black}{(i,j)}\\) that fall inside the band.

```Matlab
for j = 1:n
    alpha_1 =  max(1, j-q), alpha_2 = min(n, j+p)
    beta_1 = max(1, q+2-j), beta_2 = beta_1 + alpha_2 - alpha_1
    y(alpha_1:alpha_2) + A.band(beta_1:beta_2,j)x(j)
end
```

The algorithm involves just \\(\color{black}{2n(p+q+1)}\\) flops with the assumption that \\(\color{black}{p}\\) and \\(\color{black}{q}\\) are much smaller than \\(\color{black}{n}\\).
