---
title: Symmetric Storage Gaxpy
categories:
- Matrix Computation
---

## Introduction
A matrix \\(\color{black}{A\in\mathcal{R}^{n\times n}}\\) is *symmetric* if \\(\color{black}{A^{T}=A}\\) and *skew-symmetric* if \\(\color{black}{A^{T}=-A}\\). Likewise, a matrix \\(\color{black}{A\in\mathcal{C}^{n\times n}}\\) is *hermitian* if \\(\color{black}{A^{H}=A}\\) and *skew-hermitian* if \\(\color{black}{A^{H}=-A}\\).

For such matrices, storage requirements can be halved by simply storing the lower triangle of elements.

For general \\(\color{black}{n}\\), we set

\\[\color{black}{A.vec((n-j/2)(j-1)+i)=a_{ij} 1\le j\le i \le n}\\]

```Matlab
for j = 1:n
    for i = 1:j-1
        y(i) = y(i) + A.vec((i-1)n - i(i-1)/2 + j)x(j)
    end
    for i = j:n
        y(i) = y(i) + A.vec((j-1)n - j(j-1)/2 + i)x(j)
    end
end
```

This algorithm requires the same \\(\color{black}{2n^2}\\) flops that an ordinary gaxpy requires.