---
title: Fast Matrix-Vector Products
categories:
- Matrix Computation
---

## The Fast Fourier Transform
The ***discrete Fourier transform*** (DFT) of a vector \\(\color{black}{x \in C^n}\\) is a matrix-vector product
\\[\color{black}{y = F_nx}\\]
where the DFT matrix \\(\color{black}{F_n = (f_{ij}) \in C^{n\times n}}\\) is defined by
\\[\color{black}{f_{kj} = \omega^{(k-1)(j-1)}_n}\\]
with
\\[\color{black}{\omega_n = exp(-2\pi i/n) = cos(2\pi/n) - i\centerdot sin(2\pi/n).(\omega_n^n = 1)}\\]

The DFT is ubiquitous throughout computational science and engineering and one reason has to do with the following property:
> If \\(\color{black}{n}\\) is highly composite, then it is possible to carry out the DFT in many fewer than the \\(\color{black}{O(n^2)}\\) flops required by conventional matrix-vector multiplication.

If \\(\color{black}{n = 2m}\\), then \\(\color{black}{y = F_nx}\\) is given by
\\[\color{black}{y(1:m) = y_T+d.*y_B}\\]

\\[\color{black}{y(m+1:n) = y_T-d.*y_B}\\]
where \\(\color{black}{d = [1,\omega_n,...,\omega_n^{m-1}]^T}\\) and
\\[\color{black}{y_T = F_mx(1:2:n)}\\]

\\[\color{black}{y_B = F_mx(2:2:n)}\\]
For \\(\color{black}{n=2^t}\\), we can recur on this process until \\(\color{black}{n=1}\\), noting that \\(\color{black}{F_1x=x}\\)
```Matlab
function y = fft(x, n)
    if n = 1
        y = x
    else
        m = n/2
        y_T = fft(x(1:2:n), m)
        y_B = fft(x(2:2:n), m)
        w = exp(-2\pi i/n)
        d = [1,w,...,w^{m-1}]^T
        z = d.*y_B
        y = [y_T+z;y_T-z]
    end
```
