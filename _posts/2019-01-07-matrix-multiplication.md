---
title: Matrix Multiplication
categories:
- Matrix Computation
---

## The Notion of "Level" and the BLAS
The dot product and single alpha x plus y(saxpy) operations are example of *level-1* operations. Level-1 operations involve an amount of data and an amount of arithmetic that area linear in the dimension of the operation. An \\(\color{black}{m}\\)-by-\\(\color{black}{n}\\) outer product update or a generalized saxpy(gaxpy) operation involves a quadratic amount of data \\(\color{black}{O(mn)}\\) and a quadratic amount of work \\(\color{black}{O(mn)}\\). These are *level-2* operations. The matrix multiplication update \\(\color{black}{C=C+AB}\\) is a *level-3* operation. Level-3 operations are quadratic in data and cubic in work.

Important level-1, level-2, and level-3 operations are encapsulated in the "BLAS," an acronym that stands for <u>B</u>asic <u>L</u>inear <u>A</u>lgebra <u>S</u>ubprograms. See [OpenBLAS](http://www.openblas.net/). The design of matrix algorithm that are rich in level-3 BLAS operations is a major preoccupation of the field for reasons that have to do with data reuse.

## Structure and Efficiency
The efficiency of a given matrix algorithm depends upon several factors. Most obvious and what we treat in this section is the amount of required arithmetic and storage. How to reason about these important attributes is nicely illustrated by considering examples that involve triangular matrices, diagonal matrices, banded matrices, symmetric matrices, and permutation matrices. These are among the most important types of structured matrices that arise in practice, and various economies can be realized if they are involved in a calculation.

## Vectorization and Locality
When it comes to designing a high-performance matrix computation, it is not enough simply to minimize flops. Attention must be paid to how the arithmetic units interact with the underlying memory system. Data structures are an important part of the picture because not all matrix layouts are "architecture friendly." 
### Vector Processing
An individual floating point operation typically requires several cycles to complete. Vector processors exploit the fact that a vector operation is a very regular sequence of scalar operations. The key idea is pipelining. With pipelining, a vector processor comes with a repertoire of *vector instructions*, such as vector add, vector multiply, vector scale, dot product, and saxpy. These operations take place in vector registers with input and output handled by *vector load* and *vector store* instructions. An important attribute of a vector processor is the length \\(\color{black}{v_L}\\) of the vector registers that carry out the vector operations. A length-n vector operations of length \\(\color{black}{v_L}\\) or less. Here is how such a partitioning might be managed for a vector addition \\(\color{black}{z=x+y}\\) where \\(\color{black}{x}\\) and \\(\color{black}{y}\\) are n-vectors:
```Matlab
first = 1
while first <= n
    last = min{n, first+v_L-1}
    r_1 = x(first:last) %Vector load
    r_2 = y(first:last) %Vector load
    r_1 = r_1 + r_2 %Vector add
    z(first:last) = r_1 %Vector store
    first = last + 1
end
```

The vector addition is a register-register operation while the "flopless" movement of data to and from the vector registers. For clarity, assume that \\(\color{black}{n}\\) is very large and an integral multiple of \\(\color{black}{v_L}\\), thereby making it safe to ignore the final cleanup pass through the loop.

Regarding the vectorized addition \\(\color{black}{r_1=r_1+r_2}\\), assume it takes \\(\color{black}{\tau_{add}}\\) cycles to fill the pipeline and that once this happens, a component of \\(\color{black}{z}\\) is produced each cycle. It follows that
\\[\color{black}{N_{arith} = (\frac{n}{v_L})(\tau_{add}+v_L) = (\frac{\tau_{add}}{v_L}+1)n}\\]
accounts for the total number cycles that requires for arithmetic.

For the vector loads and stores, assume that \\(\color{black}{\tau_{data}+v_L}\\) cycles are required to transport a length-\\(\color{black}{v_L}\\) vector from memory to a register or from a register to memory, where \\(\color{black}{\tau_{data}}\\) is the number of cycles required to fill the data pipeline. With these assumptions we see that
\\[\color{black}{N_{data} = 3(\frac{n}{v_L})(\tau_{data}+v_L) = 3(\frac{\tau_{data}}{v_L}+1)n}\\]
specifies the number of cycles that are required to get data to and from the registers.

The arithmetic-to-data-motion ratio
\\[\color{black}{N_{arith}/N_{data} = \frac{\tau_{add}+v_L}{3(\tau_{add}+v_L)}}\\]
and the total cycles sum
\\[\color{black}{N_{arith}+N_{data} = (\frac{\tau_{arith}+3\tau_{data}}{v_L}+4)n}\\]
are illuminating statistics, but they are not necessarily good predictions of performance. 

### Gaxpy versus Outer Product
Two algorithms that involve the same number of flops number of flops can have substantially different data motion properties. Consider the \\(\color{black}{n}\\)-by-\\(\color{black}{n}\\) gaxpy
\\[\color{black}{y = y + Ax}\\]
and the \\(\color{black}{n}\\)-by-\\(\color{black}{n}\\) outer product update
\\[\color{black}{A = A + yx^T}\\]

Both of these level-2 operations involve \\(\color{black}{2n^2}\\) flops. However, if assume that \\(\color{black}{n=v_L}\\), then we can see that the gaxpy computation
```Matlab
r_x = x
r_y = y
for j = 1:n
    r_a = A(:,j)
    r_y = r_y + r_a * r_x(j)
end
y = r_y
``` 
requires \\(\color{black}{(3+n)}\\) load/store operations while for the outer product update
```Matlab
r_x = x
r_y = y
for j = 1:n
    r_a = A(:,j)
    r_a = r_a + r_y * r_x(j)
    A(:,j) = r_a
end
```
the corresponding count is \\(\color{black}{(2+2n)}\\). Thus, the data motion overhead for the outer product update is worse by a factor of 2, a reality that could that could be a factor in the design of a high-performance matrix computation.

### The Relevance of Stride
The time is takes to load a vector into a vector register may depend greatly on how the vector is laid out in memory. Two concepts will help frame the issue. A vector is said to have *unit stride* if its components are contiguous in memory. A matrix is said to be stored in *column-major order* if its columns have unit stride.

Let us consider the matrix multiplication update calculation
\\[\color{black}{C = C + AB}\\]
where it is assumed that the matrices \\(\color{black}{C \in R^{m \times n}, A \in R^{m \times r}}\\) and \\(\color{black}{B \in R^{r \times n}}\\) are stored in column-major order. Suppose the loading of a non-unit-stride vector. If so, then the implementation
```Matlab
for j=1:n
    for k=1:r
        C(:,j) = C(:,j) + A(:,k)*B(k,j)
    end
end
```
which accesses \\(\color{black}{C,A,}\\) and \\(\color{black}{B}\\) by column would be preferred to
```Matlab
for i=1:m
    for j=1:n
        C(i,j) = C(i,j) + A(i,:)*B(:,j)
    end
end
```
which accesses \\(\color{black}{C}\\) and \\(\color{black}{A}\\) by row. While this example points to the possible importance of stride, it is important to keep in the mind that the penalty for non-unit-stride access varies from system to system and may depend upon the value of the stride itself.

### Blocking for Data Reuse
Matrices reside in memory but memory has levels. The **cache** is a relatively small high-speed memory unit that sits just below the functional units where the arithmetic is carried out. During a matrix computation, matrix elements move up and down the memory hierarchy. The cache, which is a small heigh-speed memory situated in between the functional units and main memory, plays a particularity critical role. The overall design of the hierarchy varies from system to system. However, two maxims always apply:
- Each level in the hierarchy has a limited capacity and for economic reasons this capacity usually becomes smaller as we ascend the hierarchy.
- There is cost, sometimes relatively great, associated with the moving of data between two levels in the hierarchy.

The efficient implementation of a matrix algorithm requires an ability to reason about the flow of data between the various levels of storage.

To develop an appreciation for cache utilization we can consider the update \\(\color{black}{C=C+AB}\\) where each matrix is \\(\color{black}{n}\\)-by-\\(\color{black}{n}\\) and blocked in the same way.

Assume that these three matrices reside in main memory and that we plan to update \\(\color{black}{C}\\) block by block:
\\[\color{black}{C_{ij} = C_{ij} + \Sigma_{k=1}^pA_{ik}B{kj}}\\]

The data in the blocks must be brought up to the functional units via the cache which we assume is large enough to hold a C-block, an A-block, and a B-block. This enables to structure the computation as follows:
```Matlab
for i=1:q
    for j=1:r
        % Load C_ij from main memory into cache
        for k=1:p
            % Load A_ik from main memory into cache
            % Load B_kj from main memory into cache
            C_ij = C_ij + A_ik * B_kj
        end
        % Store C_ij in main memory
    end
end
```

Because blocking affects the amount of memory traffic in a matrix computation, it is of paramount importance when designing a high-performance implementation. Data structures are also important; storing a matrix by block rather than in column-major order could enhance performance.