---
title: Frank-Wolfe Algorithm
categories:
- Optimization
---

## Introduction
The **Frank-Worlf Algorithm** is an iterative first-order optimization algorithm for constrained convex optimization. Also known as the conditional gradient method, reduced gradient algorithm and the convex combination algorithm, the method was originally proposed by **Marguerite Frank** and **Philip Wolfe** in 1956. In each iteration, the Frank-Wolfe algorithm considers a linear approximation of the objective function, and moves towards a minimizer of the this linear function(taken over the same domain).

## Problem Statement
Consider the constrained problem
\\[\color{black}{\min_x f(x) \text{ subject to } x \in C}\\]

where \\(\color{black}{f}\\) is convex and smooth, and \\(\color{black}{C}\\) is convex. 

## Conditional gradient method
Also known as the Frank-Wolfe method, uses a local linear expansion of \\(\color{black}{f}\\):
\\[\color{black}{s^{(k-1)}\in \min_{s\in C}\triangledown f(x^{(k-1)})^Ts}\\]

\\[\color{black}{x^{(k)}=(1-\gamma_k)x^{(k-1)}+\gamma_ks^{(k-1)}}\\]

Note that there is no projection; update is solved directly over the constraint set \\(\color{black}{C}\\)

The default choice for step sizes is \\(\color{black}{\gamma_k = \frac{2}{(k+1)}, k = 1,2,3,...}\\) For any choice \\(\color{black}{0 \leq \gamma_k \leq 1,}\\) we see that \\(\color{black}{x^{(k)} \in C}\\) by convexity. Can also think of the update as
\\[\color{black}{x^{(k)} = x^{(k-1)} + \gamma_k(s^{(k-1)}-x^{(k-1)})}\\]

i.e., we are moving less and less in the direction of the linearization minimizer as the algorithm proceeds.

![A step of the Frank-Wolfe algorithm](https://upload.wikimedia.org/wikipedia/commons/thumb/e/e5/Frank-Wolfe_Algorithm.png/1024px-Frank-Wolfe_Algorithm.png)

## Norm constraints
When \\(\color{black}{C = \lbrace x:\|\|x\|\| \leq t\rbrace}\\) for a norm \\(\color{black}{\|\|\centerdot\|\|}\\). Then

\\[\color{black}{s \in \min_{\|\|s\|\|\leq t} \triangledown f(x^{(k-1)})^Ts}\\]

\\[\color{black}{= -t\centerdot (\max_{\|\|s\|\|\leq 1}\triangledown f(x^{(k-1)})^Ts)}\\]

\\[\color{black}{= -t\centerdot \partial \|\|\triangledown f(x^{(k-1)})\|\|_*}\\]

where \\(\color{black}{\|\|\centerdot\|\|_*}\\) is the corresponding dual norm. In other words, if we know how to compute **subgradients of the dual norm**, the we can easily perform Frank-Wolfe steps.

A key to Frank-Wolfe: this can often be simpler or cheaper than projection onto \\(\color{black}{C = \lbrace x: \|\|x\|\| \leq t \rbrace}\\). Also often simpler or cheaper than the proximal operator for \\(\color{black}{\|\|\centerdot\|\|}\\)

## Properties
- While competing methods such as gradient descent for constrained optimization require a projection step back to the feasible set in each iteration, the Frank–Wolfe algorithm only needs the solution of a linear problem over the same set in each iteration, and automatically stays in the feasible set.

- The convergence of the Frank–Wolfe algorithm is sublinear in general: the error in the objective function to the optimum is \\(\color{black}{O( 1 / k )}\\) after \\(\color{black}{k}\\) iterations, so long as the gradient is Lipschitz continuous with respect to some norm. The same convergence rate can also be shown if the sub-problems are only solved approximately.

- The iterates of the algorithm can always be represented as a sparse convex combination of the extreme points of the feasible set, which has helped to the popularity of the algorithm for sparse greedy optimization in machine learning and signal processing problems, as well as for example the optimization of minimum–cost flows in transportation networks.

- If the feasible set is given by a set of linear constraints, then the subproblem to be solved in each iteration becomes a linear program.

- While the worst-case convergence rate with \\(\color{black}{O ( 1 / k )}\\) can not be improved in general, faster convergence can be obtained for special problem classes, such as some strongly convex problems.