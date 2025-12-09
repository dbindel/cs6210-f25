### A Pluto.jl notebook ###
# v0.20.21

using Markdown
using InteractiveUtils

# ╔═╡ ca0206b9-e4ee-4f5f-8685-80efaf7ff0ba
using LinearAlgebra

# ╔═╡ 97e46e8a-227b-11ed-3ad4-bf493f45d214
md"""
# Final CS 6210

You may talk about problems with the TA and with me, or reference any relevant literature (including the Julia documentation!), providing attribution for any good ideas you might get. You should not consult with other people or with AI tools, and your final write-up should be your own.
"""

# ╔═╡ dd7e7d53-24ff-4c27-8344-e8ac81cdad1b
md"""
## 1. Truth or dare

State whether each of the following is true and give a short argument why or why not.
(One point each)

- In the Tikhonov-regularized least squares problem, choosing the regularization parameter $\lambda$ via leave-one-out cross validation minimizes the overall residual norm.
- If $A \in \mathbb{R}^{n \times n}$ is symmetric positive definite and $A = R^T R$ is the Cholesky factorization, then the column two-norms of $R$ satisfy $\lambda_n(A) \leq \|R_{:,j}\|^2 \leq \lambda_1(A)$, where $\lambda_1(A)$ and $\lambda_n(A)$ are the largest and smallest eigenvalues of $A$, respectively.
- Suppose $A \in \mathbb{R}^{n \times n}$ is written as $A = M-N$ where $M$ is invertible.  If the iteration $M x^{(k+1)} = N x^{(k)} + b$ converges to a limiting fixed point $x^*$, then $Ax^* = b$.
- If $A$ is diagonalizable and $p$ is a polynomial such that $p(\lambda) = \lambda^{-1}$ for every eigenvalue $\lambda$ in the spectrum of $A$, then $p(A) = A^{-1}$.
"""

# ╔═╡ 87d3be02-41fe-4222-ae6c-ae207e5649aa
md"""
## 2. Scaled iteration

Consider the stationary iteration

$$x^{(k+1)} = x^{(k)} + M^{-1} A^T (b-Ax^{(k)})$$

where $A \in \mathbb{R}^{m \times n}$ has full column rank and $M \in \mathbb{R}^{n \times n}$ is positive definite.

- (1 point): What is the fixed point?
- (1 point): Write an iteration equation for the error
- (2 points): Argue that if $\|Ax\|^2/\|x\|_M^2 < 2$ for all $x \neq 0$, then the iteration converges.
"""

# ╔═╡ 01afcd57-e86a-45ef-b86c-0a49223417ea
md"""
## 3. Factorization free

Given an initial guess $x^{(0)}$ with $|x^{(0)}-\lambda^{-1}| < |\lambda^{-1}|$
the Newton iteration

$$x^{(k+1)} = x^{(k)} + x^{(k)}(1-\lambda x^{(k)})$$

converges quadratically to $\lambda^{-1}$; more specifically, the error $e^{(k)} = x^{(k)} - \lambda^{-1}$ satisfies

$$e^{(k+1)} = -\lambda (e^{(k)})^2.$$

Similarly, the matrix iteration

$$X^{(k+1)} = X^{(k)} + X^{(k)} (I - A X^{(k)})$$

converges quadratically to $A^{-1}$ given a good initial guess.

(1 point): Argue that for the initial guess $x^{(0)} = 2-\lambda$, the scalar Newton iteration converges for $|1-\lambda| < 1$.  You may take the convergence condition given at the start of the problem as true without proof.

(1 point): Argue that if $A$ is diagonalizable and $X^{(0)} = p_0(A)$ for some polynomial $p_0$, then the matrix Newton iteration converges iff the scalar Newton iteration converges for each eigenvalue $\lambda$ of $A$ starting from $p_0(\lambda)$.

(2 points): Argue that if $\|I-A\| < 1$ for some operator norm, then the matrix Newton iteration converges starting from the initial guess $X^{(0)} = 2I-A$
"""

# ╔═╡ 7e3dfbbd-690f-43d5-8e17-7e77e13a4168
md"""
## 4. Rational extrema

Consider a rational function of the form

$$f(z) = c^T (zI-A)^{-1} b + d$$

Note that this can also be written via a Schur complement in the matrix

$$\begin{bmatrix} A-zI & b \\ c^T & d \end{bmatrix} = M-zD$$

The poles of $f$ are the eigenvalues of $A$ and the zeros are the eigenvalues of $M-zD$ (as discussed in HW 10).

(1 point): Argue that if $d \neq 0$, then the zeros of $f$ are the eigenvalues of $A-bc^T/d$.

(1 point): Write a function to evaluate the function and derivative in $O(n^2)$ time given a Hessenberg factorization object `F = hessenberg(A)`.  Your code should take $O(n^2)$ time, and should look like

    f, df = p4eval(z, F, b, c, d)

(1 point): Introducing variables $v = (zI-A)^{-1} b$ and $u = (zI-A)^{-1} v$, write a generalized eigenvalue problem whose solutions correspond to critical points of $f(z)$.

(1 point): Using `eigvals(A,B)` (which computes solutions to the generalized eigenvalue problem $Ax = \lambda Bx$), find the (real) critical points of $f$.  Your code should look like

    zs = p4critical(A, b, c, d)

Note that you should write testers for each code.
"""

# ╔═╡ ca64c90c-09da-4a34-9e81-c5ce04c37bf6
md"""
## 5. Playing with polynomials

Suppose $A = I - F$ where $\|F\|_2 < 1$.

(2 points): Find $p \in \mathcal{P}_d$ such that 
$\|I-A p(A)\|_2 \leq \|F\|_2^{d+1}$.

(2 points): Given a particular $F$, write a function `polyinv(F,d)` to compute the $d+1$ coefficients with respect to the power basis for the polynomial $p \in \mathcal{P}_d$ such that $\|I-A p(A)\|_F^2$ is minimal.  Also give an explanation for your code!
"""

# ╔═╡ 00000000-0000-0000-0000-000000000001
PLUTO_PROJECT_TOML_CONTENTS = """
[deps]
LinearAlgebra = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"
"""

# ╔═╡ 00000000-0000-0000-0000-000000000002
PLUTO_MANIFEST_TOML_CONTENTS = """
# This file is machine-generated - editing it directly is not advised

julia_version = "1.12.2"
manifest_format = "2.0"
project_hash = "f352ceee806168c8ae38887a01d7bae6ca62470b"

[[deps.Artifacts]]
uuid = "56f22d72-fd6d-98f1-02f0-08ddc0907c33"
version = "1.11.0"

[[deps.CompilerSupportLibraries_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "e66e0078-7015-5450-92f7-15fbd957f2ae"
version = "1.3.0+1"

[[deps.Libdl]]
uuid = "8f399da3-3557-5675-b5ff-fb832c97cbdb"
version = "1.11.0"

[[deps.LinearAlgebra]]
deps = ["Libdl", "OpenBLAS_jll", "libblastrampoline_jll"]
uuid = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"
version = "1.12.0"

[[deps.OpenBLAS_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Libdl"]
uuid = "4536629a-c528-5b80-bd46-f80d51c5b363"
version = "0.3.29+0"

[[deps.libblastrampoline_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "8e850b90-86db-534c-a0d3-1478176c7d93"
version = "5.15.0+0"
"""

# ╔═╡ Cell order:
# ╟─97e46e8a-227b-11ed-3ad4-bf493f45d214
# ╠═ca0206b9-e4ee-4f5f-8685-80efaf7ff0ba
# ╟─dd7e7d53-24ff-4c27-8344-e8ac81cdad1b
# ╟─87d3be02-41fe-4222-ae6c-ae207e5649aa
# ╟─01afcd57-e86a-45ef-b86c-0a49223417ea
# ╟─7e3dfbbd-690f-43d5-8e17-7e77e13a4168
# ╟─ca64c90c-09da-4a34-9e81-c5ce04c37bf6
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
