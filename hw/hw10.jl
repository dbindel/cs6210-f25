### A Pluto.jl notebook ###
# v0.20.21

using Markdown
using InteractiveUtils

# ╔═╡ ca0206b9-e4ee-4f5f-8685-80efaf7ff0ba
using LinearAlgebra

# ╔═╡ 97e46e8a-227b-11ed-3ad4-bf493f45d214
md"""
# HW10 for CS 6210

You may (and probably should) talk about problems with the each other, with the TA, and with me, providing attribution for any good ideas you might get. Your final write-up should be your own.
"""

# ╔═╡ 22f06cec-14ac-434b-865a-8f360f09fb62
md"""
## 1. Your turn

Submit a course evaluation!  I can see if you have submitted an evaluation (even if I don't see what it is).  This is worth 2 points on the final grade.
"""

# ╔═╡ 87d3be02-41fe-4222-ae6c-ae207e5649aa
md"""
## 2. Scaled iteration

Suppose $A = QR$ is an economy QR decomposition for $A \in \mathbb{R}^{m \times n}$ a tall thin matrix of full column rank.  Consider the stationary iteration

$$x^{(k+1)} = R^{-1} (Q^T b - \lambda R^{-T} x^{(k)})$$

Suppose $x^{(0)} = 0$.  Then write the other iterates $x^{(k)}$ in the form $x^{(k)} = V f_k(\Sigma) U^T b$ where $A = U \Sigma V^T$ is the economy SVD.  For what range of $\lambda$ does the iteration converge?
"""

# ╔═╡ 7e3dfbbd-690f-43d5-8e17-7e77e13a4168
md"""
## 3. Rational roots

Consider a rational function of the form

$$f(z) = c^T (zI-A)^{-1} b + d$$

Write a routine `rateval(z, F, b, c, d)` to compute $f(z)$ in $O(n^2)$ time given the precomputed Hessenberg factorization from `F = hessenberg(A)`.

Note that we can write this by as $f(z) = c^T u$ where $u = (zI-A)^{-1} b$, or equivalently

$$(A-zI) u + b = 0.$$

Putting it all together, we have

$$\begin{bmatrix} A-zI & b \\ c^T & d \end{bmatrix} \begin{bmatrix} u \\ 1 \end{bmatrix}
 = \begin{bmatrix} 0 \\ f(z) \end{bmatrix}.$$

The function call `eigvals(A,B)` finds all solutions of the eigenvalue problem $Ax = \lambda Bx$.  Using this function, write a routine `ratroots(A, b, c, d)` to compute all the *real* zeros of $f$.

You should check correctness of the roots returned by `ratroots` using the fast evaluation code from `rateval`.
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
# ╟─22f06cec-14ac-434b-865a-8f360f09fb62
# ╟─87d3be02-41fe-4222-ae6c-ae207e5649aa
# ╟─7e3dfbbd-690f-43d5-8e17-7e77e13a4168
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
