### A Pluto.jl notebook ###
# v0.20.19

using Markdown
using InteractiveUtils

# ╔═╡ ca0206b9-e4ee-4f5f-8685-80efaf7ff0ba
using LinearAlgebra

# ╔═╡ 97e46e8a-227b-11ed-3ad4-bf493f45d214
md"""
# HW5 for CS 6210

You may (and probably should) talk about problems with the each other, with the TA, and with me, providing attribution for any good ideas you might get. Your final write-up should be your own.
"""

# ╔═╡ d3eee2ee-57bc-414b-8c6e-ccd0c6a6b8f5
md"""
## 1. Collapsing columns

Suppose $A = QR$ is an economy QR factorization with a dense representation of the (tall thin) $Q$ factor.  Now consider the QR factorization of $AP = QRP$, where the action of $P$ is to remove column $j$ of $A$ (or $R$).  Via a sequence of Givens rotations, compute $AP = \tilde{Q} \tilde{R}$ where $\tilde{R}$ is the first $n-1$ rows of $G_1 G_2 \ldots G_k (RP)$ (chosen so $\tilde{R}$ is again upper triangular) and $\tilde{Q}$ is the first $n-1$ columns of $Q G_k^T G_{k-1}^T \ldots G_1^T$.  You may use the [`LinearAlgebra.givens`](https://docs.julialang.org/en/v1/stdlib/LinearAlgebra/#LinearAlgebra.Givens) functions in Julia to compute the appropriate rotations; you will want to use `lmul!` and `rmul!` to apply the rotations in place efficiently.  Your code should take $O(mn)$ time.
"""

# ╔═╡ bf5cb221-858e-487b-b619-bb7c51b5a373
function delcol_qr_naive(Q, R, j)
	m, n = size(Q)
	RP = R[:,1:n .!= j]
	F = qr(Q*RP)
	F.Q[:,1:n-1], F.R
end

# ╔═╡ ccb3f3ca-b50f-4b6c-acc3-b42504a5bc1f
let
	# Compute QR factorization of the original matrix
	A = rand(10,5)
	F = qr(A)
	Q = F.Q[:,1:5]
	R = F.R

	# Compute QR factorization of the permuted matrix
	j = 2
	AP = A[:,1:5 .!= j]
	Qt, Rt = delcol_qr_naive(Q, R, j)  # Replace with delcol_qr (your code)

	# Check that Qt^T Qt = I, AP = Qt*Rt, and Rt is upper triangular
	norm(Qt'*Qt-I), norm(Qt*Rt-AP)/norm(AP), norm(tril(Rt,-1))
end

# ╔═╡ 2dbf25c0-36cd-4d5c-8ca8-17184cd67a68
md"""
## 2. Regularized residuals

Define the Tikhonov regularized least squares solution operator as

$F(\lambda) = (A^T A + \lambda^2 I)^{-1} A^T b$

and let $r(\lambda) = (I - AF(\lambda)) b$ be the corresponding residual.  Show that

$\|r(\lambda)\|^2 = \|r(0)\|^2 + \|z(\lambda)\|^2$

where

$z(\lambda) = r(\lambda)-r(0) = \lambda^2 F(\lambda)^T x(0).$
"""

# ╔═╡ 6d15dc29-8ae2-424f-8c31-f8a66ffd3d74
md"""
## 3. A nearness problem

For a given matrix $A \in \mathbb{R}^{n \times n}$, write a short code

    nearest_sum_to_one(A)

that computes the closest matrix to $A$ whose entries sum to one.  Also provide a tester to sanity check correctness of your code -- check that the constraint is satisfied, and that a random directional derivative in a direction consistent with the constraints (i.e. a direction whose entries sum to zero) is zero.
"""

# ╔═╡ 00000000-0000-0000-0000-000000000001
PLUTO_PROJECT_TOML_CONTENTS = """
[deps]
LinearAlgebra = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"
"""

# ╔═╡ 00000000-0000-0000-0000-000000000002
PLUTO_MANIFEST_TOML_CONTENTS = """
# This file is machine-generated - editing it directly is not advised

julia_version = "1.12.0"
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
version = "5.13.1+1"
"""

# ╔═╡ Cell order:
# ╟─97e46e8a-227b-11ed-3ad4-bf493f45d214
# ╠═ca0206b9-e4ee-4f5f-8685-80efaf7ff0ba
# ╟─d3eee2ee-57bc-414b-8c6e-ccd0c6a6b8f5
# ╠═ccb3f3ca-b50f-4b6c-acc3-b42504a5bc1f
# ╠═bf5cb221-858e-487b-b619-bb7c51b5a373
# ╟─2dbf25c0-36cd-4d5c-8ca8-17184cd67a68
# ╟─6d15dc29-8ae2-424f-8c31-f8a66ffd3d74
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
