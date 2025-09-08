### A Pluto.jl notebook ###
# v0.20.17

using Markdown
using InteractiveUtils

# ╔═╡ ca0206b9-e4ee-4f5f-8685-80efaf7ff0ba
using LinearAlgebra

# ╔═╡ 97e46e8a-227b-11ed-3ad4-bf493f45d214
md"""
# HW2 for CS 6210

You may (and probably should) talk about problems with the each other, with the TA, and with me, providing attribution for any good ideas you might get. Your final write-up should be your own.
"""

# ╔═╡ 5eae7aea-ba30-45cf-98c7-1e3f77d98cef
md"""
## 1. Rewrite for range

Consider the following function for $0 < x < 1$.  Write a function `p1mine` that revises `p1ref` to maintain accuracy for $x \ll 1$.  For part 2, you may want to look up the `log1p` function in Julia.  You may check your answers using arbiitrary precision arithmetic (see `BigFloat` in Julia), but your solutions should use only ordinary floating point.
"""

# ╔═╡ 978cc177-9a87-438f-adc7-1d7104f25ce7
# Reference formulas (inaccurate for x << 1)
function p1ref(x)
	y1 = sqrt(1+x)-sqrt(1-x)
	y2 = log(sqrt(1+1/x))-log(sqrt(1/x))
	y3 = cos(x)-1
	y1, y2, y3
end

# ╔═╡ 67ee9ba1-5fe2-40bc-81ba-ff2f94c7054c
md"""
## 2. Working backward

Consider the 2-by-2 rotation matrix

$$Q = \begin{bmatrix} c & -s \\ s & c \end{bmatrix}$$

where $c^2 + s^2 = 1$.  Note that $Q^T Q = I$.  From class, we know that the code `w = Q*v` computes $\hat{w} = (Q+E) v$ where $|E| \lesssim 2 \epsilon_{\mathrm{mach}} |Q|$.

- Argue that $\hat{w} = Q\hat{v}$ where $\hat{v}-v = Q^T E v$
- Argue that therefore $\hat{v}$ has a normwise relative error (in the two norm) bounded by $4 \epsilon_{\mathrm{mach}}$
"""

# ╔═╡ 9fecc880-fb4d-4a3a-ba1d-2ff81f90f737
md"""
## 3. Connect the dots

Give an example of a dot product that is computed to low relative accuracy when using the standard algorithm in Julia.
"""

# ╔═╡ f3e1a138-c3ce-433b-8eb6-c4430ad5c471
md"""
## 4. Back it up

Running the recurrence $E_n = 1-nE_n$ forward is an unstable way to compute 
$\int_0^1 x^n e^{x-1} \, dx$.  However, we can get good results by running the recurrence backward from the estimate $E_n \approx 1/(n+1)$ starting at large enough $N$.

- Implement this iteration starting from $E_{100} \approx 1/101$ to compute $E_{20}$.  You should get good agreement with the code given (which computes the numbers in a completely different way).
- Let $\hat{E}_n$ be the floating point value computed by the backward recurrence and $f_n = \hat{E}_n-E_n$ be the error.  Write a recurrence for $f_{n-1}$, which should take the form $f_{n-1} = (f_n + \gamma_n)/n$, where $\gamma_n$ captures the roundoff.
- Using a uniform estimate of $\gamma_n$, augment your code to compute a running (estimated) error bound for the $E_n$ computed by working backward from $E_n$.  Your bound should incorporate both roundoff and the starting error.  For controlling the starting error, you may use that $(n+1)^{-1} \leq E_n \leq (n+2)^{-1}$. 
- How large should the starting point $N$ be to compute $E_{20}$ to near machine precision?

NB: `eps(Float64)` gives the distance between `1.0` and the next largest floating point number.  This is twice the quantity that we refer to as $\epsilon_{\mathrm{mach}}$ in class.
"""

# ╔═╡ b5d0294a-91d9-4e27-be95-7e3556b8023a
# Compute En by integrating a truncated Taylor series term by term
function En_ref(n)
	s = 0.0
	for j = 20:-1:1
		s = (s+1.0/(n+j+1))/j
	end
	s += 1.0/(n+1)
	s/exp(1.0)
end

# ╔═╡ 00000000-0000-0000-0000-000000000001
PLUTO_PROJECT_TOML_CONTENTS = """
[deps]
LinearAlgebra = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"
"""

# ╔═╡ 00000000-0000-0000-0000-000000000002
PLUTO_MANIFEST_TOML_CONTENTS = """
# This file is machine-generated - editing it directly is not advised

julia_version = "1.11.6"
manifest_format = "2.0"
project_hash = "ac1187e548c6ab173ac57d4e72da1620216bce54"

[[deps.Artifacts]]
uuid = "56f22d72-fd6d-98f1-02f0-08ddc0907c33"
version = "1.11.0"

[[deps.CompilerSupportLibraries_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "e66e0078-7015-5450-92f7-15fbd957f2ae"
version = "1.1.1+0"

[[deps.Libdl]]
uuid = "8f399da3-3557-5675-b5ff-fb832c97cbdb"
version = "1.11.0"

[[deps.LinearAlgebra]]
deps = ["Libdl", "OpenBLAS_jll", "libblastrampoline_jll"]
uuid = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"
version = "1.11.0"

[[deps.OpenBLAS_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Libdl"]
uuid = "4536629a-c528-5b80-bd46-f80d51c5b363"
version = "0.3.27+1"

[[deps.libblastrampoline_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "8e850b90-86db-534c-a0d3-1478176c7d93"
version = "5.11.0+0"
"""

# ╔═╡ Cell order:
# ╟─97e46e8a-227b-11ed-3ad4-bf493f45d214
# ╠═ca0206b9-e4ee-4f5f-8685-80efaf7ff0ba
# ╟─5eae7aea-ba30-45cf-98c7-1e3f77d98cef
# ╠═978cc177-9a87-438f-adc7-1d7104f25ce7
# ╟─67ee9ba1-5fe2-40bc-81ba-ff2f94c7054c
# ╟─9fecc880-fb4d-4a3a-ba1d-2ff81f90f737
# ╟─f3e1a138-c3ce-433b-8eb6-c4430ad5c471
# ╠═b5d0294a-91d9-4e27-be95-7e3556b8023a
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
