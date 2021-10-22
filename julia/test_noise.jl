### A Pluto.jl notebook ###
# v0.14.8

using Markdown
using InteractiveUtils

# ╔═╡ 9662f1cc-30bb-11ec-3315-db727be10430
using Yao, YaoExtensions,PastaQ, YaoPlots

# ╔═╡ 75cb8628-92f2-48d0-b344-caf72f355524
include("YaoPastaQ.jl")

# ╔═╡ e090f643-b236-4ed6-97ca-fa327b2f69d8
N = 2

# ╔═╡ 5c241661-9689-468e-a337-36e27c9bb641
nshots=10

# ╔═╡ 3cb7b2ce-b483-444a-ac84-629a43c21197
yao_circ = chain(N,put(1=>H),cnot(1,2))

# ╔═╡ da9ab1d5-cbbf-40e4-ad88-aa6aea2bd688
YaoPlots.plot(yao_circ)

# ╔═╡ 62f2e913-398d-444f-be37-8a5aeb3f8a3c
samples = Yao.measure(zero_state(N) |> yao_circ; nshots=nshots)
print("yao: ", samples, "\n")

# ╔═╡ 3bafcfd4-4621-4c4d-9adb-5a3012e43dd0
statev = statevec(zero_state(N) |> yao_circ)

# ╔═╡ 7c4cf3ba-9215-4bb8-a3f5-2114baeda1ac
pasta_circ = YaoPastaQ.genlist(yao_circ)

# ╔═╡ 30530e0c-2f85-4b4e-9fa1-6f85747effe6
bases = randombases(N, nshots; local_basis = ["Z"])

# ╔═╡ d28845d1-45f1-43de-8c37-c828fa84e214
ψ = runcircuit(N, pasta_circ, noise = ("amplitude_damping", (γ = 0.01,)))

# ╔═╡ 3da7ea86-857c-48ae-86b2-770f8bbd0485
pasta_data = getsamples(ψ, bases)
print("pasta: ", pasta_data, "\n")

# ╔═╡ 778ac036-cd3d-4395-8322-bc018271f454
processed_pasta = []
for pasta in eachrow(pasta_data)
    # print(pasta[1].second, pasta[2].second, " ")
    push!(processed_pasta, BitArray([pasta[1].second pasta[2].second]))
end
print("\n")

print(processed_pasta)


print("\n")

translated_bitstring = YaoExtensions.BitBasis.BitStr64{2}.(processed_pasta)

# function pasta_samples_to_yao(pasta_samples):

    
# end

# ╔═╡ Cell order:
# ╠═9662f1cc-30bb-11ec-3315-db727be10430
# ╠═75cb8628-92f2-48d0-b344-caf72f355524
# ╠═e090f643-b236-4ed6-97ca-fa327b2f69d8
# ╠═5c241661-9689-468e-a337-36e27c9bb641
# ╠═3cb7b2ce-b483-444a-ac84-629a43c21197
# ╠═da9ab1d5-cbbf-40e4-ad88-aa6aea2bd688
# ╠═62f2e913-398d-444f-be37-8a5aeb3f8a3c
# ╠═3bafcfd4-4621-4c4d-9adb-5a3012e43dd0
# ╠═7c4cf3ba-9215-4bb8-a3f5-2114baeda1ac
# ╠═30530e0c-2f85-4b4e-9fa1-6f85747effe6
# ╠═d28845d1-45f1-43de-8c37-c828fa84e214
# ╠═3da7ea86-857c-48ae-86b2-770f8bbd0485
# ╠═778ac036-cd3d-4395-8322-bc018271f454
