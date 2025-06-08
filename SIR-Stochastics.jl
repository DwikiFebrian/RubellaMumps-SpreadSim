### A Pluto.jl notebook ###
# v0.20.8

using Markdown
using InteractiveUtils

# ╔═╡ 3ac52b9f-474f-4b06-8cfa-3cee29d56cae
using Pkg

# ╔═╡ 8fde9a00-273c-47e5-92ff-597dc4071eb2
Pkg.add(url="https://github.com/sdwfrost/Gillespie.jl")

# ╔═╡ 770fc43e-6e55-46fa-8f9a-1465e6bc356d
Pkg.add(["Plots", "DataFrames", "Gadfly", "Interact", "Statistics", "Random"])

# ╔═╡ 2f9f7d92-5c99-4475-8cf1-4b9751ebdd47
using Gadfly

# ╔═╡ ae7f5cca-d6e1-4366-87d9-0b1100130e1e
using Interact

# ╔═╡ 6d3df3f7-4dee-4263-b447-a38b12c8082d
using Statistics

# ╔═╡ 22ac5aea-0e87-4f02-811d-2afa1232fe33
using Random

# ╔═╡ 60d62e55-3472-4bac-b9d9-02bc455cbe5a
using Gillespie

# ╔═╡ 9aa222c2-9300-40ec-8f81-b7fb56a34d23
using DataFrames

# ╔═╡ b587f1c0-37e1-11f0-22a9-f31874e613cf
md"# Stochastic SIR Model"

# ╔═╡ 0182a884-b92b-4db5-9499-04d48da45687
md"## Siapin Package"

# ╔═╡ ab784787-2d89-4296-9289-717dff1bd7a1
md"## Model"

# ╔═╡ 0c985392-dbc4-4258-99de-e3e440a22cb8
md""" 
Model yang digunakan adalah SIR, yakni model sederhana yang biasa digunakan untuk menjelaskan persebaran penyakit. \

Dalam diskrit:

$\begin{aligned}
	S_{t+1} &= S_t -\beta S_t I_t \\
	I_{t+1} &= I_t + \beta S_t I_t - \gamma I_t \\
	R_{t+1} &= R_t + \gamma I_t
\end{aligned}$

Dalam kontinu: 

$\begin{aligned}
	\frac{dS(t)}{dt} &= -\beta S(t) I(t) \\
	\frac{dI(t)}{dt} &= \beta S(t) I(t) - \gamma I(t) \\
	\frac{dR(t)}{dt} &= \gamma I(t)
\end{aligned}$

dengan: \
S(t) adalah banyaknya orang yang "mungkin" terinfeksi. \
I(t) adalah banyaknya orang yang terinfeksi. \
R(t) adalah banyaknya orang yang sembuh. \


"""

# ╔═╡ 79f0ef78-20b9-48fa-9964-8056a8f2d32b
md" Model SIR punya ilustrasi sebagai berikut: 

$\begin{aligned}
	S \rightarrow I \rightarrow R
\end{aligned}$
atau:

$\begin{aligned}
	S \rightarrow S - 1; I \rightarrow I + 1 
\end{aligned}$ dengan rate $\beta S(t) I(t)/N(t)$

$\begin{aligned}
	I \rightarrow I - 1; I \rightarrow R + 1 
\end{aligned}$ dengan rate $\gamma I(t)$
"

# ╔═╡ e5e413fc-cca6-4746-949e-6441690121f1
md"## Definisikan Fungsi"

# ╔═╡ 96689eb5-217e-4c5d-accb-3e67e3fd4c1f
function F(x,parms)
  (S,I,R) = x
  (β,γ) = parms
  infection = β*S*I/(S+I+R)
  recovery = γ*I
  [infection,recovery]
    end;

# ╔═╡ db6bbd49-fc28-473b-833d-b683dc3c93d4
md"## Set Nilai Awal"

# ╔═╡ a1eb520b-569e-4fee-9cac-349c22df1c5a
md"x0 adalah kondisi orang di t = 0"

# ╔═╡ 8c197a89-27a8-4e59-89d9-cf51cfa5fc30
x0 = [999, 1, 0] # Ada 1000 orang, 999 sehat, 1 terjangkit.

# ╔═╡ 2e9076e1-a18d-40c8-ab04-fb6437f8d8ed
md" nu adalah kejadian pada model, yakni - 1 orang sehat +1 orang terinfeksi, dan - 1 oeang terinfeksi + 1 orang sembuh."

# ╔═╡ 2f588020-63eb-46bc-b8ac-6e9c275e58db
nu = [[-1 1 0];[0 -1 1]] # Kondisi perubahan populasi

# ╔═╡ c9810c23-7a36-4bf1-aa6b-61c8d39301e6
md" $\beta$ adalah peluang seorang yang sembuh bisa tertular. $\beta$ bisa bisa beda tergantung penyakit."

# ╔═╡ df74db5a-5dc8-4d29-adf0-571b15746bfc
# ╠═╡ disabled = true
#=╠═╡
beta = 0.1/10000
  ╠═╡ =#

# ╔═╡ 9adb34f8-71cb-4bd8-bf08-44da1b52d0da
md" $\gamma$ adalah peluang rata-rata seorang yang terinsfeksi dapat sembuh. $\gamma$ bisa beda tergantung apa penyakitnya."

# ╔═╡ 4cd27a3d-2300-49eb-bc8d-89d339a1c718
# ╠═╡ disabled = true
#=╠═╡
gamma = 0.05
  ╠═╡ =#

# ╔═╡ 8910d1e3-435e-43e6-b65d-8c114ad5bc70
tf = 500.0 #rentang waktu simulasi

# ╔═╡ 894dffa6-0fe3-45b0-81d9-76b65dd67c4b
md"## Basic Reproduction Number/R0"

# ╔═╡ 948683e6-2d83-48eb-a4fd-77ffc9175d1f
md" Jika R₀ > 1, infeksi berpotensi menyebar. \
Jika R₀ < 1, infeksi berpotensi berhenti."

# ╔═╡ a44a4661-5ef7-457e-a47b-de0bc93d14b1
md"## Simulasi Rubella"

# ╔═╡ 300d5e53-cccf-4808-a496-70799dd8d5b5
paramrub = [0.5, 0.077] #[betha, gamma]

# ╔═╡ d53e56a0-f85e-45fa-afb5-24010fbb4561
begin
	Random.seed!(1478)
	resultrub = ssa(x0,F,nu,paramrub,tf)
end

# ╔═╡ 3bd84940-6dbd-4e79-b3e6-4e83616f0d7b
timerub = resultrub.time;

# ╔═╡ 0f6d8a94-cae0-4875-82da-6e2009b92ed7
simrub = resultrub.data;

# ╔═╡ 487d56ab-fcca-4906-870f-329973e9db00
df_rub = DataFrame(time = timerub, S = simrub[:,1], I = simrub[:,2], R = simrub[:,3])

# ╔═╡ d509d538-3d23-4513-9960-38dafada053e
md"## Plot Simulasi Rubella"

# ╔═╡ 90fdc624-324c-4c2d-82ca-71999ca60a59
plot(df_rub,
    layer(x=:time, y=:S, Geom.step, Theme(default_color=colorant"red")),
    layer(x=:time, y=:I, Geom.step, Theme(default_color=colorant"blue")),
    layer(x=:time, y=:R, Geom.step, Theme(default_color=colorant"green")),
    Guide.xlabel("Time"),
    Guide.ylabel("Population"),
    Guide.title("SIR Epidemiological Model for Rubella"),
    Guide.manual_color_key("Compartment", ["S", "I", "R"], ["red", "blue", "green"])
)

# ╔═╡ 52c4328d-eca1-4e6f-9c52-ff4f8535ccad
md"## Analisis Simulasi Rubella"

# ╔═╡ 9c52cf95-e36c-40a5-8b75-f5e9bc230e11
begin
	R0_rub = paramrub[1] / paramrub[2]
	println("Basic reproduction number R₀ untuk Rubella = ", R0_rub)
end

# ╔═╡ cc41c637-d3af-40aa-a2db-dce9898278f6
begin
	peak_I_rub = maximum(simrub[:,2])                   # maksimum I(t)
	idx_peak_rub = findfirst(x -> x == peak_I_rub, simrub[:,2])
	time_peak_rub = timerub[idx_peak_rub]
	
	println("Waktu puncak infeksi: ", time_peak_rub, " satuan waktu")
	println("Jumlah puncak infeksi: ", peak_I_rub)
end

# ╔═╡ eb7b5728-652d-49ce-aef9-f6ef54efcc0a
begin
	# cari indeks pertama setelah peak di mana I = 0
	idx_recovery_rub = findfirst(i -> i >= idx_peak_rub && simrub[i,2] == 0, 1:length(timerub))
	
	if idx_recovery_rub !== nothing
	    time_outbreak_end_rub = timerub[idx_recovery_rub]
	    println("Waktu outbreak selesai (I=0): ", time_outbreak_end_rub, " satuan waktu")
	else
	    println("Outbreak belum selesai pada waktu simulasi terakhir")
	end
end

# ╔═╡ d18cdc71-701a-40ab-bff2-e224e99493ee
md"## Simulasi Gondogan"

# ╔═╡ ff5747d5-4e0e-4734-9194-626f27c68338
paramgon = [0.7, 0.1429] #[betha, gamma]

# ╔═╡ 6369b7aa-fb79-4c3f-87b3-dead8c5cfcb1
begin
	Random.seed!(1478)
	resultgon = ssa(x0,F,nu,paramgon,tf)
end

# ╔═╡ 4d52fe5f-eb5e-4635-8d22-8862453f06ed
timegon = resultgon.time;

# ╔═╡ 2cd2e889-cf65-45bd-bf6c-0989072ae5eb
simgon = resultgon.data;

# ╔═╡ 1759923e-1490-4a2e-845d-576215ba2ddf
df_gon = DataFrame(time = timegon, S = simgon[:,1], I = simgon[:,2], R = simgon[:,3])

# ╔═╡ 20de780a-d28f-4833-baba-cc384619ba4b
md"## Plot Simulasi untuk Gondogan"

# ╔═╡ 5d192d21-1ea7-4793-9287-7ac4a1288bad
plot(df_gon,
    layer(x=:time, y=:S, Geom.step, Theme(default_color=colorant"red")),
    layer(x=:time, y=:I, Geom.step, Theme(default_color=colorant"blue")),
    layer(x=:time, y=:R, Geom.step, Theme(default_color=colorant"green")),
    Guide.xlabel("Time"),
    Guide.ylabel("Population"),
    Guide.title("SIR Epidemiological Model for Mumps"),
    Guide.manual_color_key("Compartment", ["S", "I", "R"], ["red", "blue", "green"])
)

# ╔═╡ 17b9e760-e7b0-49d5-b1d1-2cb0d7ac83e6
md"## Analisis Hasil Simulasi"

# ╔═╡ 9f375394-5799-48fc-901c-117af9acd04c
begin
	R0_gon = paramgon[1] / paramgon[2]
	println("Basic reproduction number R₀ untuk Gondogan = ", R0_gon)
end

# ╔═╡ 73c2ee3c-52b2-48d3-843a-fe426a358d4c
begin
	peak_I_gon = maximum(simgon[:,2])                   # maksimum I(t)
	idx_peak_gon = findfirst(x -> x == peak_I_gon, simgon[:,2])
	time_peak_gon = timerub[idx_peak_gon]
	
	println("Waktu puncak infeksi: ", time_peak_gon, " satuan waktu")
	println("Jumlah puncak infeksi: ", peak_I_gon)
end

# ╔═╡ 39af1d91-2529-43af-bad5-e2694c48d47c
begin
	# cari indeks pertama setelah peak di mana I = 0
	idx_recovery_gon = findfirst(i -> i >= idx_peak_gon && simgon[i,2] == 0, 1:length(timegon))
	
	if idx_recovery_gon !== nothing
	    time_outbreak_end_gon = timerub[idx_recovery_gon]
	    println("Waktu outbreak selesai (I=0): ", time_outbreak_end_gon, " satuan waktu")
	else
	    println("Outbreak belum selesai pada waktu simulasi terakhir")
	end
end

# ╔═╡ Cell order:
# ╠═b587f1c0-37e1-11f0-22a9-f31874e613cf
# ╠═0182a884-b92b-4db5-9499-04d48da45687
# ╠═3ac52b9f-474f-4b06-8cfa-3cee29d56cae
# ╠═8fde9a00-273c-47e5-92ff-597dc4071eb2
# ╠═770fc43e-6e55-46fa-8f9a-1465e6bc356d
# ╠═2f9f7d92-5c99-4475-8cf1-4b9751ebdd47
# ╠═ae7f5cca-d6e1-4366-87d9-0b1100130e1e
# ╠═6d3df3f7-4dee-4263-b447-a38b12c8082d
# ╠═22ac5aea-0e87-4f02-811d-2afa1232fe33
# ╠═60d62e55-3472-4bac-b9d9-02bc455cbe5a
# ╠═9aa222c2-9300-40ec-8f81-b7fb56a34d23
# ╠═ab784787-2d89-4296-9289-717dff1bd7a1
# ╠═0c985392-dbc4-4258-99de-e3e440a22cb8
# ╠═79f0ef78-20b9-48fa-9964-8056a8f2d32b
# ╠═e5e413fc-cca6-4746-949e-6441690121f1
# ╠═96689eb5-217e-4c5d-accb-3e67e3fd4c1f
# ╠═db6bbd49-fc28-473b-833d-b683dc3c93d4
# ╠═a1eb520b-569e-4fee-9cac-349c22df1c5a
# ╠═8c197a89-27a8-4e59-89d9-cf51cfa5fc30
# ╠═2e9076e1-a18d-40c8-ab04-fb6437f8d8ed
# ╠═2f588020-63eb-46bc-b8ac-6e9c275e58db
# ╠═c9810c23-7a36-4bf1-aa6b-61c8d39301e6
# ╠═df74db5a-5dc8-4d29-adf0-571b15746bfc
# ╠═9adb34f8-71cb-4bd8-bf08-44da1b52d0da
# ╠═4cd27a3d-2300-49eb-bc8d-89d339a1c718
# ╠═8910d1e3-435e-43e6-b65d-8c114ad5bc70
# ╠═894dffa6-0fe3-45b0-81d9-76b65dd67c4b
# ╠═948683e6-2d83-48eb-a4fd-77ffc9175d1f
# ╠═a44a4661-5ef7-457e-a47b-de0bc93d14b1
# ╠═300d5e53-cccf-4808-a496-70799dd8d5b5
# ╠═d53e56a0-f85e-45fa-afb5-24010fbb4561
# ╠═3bd84940-6dbd-4e79-b3e6-4e83616f0d7b
# ╠═0f6d8a94-cae0-4875-82da-6e2009b92ed7
# ╠═487d56ab-fcca-4906-870f-329973e9db00
# ╠═d509d538-3d23-4513-9960-38dafada053e
# ╠═90fdc624-324c-4c2d-82ca-71999ca60a59
# ╠═52c4328d-eca1-4e6f-9c52-ff4f8535ccad
# ╠═9c52cf95-e36c-40a5-8b75-f5e9bc230e11
# ╠═cc41c637-d3af-40aa-a2db-dce9898278f6
# ╠═eb7b5728-652d-49ce-aef9-f6ef54efcc0a
# ╠═d18cdc71-701a-40ab-bff2-e224e99493ee
# ╠═ff5747d5-4e0e-4734-9194-626f27c68338
# ╠═6369b7aa-fb79-4c3f-87b3-dead8c5cfcb1
# ╠═4d52fe5f-eb5e-4635-8d22-8862453f06ed
# ╠═2cd2e889-cf65-45bd-bf6c-0989072ae5eb
# ╠═1759923e-1490-4a2e-845d-576215ba2ddf
# ╠═20de780a-d28f-4833-baba-cc384619ba4b
# ╠═5d192d21-1ea7-4793-9287-7ac4a1288bad
# ╠═17b9e760-e7b0-49d5-b1d1-2cb0d7ac83e6
# ╠═9f375394-5799-48fc-901c-117af9acd04c
# ╠═73c2ee3c-52b2-48d3-843a-fe426a358d4c
# ╠═39af1d91-2529-43af-bad5-e2694c48d47c
