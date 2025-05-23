include("files_path.jl")
using Pkg
Pkg.activate(".")
Pkg.instantiate()
Pkg.add("Plots")
Pkg.add("XLSX")
Pkg.add("DataFrames")
Pkg.add("Dates")
Pkg.add("Distributions")
Pkg.status()
using XLSX
using DataFrames
using Dates
using LinearAlgebra
using Distributions
using Plots
using Statistics

xf = XLSX.readxlsx("$input_raw_data/BQ1989_Data.xlsx")
XLSX.sheetnames(xf)
sh = xf["Sheet1"] # Get a reference to worksheet 1, as in the replication package
X = Matrix(sh[3:end, 2:end]) # We create a matrix of data
dates = sh[3:end, 1] # Then a matrix of names for variables
vnames_long = sh[1, 2:end]
vnames = sh[2, 2:end]
sh[1:5, 1:end]  # Debug: we check if the imported dataframe corresponds to the excel file
nvar = length(vnames) # We have exactly two variables
# Converti in matrice numerica (Float64) e sostituisci i missing
data = Array{Float64}(undef, size(X))
# Vediamo se ci sono missing values o stringhe
for i in eachindex(X)
    val = X[i]
    if val isa Float64
        data[i] = val
    elseif val isa Missing
        data[i] = NaN
    elseif val isa AbstractString
        parsed = tryparse(Float64, val)
        data[i] = isnothing(parsed) ? NaN : parsed
    else
        data[i] = NaN
    end
end
data[1:5, 1:end]  # Debug: we check if the imported dataframe corresponds to the excel file
first_date = string(dates[1])      # "1980:1"
year = parse(Int, first_date[1:4]) # "1980" -> 1980
quarter = parse(Int, first_date[end]) # "1" -> 1
nobs = size(data, 1)  # Numero di righe
det = 1
(nobs, nvar) = size(data)
nlags = 8


# Funzione che restituisce un oggetto VAROptions con i valori di default

# Costruisci DataFrame con nomi di colonna
vnames_1d = vec(vnames)  # converte da matrice 1x2 a vettore [ "y", "u" ]
df = DataFrame(data, Symbol.(vnames_1d))

# Numero di osservazioni e variabili
(nobs, nvar) = size(data)

# Definisci la struct che rappresenta le opzioni del VAR
mutable struct VAROptions
    vnames::Union{Vector{String}, Nothing}
    vnames_ex::Union{Vector{String}, Nothing}
    snames::Union{Vector{String}, Nothing}
    nsteps::Int
    impact::Int
    shut::Int
    ident::String
    recurs::String
    ndraws::Int
    mult::Int
    pctg::Int
    method::String
    sr_hor::Int
    sr_rot::Int
    sr_draw::Int
    sr_mod::Int
    pick::Int
    quality::Int
    suptitle::Int
    datesnum::Union{Vector{Float64}, Nothing}
    datestxt::Union{Vector{String}, Nothing}
    datestype::Int
    firstdate::Union{Float64, Nothing}
    frequency::String
    figname::Union{String, Nothing}
    FigSize::Tuple{Int, Int}
end

function VARoption()
    return VAROptions(
        nothing,          # vnames
        nothing,          # vnames_ex
        nothing,          # snames
        40,               # nsteps
        0,                # impact
        0,                # shut
        "short",          # ident
        "wold",           # recurs
        1000,             # ndraws
        10,               # mult
        95,               # pctg
        "bs",             # method
        1,                # sr_hor
        500,              # sr_rot
        100000,           # sr_draw
        1,                # sr_mod
        0,                # pick
        2,                # quality
        0,                # suptitle
        nothing,          # datesnum
        nothing,          # datestxt
        1,                # datestype
        nothing,          # firstdate
        "q",              # frequency
        nothing,          # figname
        (26, 24)          # FigSize
    )
end

struct OLSResult
    meth::String
    y::Vector{Float64}
    x::Matrix{Float64}
    nobs::Int
    nvar::Int
    beta::Vector{Float64}
    yhat::Vector{Float64}
    resid::Vector{Float64}
    sige::Float64
    sigbeta::Matrix{Float64}
    bstd::Vector{Float64}
    bint::Matrix{Float64}
    tstat::Vector{Float64}
    tprob::Vector{Float64}
    rsqr::Float64
    rbar::Float64
    dw::Float64
    add_const::Int
    F::Union{Float64, Nothing}
end

function VARmodel(ENDO::Matrix{Float64}, nlag::Int; add_const::Int=1, EXOG=nothing, nlag_ex::Int=0)
    nobs, nvar = size(ENDO)
    VARopt = VARoption()

    VAR = Dict{Symbol, Any}()
    VAR[:ENDO] = ENDO
    VAR[:nlag] = nlag
    VAR[:add_const] = add_const

    if EXOG !== nothing
        nobs2, nvar_ex = size(EXOG)
        if nobs2 != nobs
            error("var: nobs in EXOG not same as ENDO")
        end
        VAR[:EXOG] = EXOG
    else
        nvar_ex = 0
        nlag_ex = 0
        VAR[:EXOG] = nothing
    end

    nobse = nobs - max(nlag, nlag_ex)
    VAR[:nobs] = nobse
    VAR[:nvar] = nvar
    VAR[:nvar_ex] = nvar_ex
    VAR[:nlag] = nlag
    VAR[:nlag_ex] = nlag_ex

    ncoeff = nvar * nlag
    VAR[:ncoeff] = ncoeff
    ncoeff_ex = nvar_ex * (nlag_ex + 1)
    ntotcoeff = ncoeff + ncoeff_ex + add_const
    VAR[:ntotcoeff] = ntotcoeff

    # Supponiamo che VARmakexy sia definita altrove e funzioni correttamente
    Y, X = VARmakexy(ENDO, nlag, add_const)

    # Gestione EXOG
    if nvar_ex > 0
        X_EX = VARmakelags(EXOG, nlag_ex)
        if nlag == nlag_ex
            X = hcat(X, X_EX)
        elseif nlag > nlag_ex
            diff = nlag - nlag_ex
            X_EX = X_EX[(diff+1):end, :]
            X = hcat(X, X_EX)
        elseif nlag < nlag_ex
            diff = nlag_ex - nlag
            Y = Y[(diff+1):end, :]
            X = hcat(X[(diff+1):end, :], X_EX)
        end
    end

    # Qui servono OLSmodel e tdis_prb definiti altrove
    for j in 1:nvar
        Yvec = Y[:, j]
        OLSout = OLSmodel(Yvec, X, 0)
        tprob = tdis_prb(OLSout.tstat, nobse - ncoeff)
        VAR[Symbol("eq$j")] = Dict(
            :beta => OLSout.beta,
            :tstat => OLSout.tstat,
            :bstd => OLSout.bstd,
            :tprob => tprob,
            :resid => OLSout.resid,
            :yhat => OLSout.yhat,
            :y => Yvec,
            :rsqr => OLSout.rsqr,
            :rbar => OLSout.rbar,
            :sige => OLSout.sige,
            :dw => OLSout.dw
        )
    end

    Ft = (X' * X) \ (X' * Y)
    VAR[:Ft] = Ft
    VAR[:F] = Ft'

    SIGMA = (1/(nobse - ntotcoeff)) * (Y - X * Ft)' * (Y - X * Ft)
    VAR[:sigma] = SIGMA
    VAR[:resid] = Y - X * Ft
    VAR[:X] = X
    VAR[:Y] = Y
    if nvar_ex > 0
        VAR[:X_EX] = X_EX
    end

    # FIX: definisci F prima dell’uso
    F = VAR[:F]
    ncoef_dyn = nvar * nlag
    top = F[:, end - ncoef_dyn + 1:end]

    bottom = hcat(
    Matrix{Float64}(I, nvar*(nlag-1), nvar*(nlag-1)),
    zeros(nvar*(nlag-1), nvar)
    )

    Fcomp = vcat(top, bottom)
    VAR[:Fcomp] = Fcomp
VAR[:maxEig] = maximum(abs.(eigvals(Fcomp)))
    # Inizializzazione placeholder
    VAR[:B] = nothing
    VAR[:Biv] = nothing
    VAR[:PSI] = nothing
    VAR[:Fp] = nothing
    VAR[:IV] = nothing

    return VAR, VARopt
end

function VARmakexy(DATA::Matrix{Float64}, lags::Int, add_const::Int)
    nobs, _ = size(DATA)
    
    # Y: dati da riga (lags+1) fino a fine
    Y = DATA[(lags+1):end, :]
    
    # Costruzione matrice X
    X = Array{Float64,2}(undef, nobs - lags, 0)  # matrice vuota con righe corrette
    
    # Concatenazione delle lagged variables
    for jj in 0:(lags-1)
        # seleziona dati ritardati
        lagged = DATA[(jj+1):(nobs - lags + jj), :]
        # concatena a sinistra (come in Matlab con [DATA(...), X])
        X = hcat(lagged, X)
    end
    
    # Aggiungi costante e trend se richiesto
    if add_const == 1
        X = hcat(ones(nobs - lags), X)  # costante
    elseif add_const == 2
        trend = collect(1:(nobs - lags))
        X = hcat(ones(nobs - lags), trend, X)
    elseif add_const == 3
        trend = collect(1:(nobs - lags))
        X = hcat(ones(nobs - lags), trend, trend.^2, X)
    end
    
    return Y, X
end

function VARmakelags(DATA::Matrix{Float64}, lag::Int)
    nobs, nvar = size(DATA)
    OUT = Array{Float64, 2}(undef, 0, 0)  # matrice vuota

    # Costruisci la matrice laggata (attenzione all'ordine)
    for jj in 0:(lag - 1)
        part = DATA[(jj+1):(nobs - lag + jj), :]
        if size(OUT, 1) == 0
            OUT = part
        else
            OUT = hcat(part, OUT)
        end
    end

    # Prendi la parte non laggata (dati dalla riga lag+1 in poi)
    aux = DATA[(lag+1):end, :]

    # Concatenale orizzontalmente con la matrice laggata
    OUT = hcat(aux, OUT)

    return OUT
end

function OLSmodel(y::Vector{Float64}, x::Matrix{Float64}, add_const::Int=0)
    if isempty(x)
        nobs = length(y)
        nvar = 0
    else
        nobs, nvar = size(x)
        if nobs != length(y)
            error("x and y must have same number of observations")
        end
    end

    # Add constant, trend, or trend^2
    if add_const == 1
        x = hcat(ones(nobs), x)
        nvar += 1
    elseif add_const == 2
        trend = collect(1:nobs)
        x = hcat(ones(nobs), trend, x)
        nvar += 2
    elseif add_const == 3
        trend = collect(1:nobs)
        x = hcat(ones(nobs), trend.^2, x)
        nvar += 3
    end

    # OLS estimation
    xpxi = nobs < 10000 ? inv(qr(x).R' * qr(x).R) : inv(x' * x)
    beta = xpxi * (x' * y)
    yhat = x * beta
    resid = y - yhat

    # Variance and standard errors
    sigu = resid' * resid
    sige = sigu / (nobs - nvar)
    sigbeta = sige * xpxi
    bstd = sqrt.(diag(sigbeta))

    # t-stats and confidence intervals
    tcrit = quantile(TDist(nobs - nvar), 0.975)
    bint = hcat(beta .- tcrit .* bstd, beta .+ tcrit .* bstd)
    tstat = beta ./ bstd
    tprob = 2.0 .* (1.0 .- cdf(TDist(nobs - nvar), abs.(tstat)))

    # R-squared and adjusted R²
    rsqr1 = sigu
    rsqr2 = sum((y .- mean(y)).^2)
    rsqr = 1.0 - rsqr1 / rsqr2
    rsqr1 /= (nobs - nvar)
    rsqr2 /= (nobs - 1)
    rbar = rsqr2 != 0 ? 1.0 - rsqr1 / rsqr2 : rsqr

    # Durbin-Watson statistic
    ediff = resid[2:end] .- resid[1:end-1]
    dw = sum(ediff.^2) / sigu

    # F-statistic (if constant present)
    F = nothing
    if add_const > 0
        fx = x[:, 1]
        fxpxi = inv(fx' * fx)
        fbeta = fxpxi * (fx' * y)
        fyhat = fx * fbeta
        fresid = y - fyhat
        fsigu = fresid' * fresid
        frsqr = 1.0 - fsigu / sum((y .- mean(y)).^2)
        F = ((frsqr - rsqr) / (1 - nvar)) / ((1 - rsqr) / (nobs - nvar))
    end

    return OLSResult(
        "ols", y, x, nobs, nvar,
        beta, yhat, resid, sige, sigbeta,
        bstd, bint, tstat, tprob,
        rsqr, rbar, dw, add_const, F
    )
end

# Metodo per singolo valore t-stat
function tdis_prb(tstat::Float64, dof::Int)
    dist = TDist(dof)
    return 2 * (1 - cdf(dist, abs(tstat)))
end

# Metodo per array di t-stat
function tdis_prb(tstat::AbstractVector{<:Real}, dof::Int)
    dist = TDist(dof)
    return 2 .* (1 .- cdf.(dist, abs.(tstat)))
end


# Stima del VAR
VAR, VARopt = VARmodel(data, nlags; add_const=det)

# Carica le opzioni
VARopt = VARoption()

VARopt.nsteps = 40
VARopt.ident = "long"
vnames_vec = vec(vnames_long)
VARopt.vnames = vnames_vec
VARopt.FigSize = (26, 12)


function nearest_posdef(A)
    B = (A + A') / 2
    _, S, V = svd(B)
    H = V * Diagonal(S) * V'
    A2 = (B + H) / 2
    A3 = (A2 + A2') / 2

    if isposdef(A3)
        return A3
    end

    spacing = eps(Float64)
    I = Matrix{Float64}(I, size(A,1), size(A,2))
    k = 1
    while !isposdef(A3)
        mineig = minimum(eigvals(A3))
        A3 += I * (-mineig * k^2 + spacing)
        k += 1
    end

    return A3
end


function VARir(VAR::Dict{Symbol,Any}, VARopt)
    # Controllo input
    if !haskey(VAR, :IV) && VARopt.ident == "iv"
        error("You need to provide the data for the instrument in VAR (IV)")
    end

    # Parametri
    nsteps = VARopt.nsteps
    impact = VARopt.impact
    shut = VARopt.shut
    recurs = VARopt.recurs
    Fcomp = copy(VAR[:Fcomp])  # meglio copiare per non modificare in-place
    nvar = VAR[:nvar]
    nlag = VAR[:nlag]
    sigma = VAR[:sigma]

    println("DEBUG: nsteps=$nsteps, nvar=$nvar, nlag=$nlag")
    println("DEBUG: size(VAR[:F]) = ", size(VAR[:F]))
    println("DEBUG: size(Fcomp) = ", size(Fcomp))

    IR = Array{Float64}(undef, nsteps, nvar, nvar)  # (horizons, variables, shocks)

    # Calcolo rappresentazione di Wold (PSI multipliers)
    PSI = zeros(nvar, nvar, nsteps)

    # Inizializza Fp (nvar x nvar x nsteps)
    Fp = zeros(nvar, nvar, nsteps)
    I_idx = VAR[:add_const] + 1

    println("DEBUG: Inizio caricamento Fp")
    for ii in 1:nlag
        start_col = I_idx
        end_col = I_idx + nvar - 1
        println("  ii=$ii, I_idx range: $start_col:$end_col")
        if end_col > size(VAR[:F], 2)
            error("Index out of bounds when accessing VAR[:F] at columns $start_col:$end_col")
        end
        Fp[:, :, ii] = VAR[:F][:, start_col:end_col]
        I_idx += nvar
    end

    for ii in nlag+1:nsteps
        Fp[:, :, ii] .= 0
    end
    println("DEBUG: Fp caricato correttamente")

    # Riempi PSI con la matrice identità per il primo step
    PSI[:, :, 1] = Matrix(I, nvar, nvar)

    for ii in 2:nsteps
        aux = zeros(nvar, nvar)
        for jj in 1:(ii-1)
            aux += PSI[:, :, ii-jj] * Fp[:, :, jj]
        end
        PSI[:, :, ii] = aux
    end
    println("DEBUG: Calcolo PSI completato")

    VAR[:PSI] = PSI

    # Identificazione: matrice B
    B = zeros(nvar, nvar)
    if VARopt.ident == "short"
        chol_result = cholesky(sigma)
        B = chol_result.L
    elseif VARopt.ident == "long"
     Finf_big = inv(Matrix(I, nvar * nlag, nvar * nlag) - Fcomp)
     Finf = Finf_big[1:nvar, 1:nvar]

     M = Finf * sigma * Finf'
     M = (M + M') / 2
     M_posdef = nearest_posdef(M)
     chol_factor = cholesky(M_posdef)
     D = chol_factor.L'   # trasposta di L
     B = Finf \ D
    elseif VARopt.ident == "sign"
        if !haskey(VAR, :B)
            error("You need to provide the B matrix with SR or SignRestrictions")
        else
            B = VAR[:B]
        end
    elseif VARopt.ident == "iv"
        error("IV identification not implemented yet in this Julia version")
    elseif VARopt.ident == "exog"
        B[:, 1] = VAR[:F][:, VAR[:const] + VAR[:ncoeff] + 1]
    else
        error("Identification incorrectly specified")
    end
    println("DEBUG: Matrice B calcolata")

    # Calcolo IRFs
    for mm in 1:nvar
        if shut != 0
            Fcomp[shut, :] .= 0
        end

        response = zeros(nvar, nsteps)
        impulse = zeros(nvar)

        if impact == 0
            impulse[mm] = 1.0
        elseif impact == 1
            impulse[mm] = 1 / B[mm, mm]
        else
            error("Impact must be either 0 or 1")
        end

        response[:, 1] = B * impulse

        if shut != 0
            response[shut, 1] = 0
        end

        if recurs == "wold"
            for kk in 2:nsteps
                response[:, kk] = PSI[:, :, kk] * B * impulse
            end
        elseif recurs == "comp"
            for kk in 2:nsteps
                FcompN = Fcomp^(kk-1)
                response[:, kk] = FcompN[1:nvar, 1:nvar] * B * impulse
            end
        else
            error("recurs option must be 'wold' or 'comp'")
        end

        # Salva in IR: dimensioni (nsteps x nvar x nvar)
        IR[:, :, mm] = transpose(response)
        println("DEBUG: IRF calcolata per shock $mm")
    end

    VAR[:B] = B

    println("DEBUG: Funzione VARir completata correttamente")
    return IR, VAR
end

# Calcolo delle funzioni di risposta all'impulso (IRF)
IR, VAR = VARir(VAR, VARopt)

function VARirband(VAR::Dict{Symbol,Any}, VARopt)
    # ===================================================================
    # Calcola bande di confidenza per IRF usando bootstrap e deviazione std
    #
    # INPUT:
    #  - VAR: dizionario risultati VAR
    #  - VARopt: opzioni VAR
    #
    # OUTPUT:
    #  - INF: lower band (nsteps x nvar x nvar)
    #  - SUP: upper band (nsteps x nvar x nvar)
    #  - MED: mean IRF (nsteps x nvar x nvar)
    #  - BAR: mean IRF (idem)
    #
    # ===================================================================

    if VAR === nothing
        error("Devi fornire struttura VAR")
    end
    if VARopt === nothing
        error("Devi fornire VARopt")
    end

    nsteps = VARopt.nsteps
    ndraws = VARopt.ndraws
    pctg = VARopt.pctg
    method = VARopt.method

    Ft = VAR[:Ft]
    nvar = VAR[:nvar]
    nvar_ex = VAR[:nvar_ex]
    nlag = VAR[:nlag]
    nlag_ex = VAR[:nlag_ex]
    add_const = VAR[:add_const]
    nobs = VAR[:nobs]
    resid = VAR[:resid]
    ENDO = VAR[:ENDO]
    EXOG = VAR[:EXOG]
    IV = haskey(VAR, :IV) ? VAR[:IV] : nothing

    INF = zeros(Float64, nsteps, nvar, nvar)
    SUP = zeros(Float64, nsteps, nvar, nvar)
    MED = zeros(Float64, nsteps, nvar, nvar)
    BAR = zeros(Float64, nsteps, nvar, nvar)

    y_artificial = zeros(nobs + nlag, nvar)

    tt = 1
    ww = 1

    IR = zeros(Float64, nsteps, nvar, nvar, ndraws)  # IRF bootstrapate

    while tt <= ndraws
        if tt == VARopt.mult * ww
            println("Loop $tt / $ndraws draws")
            ww += 1
        end

        # Step 1: Genera residui bootstrap
        if method == "bs"
            idxs = rand(1:size(resid,1), nobs)
            u = resid[idxs, :]
        elseif method == "wild"
            if VARopt.ident == "iv"
                rr = 1 .- 2 .* (rand(nobs, size(IV,2)) .> 0.5)
                u = resid .* (rr * ones(size(IV,2), nvar))'
                Z = vcat(IV[1:nlag, :], IV[nlag+1:end, :] .* rr)
            else
                rr = 1 .- 2 .* (rand(nobs) .> 0.5)
                u = resid .* (rr * ones(1, nvar))
            end
        else
            error("Metodo $method non disponibile")
        end

        # Step 2: dati artificiali
        LAG = Float64[]
        for jj in 1:nlag
            y_artificial[jj, :] = ENDO[jj, :]
            LAG = vcat(y_artificial[jj, :], LAG)
        end

        T = collect(1:nobs)

        if add_const == 0
            LAGplus = LAG
        elseif add_const == 1
            LAGplus = vcat(1.0, LAG)
        elseif add_const == 2
            LAGplus = vcat(1.0, T[1], LAG)
        elseif add_const == 3
            LAGplus = vcat(1.0, T[1], T[1]^2, LAG)
        else
            error("Const option must be 0,1,2, or 3")
        end

        if nvar_ex != 0
            LAGplus = vcat(LAGplus, VAR[:X_EX][1, :])
        end

        for jj in (nlag+1):(nobs+nlag)
            for mm in 1:nvar
                y_artificial[jj, mm] = dot(LAGplus, Ft[:, mm]) + u[jj-nlag, mm]
            end

            if jj < nobs + nlag
                LAG = vcat(y_artificial[jj, :], LAG[1:end - nvar])
                if add_const == 0
                    LAGplus = LAG
                elseif add_const == 1
                    LAGplus = vcat(1.0, LAG)
                elseif add_const == 2
                    LAGplus = vcat(1.0, T[jj-nlag], LAG)
                elseif add_const == 3
                    LAGplus = vcat(1.0, T[jj-nlag], T[jj-nlag]^2, LAG)
                end
                if nvar_ex != 0
                    LAGplus = vcat(LAGplus, VAR[:X_EX][jj - nlag, :])
                end
            end
        end

        # Step 3: Stima VAR sui dati artificiali
        VAR_draw = nothing
        if nvar_ex != 0
            VAR_draw, _ = VARmodel(y_artificial, nlag; add_const=add_const, EXOG=EXOG, nlag_ex=nlag_ex)
        else
            VAR_draw, _ = VARmodel(y_artificial, nlag; add_const=add_const)
        end

        if @isdefined(Z)
            VAR_draw[:IV] = Z
        end

        # Step 4: Calcola IRF
        IR_draw, VAR_draw = VARir(VAR_draw, VARopt)

        if VAR_draw[:maxEig] < 0.9999
            IR[:, :, :, tt] = IR_draw
            tt += 1
        end
    end

    println("-- Done!\n")

    # Calcola media e deviazione standard sulle IRF bootstrapate
    B = mean(IR, dims=4)
    BAR = dropdims(B, dims=4)
    tmp = std(IR, dims=4)
    STD = dropdims(tmp, dims=4)

    # Calcola z critico per intervallo di confidenza (es. 95%)
    alpha = (100 - pctg)/100
    z = quantile(Normal(0,1), 1 - alpha/2)  # Requires using Distributions

    # Calcola bande simmetriche intorno alla media
    INF = BAR .- z .* STD
    SUP = BAR .+ z .* STD

    MED = BAR  # mediana qui la sostituiamo con la media

    return INF, SUP, MED, BAR
end


# Calcolo delle bande di errore per l'IRF
IRinf, IRsup, IRmed, IRbar = VARirband(VAR, VARopt)


function VARirplot(IR, VARopt; INF=nothing, SUP=nothing)
    
    if !(hasproperty(VARopt, :vnames) && VARopt.vnames !== nothing && !isempty(VARopt.vnames))
        error("Devi fornire i nomi delle variabili (vnames) in VARopt")
    end
    vnames = VARopt.vnames

    snames = if hasproperty(VARopt, :snames) && VARopt.snames !== nothing && !isempty(VARopt.snames)
        VARopt.snames
    else
        vnames
    end

    nsteps, nvars, nshocks = size(IR)
    pick = hasproperty(VARopt, :pick) ? VARopt.pick : 0
    if !(isa(pick, Integer))
        pick = 0
    end

    if pick < 0 || pick > nvars
        error("Shock selezionato non valido")
    elseif pick == 0
        pick = 1
    else
        nshocks = pick
    end

    for jj in pick:nshocks
        for ii in 1:nvars
            # Qui il tuo codice di plotting
            # esempio: println("Plot var $(vnames[ii]) shock $(snames[jj])")
        end
    end
end

# Chiamata corretta con keyword arguments per INF e SUP
VARirplot(IRbar, VARopt; INF=IRinf, SUP=IRsup)


# Supponiamo che IR sia un array (nsteps, nvar, nshocks)
# e VARopt.nsteps sia definito

cmap = [:blue, :red]  # Colori per le serie

# Imposta dimensione figura
default(size=(1300, 400))  # equivalente a FigSize(26,8) in pollici * 50 (circa)

# Supply shock (subplot 1)
p1 = plot()
plot!(p1, cumsum(IR[:,1,1]), linewidth=2.5, color=cmap[1], label="GDP Level")
plot!(p1, IR[:,2,1], linewidth=2.5, color=cmap[2], label="Unemployment")
plot!(p1, zeros(VARopt.nsteps), linestyle=:dash, color=:black, label="")
title!(p1, "Supply shock")

# Demand shock (subplot 2)
p2 = plot()
plot!(p2, cumsum(-IR[:,1,2]), linewidth=2.5, color=cmap[1], label="GDP Level")
plot!(p2, -IR[:,2,2], linewidth=2.5, color=cmap[2], label="Unemployment")
plot!(p2, zeros(VARopt.nsteps), linestyle=:dash, color=:black, label="")
title!(p2, "Demand shock")

# Combina i due subplot in una riga 1x2
plot(p1, p2, layout = (1, 2), legend = :topright)

# Salvataggio figura (opzionale, sostituisci il path se vuoi)
savefig("$output_figures/Fig_1_2.pdf")

# Pulisce la figura (come clf('reset'))
closeall()

# Colori e stile
c_gdp = :black
line_style = (:solid, :dash, :dash)  # median, inf, sup
lw = 2.5
steps = 1:VARopt.nsteps

# Estrai IRFs
# IRmed, IRinf, IRsup: (nsteps, nvar, nshocks)
# GDP = 1, Unemployment = 2
# Shock 1 = Supply, Shock 2 = Demand

# Output response to Demand (Fig 3)
p1 = plot(steps, IRmed[:, 1, 2], color=c_gdp, lw=lw, linestyle=:solid, label="")
plot!(p1, steps, IRinf[:, 1, 2], color=c_gdp, lw=1.5, linestyle=:dash, label="")
plot!(p1, steps, IRsup[:, 1, 2], color=c_gdp, lw=1.5, linestyle=:dash, label="")
title!(p1, "Figure 3. Output Response to Demand")
xlims!(p1, (0, 40))
ylims!(p1, (-0.2, 1.4))
hline!(p1, [0.0], linestyle=:dash, color=:gray, label="")

# Unemployment response to Demand (Fig 5)
p2 = plot(steps, IRmed[:, 2, 2], color=c_gdp, lw=lw, linestyle=:solid, label="")
plot!(p2, steps, IRinf[:, 2, 2], color=c_gdp, lw=1.5, linestyle=:dash, label="")
plot!(p2, steps, IRsup[:, 2, 2], color=c_gdp, lw=1.5, linestyle=:dash, label="")
title!(p2, "Figure 5. Unemployment Response to Demand")
xlims!(p2, (0, 40))
ylims!(p2, (-0.6, 0.1))
hline!(p2, [0.0], linestyle=:dash, color=:gray, label="")

# Output response to Supply (Fig 4)
p3 = plot(steps, IRmed[:, 1, 1], color=c_gdp, lw=lw, linestyle=:solid, label="")
plot!(p3, steps, IRinf[:, 1, 1], color=c_gdp, lw=1.5, linestyle=:dash, label="")
plot!(p3, steps, IRsup[:, 1, 1], color=c_gdp, lw=1.5, linestyle=:dash, label="")
title!(p3, "Figure 4. Output Response to Supply")
xlims!(p3, (0, 40))
ylims!(p3, (-0.2, 1.4))
hline!(p3, [0.0], linestyle=:dash, color=:gray, label="")

# Unemployment response to Supply (Fig 6)
p4 = plot(steps, IRmed[:, 2, 1], color=c_gdp, lw=lw, linestyle=:solid, label="")
plot!(p4, steps, IRinf[:, 2, 1], color=c_gdp, lw=1.5, linestyle=:dash, label="")
plot!(p4, steps, IRsup[:, 2, 1], color=c_gdp, lw=1.5, linestyle=:dash, label="")
title!(p4, "Figure 6. Unemployment Response to Supply")
xlims!(p4, (0, 40))
ylims!(p4, (-0.2, 0.5))
hline!(p4, [0.0], linestyle=:dash, color=:gray, label="")

# Metti tutto in un layout 2x2
plot(p1, p3, p2, p4, layout=(2,2), size=(1000,700))

# (Opzionale) salva
savefig("$output_figures/Fig_3_4_5_6.pdf")
