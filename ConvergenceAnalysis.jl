module ConvergenceAnalysis
using LinearAlgebra, Statistics

export ConvergenceAnalysis

# Erro Médio Quadrático (MSE)
# O MSE é calculado entre múltiplos sinais, considerando o erro quadrático médio em cada índice.
function mse(signals::Vector{Vector{Float64}})
    min_length = minimum(map(length, signals))

    # Calcular o erro médio quadrático
    error = 0.0
    for i in 1:min_length
        values = [signal[i] for signal in signals]
        mean_value = mean(values)
        error += sum((values .- mean_value) .^ 2)
    end
    return error / min_length
end

# Correlação Cruzada
# A correlação cruzada é calculada como a média da correlação entre todos os pares de sinais.
function correlation(signals::Vector{Vector{Float64}})
    min_length = minimum(map(length, signals))

    # Normalizar os sinais e truncar ao menor comprimento
    truncated_signals = [signal[1:min_length] for signal in signals]

    # Calcular correlação cruzada média
    num_signals = length(signals)
    total_correlation = 0.0
    count = 0
    for i in 1:num_signals-1
        for j in i+1:num_signals
            total_correlation += cor(truncated_signals[i], truncated_signals[j])
            count += 1
        end
    end
    return total_correlation / count
end

# Distância Euclidiana
# A distância Euclidiana é calculada como a norma da diferença ponto a ponto entre os sinais.
function euclidean_distance(signals::Vector{Vector{Float64}})
    min_length = minimum(map(length, signals))
    #lengths = map(length, signals)
    # if length(unique(lengths)) != 1
    #     throw(ArgumentError("Todos os sinais devem ter o mesmo comprimento"))
    # end
    distances = zeros(min_length)

    # Calcular a distância Euclidiana para cada índice
    for i in 1:min_length
        values = [signal[i] for signal in signals]
        distances[i] = norm(values .- mean(values))
        # distances[i] = norm(values)  # Distância euclidiana dos valores
    end
    return distances
end

# Convergência ao Longo do Tempo
# Calcular a convergência ao longo do tempo usando uma janela deslizante.
function sliding_window_convergence(signals::Vector{Vector{Float64}}, window_size::Int)
    # Determinar o menor comprimento entre os sinais
    min_length = minimum(map(length, signals))

    # Ajustar o comprimento para a análise com janela
    effective_length = min_length - window_size + 1
    convergence = zeros(effective_length)

    # Calcular a convergência em cada janela
    for start in 1:effective_length
        window_values = [signal[start:start+window_size-1] for signal in signals]
        # Distância Euclidiana na janela
        mean_window = mean(hcat(window_values...), dims=2)
        convergence[start] = sum(norm(signal .- mean_window) for signal in window_values)
    end

    return convergence
end

function convergence_rate(signals::Vector{Vector{Float64}})
    min_length = minimum(map(length, signals))

    # Calcular as distâncias Euclidianas em cada índice
    distances = [norm([signal[i] for signal in signals]) for i in 1:min_length]

    # Ajustar uma curva exponencial: d(t) ≈ d₀ * exp(-λ * t)
    times = 1:min_length
    log_distances = log.(distances)
    coeffs = polyfit(times, log_distances, 1)  # Ajusta uma reta: log(d) ≈ -λt + const
    lambda = -coeffs[1]  # Coeficiente negativo indica a taxa de decaimento

    return lambda
end

function correlation_convergence(signals::Vector{Vector{Float64}})
    min_length = minimum(map(length, signals))

    # Calcular as correlações médias em cada instante
    correlations = zeros(min_length)
    for t in 1:min_length
        values = [signal[t] for signal in signals]
        correlations[t] = cor(values, ones(length(values)))  # Correlação com vetor constante
    end

    # Ajustar uma curva exponencial: corr(t) ≈ 1 - exp(-λ * t)
    times = 1:min_length
    coeffs = polyfit(times, log.(1 .- correlations), 1)  # Ajuste exponencial
    lambda = -coeffs[1]

    return lambda
end

function variance_convergence(signals::Vector{Vector{Float64}})
    min_length = minimum(map(length, signals))

    # Calcular a variância em cada instante
    variances = [var([signal[t] for signal in signals]) for t in 1:min_length]

    # Ajustar uma curva exponencial: var(t) ≈ var₀ * exp(-λ * t)
    times = 1:min_length
    coeffs = polyfit(times, log.(variances), 1)  # Ajusta uma reta: log(var) ≈ -λt + const
    lambda = -coeffs[1]

    return lambda
end




end