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

end