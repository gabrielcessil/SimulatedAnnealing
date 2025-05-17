function DELTA_E = calculate_DeltaE(Xinit, Fcusto, Imax)     
    positive_counter = 0; % Contador de transicoes positivas de energia
    DELTAS = zeros(1,Imax); % Transicoes de energia
    N = length(Xinit); % Numero de parametros

    X_new = Xinit; 
    for i = 1:Imax
        fprintf(' Calculo de DELTA_E: %d/%d\n', i, Imax);
        X_old = X_new; 
        X_new = PerturbaSolucao(X_old, N);
        
        custo_old = Fcusto(X_old);
        try
            custo_new = Fcusto(X_new);
        catch
            % Desconsidera a nova solucao em caso de erro
            custo_new = custo_old;
            X_new = X_old; 
            continue
        end
        
        % Se o custo aumentou: anota diferencial e conta transicao
        if(custo_new - custo_old > 0 )
            DELTAS(i) = custo_new - custo_old;
            positive_counter = positive_counter+1;
        end
            
    end
    
    
    DELTA_E = sum(DELTAS)/positive_counter;  % Media aritimetica(para um numero aleatorio de perturbacoes) dos incrementos na funcao objetio
end