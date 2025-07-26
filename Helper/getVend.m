function Vend = getVend(sol)
    % Obtains the applied potential VEND at the final time point of SOL

    Vappt = df_analysis.calcVapp(sol);
    Vend = Vappt(end);

end
