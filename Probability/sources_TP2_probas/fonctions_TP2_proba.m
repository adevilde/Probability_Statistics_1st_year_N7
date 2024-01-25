
% TP2 de Probabilites : fonctions a completer et rendre sur Moodle
% Nom : Devilder
% Prenom : Alice
% Groupe : 1SN-M

function varargout = fonctions_TP2_proba(varargin)

    switch varargin{1}
        
        case {'calcul_frequences_caracteres','determination_langue','coeff_compression','gain_compression','partitionnement_frequences'}

            varargout{1} = feval(varargin{1},varargin{2:end});
            
        case {'recuperation_caracteres_presents','tri_decroissant_frequences','codage_arithmetique'}
            
            [varargout{1},varargout{2}] = feval(varargin{1},varargin{2:end});
            
        case 'calcul_parametres_correlation'
            
            [varargout{1},varargout{2},varargout{3}] = feval(varargin{1},varargin{2:end});
            
    end

end

% Fonction calcul_frequences_caracteres (exercice_0.m) --------------------
function frequences = calcul_frequences_caracteres(texte,alphabet)
    
    freq_abs = zeros(size(alphabet));

    for i = 1:length(alphabet)
        lettre = find(texte == alphabet(i));
        freq_abs(i) = freq_abs(i) + length(lettre);
    end 

    frequences = (1/length(texte))*freq_abs;

end

% Fonction recuperation_caracteres_presents (exercice_0.m) ----------------
function [selection_frequences,selection_alphabet] = recuperation_caracteres_presents(frequences,alphabet)
    
    selection_frequences = frequences(frequences ~= 0);
    selection_alphabet = alphabet(frequences ~= 0);

end

% Fonction tri_decremental_frequences (exercice_0.m) ----------------------
function [frequences_triees,indices_frequences_triees] = ...
                           tri_decroissant_frequences(selection_frequences)
    [frequences_triees,indices_frequences_triees] = sort(selection_frequences,'descend');

end

% Fonction determination_langue (exercice_1.m) ----------------------------
function indice_langue = determination_langue(frequences_texte, frequences_langue, nom_norme)
    % Note : la variable nom_norme peut valoir 'L1', 'L2' ou 'Linf'
    switch nom_norme 
        case 'L2'
            erreur = sum((frequences_texte - frequences_langue).^2,2);
        case 'L1'
            erreur = sum(abs(frequences_texte - frequences_langue),2);
        case 'Linf'
            erreur = max(abs(frequences_texte - frequences_langue), [], 2);
    end
    
    [~,indice_langue] = min(erreur);
        
end

% Fonction coeff_compression (exercice_2.m) -------------------------------
function coeff_comp = coeff_compression(signal_non_encode,signal_encode)

    coeff_comp = length(signal_non_encode)/length(signal_encode);

end

% Fonction coeff_compression (exercice_2bis.m) -------------------------------
function gain_comp = gain_compression(coeff_comp_avant,coeff_comp_apres)

    gain_comp = length(coeff_comp_apres)/length(coeff_comp_avant);

end

% Fonction partitionnement_frequences (exercice_3.m) ----------------------
function bornes = partitionnement_frequences(selection_frequences)

    borne_inf = 0;
    borne_sup = 1;
    bornes = zeros(1.length(selection_frequences);
    for i = 1:length(selection_frequences)
        j = alphabet(i);
        largeur = borne_sup - borne_inf;
        borne_sup = borne_inf + largeur * bornes(2,j);
        borne_inf = borne_inf + largeur * bornes(1,j);
        nb = randi[borne_inf,borne_sup];
        bornes(i) = nb;
end

% Fonction codage_arithmetique (exercice_3.m) -----------------------------
function [borne_inf,borne_sup] = ...
                       codage_arithmetique(texte,selection_alphabet,bornes)


    
end