
% TP1 de Statistiques : fonctions a completer et rendre sur Moodle
% Nom : Devilder
% Prénom : Alice
% Groupe : 1SN-M

function varargout = fonctions_TP1_stat(varargin)

    [varargout{1},varargout{2}] = feval(varargin{1},varargin{2:end});

end

% Fonction G_et_R_moyen (exercice_1.m) ----------------------------------
function [G, R_moyen, distances] = ...
         G_et_R_moyen(x_donnees_bruitees,y_donnees_bruitees)

    G = [mean(x_donnees_bruitees) mean(y_donnees_bruitees)];
    distances = sqrt((x_donnees_bruitees - G(:,1)).^2 + (y_donnees_bruitees - G(:,2)).^2);
    R_moyen = mean(distances);


end

% Fonction estimation_C_uniforme (exercice_1.m) ---------------------------
function [C_estime, R_moyen] = ...
         estimation_C_uniforme(x_donnees_bruitees,y_donnees_bruitees,n_tests)
     
    [G, R_moyen, ~] = G_et_R_moyen(x_donnees_bruitees,y_donnees_bruitees);

    C = (rand(n_tests,2)*2 - 1)*R_moyen + G;

    Cx = repmat(C(:,1),1,length(x_donnees_bruitees));
    Cy = repmat(C(:,2),1,length(y_donnees_bruitees));

    x = repmat(x_donnees_bruitees,n_tests,1);
    y = repmat(y_donnees_bruitees,n_tests,1);

    distance_Pi_C = sqrt((x - Cx).^2 + (y - Cy).^2);

    somme = sum((distance_Pi_C - R_moyen).^2 , 2);

    [~,indice] = min(somme);

    C_estime = C(indice,:);

end

% Fonction estimation_C_et_R_uniforme (exercice_2.m) ----------------------
function [C_estime, R_estime] = ...
         estimation_C_et_R_uniforme(x_donnees_bruitees,y_donnees_bruitees,n_tests)

    [G, R_moyen, ~] = G_et_R_moyen(x_donnees_bruitees,y_donnees_bruitees);

    C = (rand(n_tests,2)*2 - 1)*R_moyen + G;
    R = (rand(n_tests,1) + 0.5)*R_moyen;

    Cx = repmat(C(:,1),1,length(x_donnees_bruitees));
    Cy = repmat(C(:,2),1,length(y_donnees_bruitees));

    x = repmat(x_donnees_bruitees,n_tests,1);
    y = repmat(y_donnees_bruitees,n_tests,1);

    distance_Pi_C = sqrt((x - Cx).^2 + (y - Cy).^2);

    R_dist = repmat(R,1,length(x_donnees_bruitees));

    somme = sum((distance_Pi_C - R_dist).^2 , 2);

    [~,indice] = min(somme);

    C_estime = C(indice,:);
    R_estime = R(indice);

end

% Fonction occultation_donnees (donnees_occultees.m) ----------------------
function [x_donnees_bruitees, y_donnees_bruitees] = ...
         occultation_donnees(x_donnees_bruitees, y_donnees_bruitees, theta_donnees_bruitees)
    
    theta_1 = rand()*2*pi;
    theta_2 = rand()*2*pi;
    
    if theta_2 >= theta_1 
        indice = find(theta_donnees_bruitees >= theta_1 & theta_donnees_bruitees <= theta_2);
    else
        indice = find(theta_donnees_bruitees >= theta_1 | theta_donnees_bruitees <= theta_2);
    end

    if length(indice) == 0
        indice = 1
    end 
    x_donnees_bruitees = x_donnees_bruitees(indice);
    y_donnees_bruitees = y_donnees_bruitees(indice);    

end

% Fonction estimation_C_et_R_normale (exercice_4.m, bonus) ----------------
function [C_estime, R_estime] = ...
         estimation_C_et_R_normale(x_donnees_bruitees,y_donnees_bruitees,n_tests)



end
