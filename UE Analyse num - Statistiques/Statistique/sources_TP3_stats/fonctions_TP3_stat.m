
% TP3 de Statistiques : fonctions a completer et rendre sur Moodle
% Nom : Devilder
% Prenom : Alice
% Groupe : 1SN-M

function varargout = fonctions_TP3_stat(varargin)

    [varargout{1},varargout{2}] = feval(varargin{1},varargin{2:end});

end

% Fonction estimation_F (exercice_1.m) ------------------------------------
function [rho_F,theta_F,ecart_moyen] = estimation_F(rho,theta)

    A = [cos(theta) sin(theta)];
    B = rho;

    X = A\B;

    x_F = X(1);
    y_F = X(2);

    rho_F = sqrt(x_F.^2 + y_F.^2);
    theta_F = atan2(y_F, x_F);

    % A modifier lors de l'utilisation de l'algorithme RANSAC (exercice 2)
    m = length(rho);
    donnees_conformes = abs(rho - rho_F*cos(theta -theta_F));
    ecart_moyen = 1/m*sum(donnees_conformes);

end

% Fonction RANSAC_2 (exercice_2.m) ----------------------------------------
function [rho_F_estime,theta_F_estime] = RANSAC_2(rho,theta,parametres)

    S1 = parametres(1);
    S2 = parametres(2);
    k_max = parametres(3);
    n = length(rho);
    critere_min = Inf;
   
    for i = 1:k_max
        D = randperm(n,2);
        [rho_I,theta_I,~] = estimation_F(rho(D),theta(D));
        donnees_conformes = abs(rho - rho_I*cos(theta -theta_I)) < S1;
        m = sum(donnees_conformes);
        if m/n > S2
            [rho_F,theta_F,ecart_moyen] = estimation_F(rho(donnees_conformes),theta(donnees_conformes));
            if ecart_moyen < critere_min
                rho_F_estime = rho_F;
                theta_F_estime = theta_F;
                critere_min = ecart_moyen;
            end
        end 
    end 

end

% Fonction G_et_R_moyen (exercice_3.m, bonus, fonction du TP1) ------------
function [G, R_moyen, distances] = ...
         G_et_R_moyen(x_donnees_bruitees,y_donnees_bruitees)

    G = [mean(x_donnees_bruitees) mean(y_donnees_bruitees)];
    distances = sqrt((x_donnees_bruitees - G(:,1)).^2 + (y_donnees_bruitees - G(:,2)).^2);
    R_moyen = mean(distances);

end

% Fonction estimation_C_et_R (exercice_3.m, bonus, fonction du TP1) -------
function [C_estime,R_estime,critere] = ...
         estimation_C_et_R(x_donnees_bruitees,y_donnees_bruitees,n_tests,C_tests,R_tests)
     
    % Attention : par rapport au TP1, le tirage des C_tests et R_tests est 
    %             considere comme etant deje effectue 
    %             (il doit etre fait au debut de la fonction RANSAC_3)

    Cx = repmat(C_tests(:,1),1,length(x_donnees_bruitees));
    Cy = repmat(C_tests(:,2),1,length(y_donnees_bruitees));

    x = repmat(x_donnees_bruitees,n_tests,1);
    y = repmat(y_donnees_bruitees,n_tests,1);

    distance_Pi_C = sqrt((x - Cx).^2 + (y - Cy).^2);

    R_dist = repmat(R_tests,1,length(x_donnees_bruitees));

    somme = sum((distance_Pi_C - R_dist).^2 , 2);

    [~,indice] = min(somme);

    C_estime = C_tests(indice,:);
    R_estime = R_tests(indice);

    % critere =

end

% Fonction RANSAC_3 (exercice_3, bonus) -----------------------------------
function [C_estime,R_estime] = ...
         RANSAC_3(x_donnees_bruitees,y_donnees_bruitees,parametres)
     
    % Attention : il faut faire les tirages de C_tests et R_tests ici
    
    [G, R_moyen, ~] = G_et_R_moyen(x_donnees_bruitees,y_donnees_bruitees);

    n_tests = length(x_donnees_bruitees);

    C_tests = (rand(n_tests,2)*2 - 1)*R_moyen + G;
    R_tests = (rand(n_tests,1) + 0.5)*R_moyen;

    critere_min = Inf;
    S1 = parametres(1);
    S2 = parametres(2);
    k_max = parametres(3);
    
    for i = 1:k_max
        cercle = randperm(n_tests,3);
        [C,R] = cercle_3_points(x_donnees_bruitees(cercle),y_donnees_bruitees(cercle));

        donnees_conformes = ((C(1) - x_donnees_bruitees)^2 + (C(2) - y_donnees_bruitees)^2)/R^2 < S1;
        x_C = x_donnees_bruitees(C);
        y_C = y_donnees_bruitees(C);
        m = sum(donnees_conformes);

        if m/n_tests > S2
            [C_estime_test,R_estime_test,critere] = ...
         estimation_C_et_R(x_C,y_C,n_tests,C_tests,R_tests);
            if critere < critere_min
                C_estime = C_estime_test;
                R_estime = R_estime_test;
                critere_min = critere;
            end
        end 
    end 

end
