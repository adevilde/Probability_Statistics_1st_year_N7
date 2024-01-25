
% TP2 de Statistiques : fonctions a completer et rendre sur Moodle
% Nom : Devilder
% Prénom : Alice
% Groupe : 1SN-M

function varargout = fonctions_TP2_stat(varargin)

    [varargout{1},varargout{2}] = feval(varargin{1},varargin{2:end});

end

% Fonction centrage_des_donnees (exercice_1.m) ----------------------------
function [x_G, y_G, x_donnees_bruitees_centrees, y_donnees_bruitees_centrees] = ...
                centrage_des_donnees(x_donnees_bruitees,y_donnees_bruitees)
    
     x_G = mean(x_donnees_bruitees);
     y_G = mean(y_donnees_bruitees);

     x_donnees_bruitees_centrees = x_donnees_bruitees - x_G;
     y_donnees_bruitees_centrees = y_donnees_bruitees - y_G;
     
end

% Fonction estimation_Dyx_MV (exercice_1.m) -------------------------------
function [a_Dyx,b_Dyx] = ...
           estimation_Dyx_MV(x_donnees_bruitees,y_donnees_bruitees,n_tests)

    [x_G, y_G, x_donnees_bruitees_centrees, y_donnees_bruitees_centrees] = ...
        centrage_des_donnees(x_donnees_bruitees,y_donnees_bruitees);

    psi_alea = (rand(n_tests,1)*2 - 1)* pi/2;

    x = repmat(x_donnees_bruitees_centrees,n_tests,1);
    y = repmat(y_donnees_bruitees_centrees,n_tests,1);
    psi = repmat(psi_alea, 1, length(x_donnees_bruitees_centrees));

    somme = sum((y - tan(psi).*x).^2, 2);
    [~,indice] = min(somme);

    psi_estime = psi(indice);
    a_Dyx = tan(psi_estime);
    b_Dyx = y_G - a_Dyx.*x_G;

end

% Fonction estimation_Dyx_MC (exercice_2.m) -------------------------------
function [a_Dyx,b_Dyx] = ...
                   estimation_Dyx_MC(x_donnees_bruitees,y_donnees_bruitees)

    A = [x_donnees_bruitees; ones(1,length(x_donnees_bruitees))]';
    B = y_donnees_bruitees';

    X_estime = (inv(A'*A)*A')*B;
    a_Dyx = X_estime(1);
    b_Dyx = X_estime(2);
    
end

% Fonction estimation_Dorth_MV (exercice_3.m) -----------------------------
function [theta_Dorth,rho_Dorth] = ...
         estimation_Dorth_MV(x_donnees_bruitees,y_donnees_bruitees,n_tests)

    [x_G, y_G, x_donnees_bruitees_centrees, y_donnees_bruitees_centrees] = ...
                centrage_des_donnees(x_donnees_bruitees,y_donnees_bruitees);
    
    theta_alea = (rand(n_tests,1)*2 - 1)* pi/2;
    
    x = repmat(x_donnees_bruitees_centrees,n_tests,1);
    y = repmat(y_donnees_bruitees_centrees,n_tests,1);
    theta = repmat(theta_alea, 1, length(x_donnees_bruitees_centrees));

    somme = sum((x.*cos(theta) + y.*sin(theta)).^2,2);

    [~, indice] = min(somme);

    theta_Dorth = theta(indice);
    rho_Dorth = x_G*cos(theta_Dorth) + y_G*sin(theta_Dorth);

end

% Fonction estimation_Dorth_MC (exercice_4.m) -----------------------------
function [theta_Dorth,rho_Dorth] = ...
                 estimation_Dorth_MC(x_donnees_bruitees,y_donnees_bruitees)

    [x_G, y_G, x_donnees_bruitees_centrees, y_donnees_bruitees_centrees] = ...
                centrage_des_donnees(x_donnees_bruitees,y_donnees_bruitees);
    
    C = [x_donnees_bruitees_centrees; y_donnees_bruitees_centrees]';
   
    [Y_estime, ~] = eig(C'*C);

    theta_Dorth = atan(Y_estime(2)/Y_estime(1));
    rho_Dorth = x_G*cos(theta_Dorth) + y_G*sin(theta_Dorth);
    
end
