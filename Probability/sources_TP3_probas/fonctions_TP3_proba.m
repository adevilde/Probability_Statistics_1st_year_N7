
% TP3 de Probabilites : fonctions a completer et rendre sur Moodle
% Nom : Devilder
% Prénom : Alice
% Groupe : 1SN-M

function varargout = fonctions_TP3_proba(varargin)

    switch varargin{1}
        
        case 'matrice_inertie'
            
            [varargout{1},varargout{2}] = feval(varargin{1},varargin{2:end});
            
        case {'ensemble_E_recursif','calcul_proba'}
            
            [varargout{1},varargout{2},varargout{3}] = feval(varargin{1},varargin{2:end});
    
    end
end

% Fonction ensemble_E_recursif (exercie_1.m) ------------------------------
function [E,contour,G_somme] = ...
    ensemble_E_recursif(E,contour,G_somme,i,j,voisins,G_x,G_y,card_max,cos_alpha)
    
    contour(i,j)=0;

    for k = 1 : size(voisins,1)

        i_voisin = i + voisins(k,1);
        j_voisin = j + voisins(k,2);

        if contour(i_voisin,j_voisin)==1

            if size(E,1)<=card_max

                G_somme_norme = norm(G_somme);
                G_ij = [G_x(i_voisin, j_voisin) G_y(i_voisin, j_voisin)];
                G_ij_norme = norm(G_ij);

                if (G_ij/G_ij_norme)*(G_somme/G_somme_norme)' >= cos_alpha
                    E = [E ; [i_voisin j_voisin]];
                    G_somme=G_somme+G_ij;
                    [E, contour, G_somme] = ensemble_E_recursif(E,contour,G_somme,i_voisin,j_voisin,voisins,G_x,G_y,card_max,cos_alpha);
                end
            end
        end
    end
end

% Fonction matrice_inertie (exercice_2.m) ---------------------------------
function [M_inertie,C] = matrice_inertie(E,G_norme_E)
    
    x = E(:,2);
    y = E(:,1);
    pi = sum(G_norme_E);

    x_c = (1/pi)*sum(G_norme_E.*x);
    y_c = (1/pi)*sum(G_norme_E.*y);

    M_11 = (1/pi)*sum (G_norme_E.*((x-x_c).^2));
    M_12 = (1/pi)*sum (G_norme_E.*(x-x_c).*(y-y_c));
    M_21 = (1/pi)*sum (G_norme_E.*(y-y_c).*(x-x_c));
    M_22 = (1/pi)*sum (G_norme_E.*((y-y_c).^2));

    M_inertie = [M_11 M_12; M_21 M_22];
    C = [x_c y_c];

end

% Fonction calcul_proba (exercice_2.m) ------------------------------------
function [x_min,x_max,probabilite] = calcul_proba(E_nouveau_repere,p)

    n = size(E_nouveau_repere,1);

    x_min = min(E_nouveau_repere(:,1));
    x_max = max(E_nouveau_repere(:,1));
    y_min = min(E_nouveau_repere(:,2));
    y_max = max(E_nouveau_repere(:,2));

    N = round((x_max - x_min)*(y_max - y_min));

    probabilite = 1 - binocdf(n-1,N,p);

    
end
