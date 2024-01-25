
% TP1 de Probabilites : fonctions a completer et rendre sur Moodle
% Nom : Devilder
% Prénom : Alice
% Groupe : 1SN-M

function varargout = fonctions_TP1_proba(varargin)

    switch varargin{1}     
        case 'ecriture_RVB'
            varargout{1} = feval(varargin{1},varargin{2:end});
        case {'vectorisation_par_colonne','decorrelation_colonnes'}
            [varargout{1},varargout{2}] = feval(varargin{1},varargin{2:end});
        case 'calcul_parametres_correlation'
            [varargout{1},varargout{2},varargout{3}] = feval(varargin{1},varargin{2:end}); 
    end

end

% Fonction ecriture_RVB (exercice_0.m) ------------------------------------
% (Copiez le corps de la fonction ecriture_RVB du fichier du meme nom)
function image_RVB = ecriture_RVB(image_originale)

nb_lignes = size(image_originale,1)/2;
nb_colonnes = size(image_originale,2)/2;

image_RVB = zeros(nb_lignes,nb_colonnes,3);

image_RVB(:,:,1) = image_originale(1:2:end,2:2:end);
image_RVB(:,:,3) = image_originale(2:2:end,1:2:end);

V1 = zeros(nb_lignes,nb_colonnes);
V1 = image_originale(1:2:end,1:2:end);

V2 = zeros(nb_lignes,nb_colonnes);
V2 = image_originale(2:2:end,2:2:end);

image_RVB(:,:,2) = (V1+V2)/2;

end

% Fonction vectorisation_par_colonne (exercice_1.m) -----------------------
function [Vd,Vg] = vectorisation_par_colonne(I)

vd = I(1:end,2:end);
vg = I(1:end,1:end-1);
Vg = vg(:);
Vd = vd(:);

end

% Fonction calcul_parametres_correlation (exercice_1.m) -------------------
function [r,a,b] = calcul_parametres_correlation(Vd,Vg)

moy_vd = sum(Vd)/size(Vd,1);
moy_vg = sum(Vg)/size(Vg,1);

sigma_vd2 = 1/size(Vd,1) * Vd'*Vd - moy_vd*moy_vd;
sigma_vg2 = 1/size(Vg,1) * Vg'*Vg - moy_vg*moy_vg;

sigma_vd = sqrt(sigma_vd2);
sigma_vg = sqrt(sigma_vg2);

cov = 1/size(Vd,1) * Vd'*Vg - moy_vd*moy_vg;

r = cov / (sigma_vd*sigma_vg);

a = cov/(sigma_vd*sigma_vd);
b = moy_vg - a*moy_vd;

end

% Fonction decorrelation_colonnes (exercice_2.m) --------------------------
function [I_decorrelee,I_min] = decorrelation_colonnes(I,I_max)

I_decorrelee = I;
I_decorrelee(:,2:end) = I(:,2:end) - I(:,1:end-1);
I_min = -I_max;

end



