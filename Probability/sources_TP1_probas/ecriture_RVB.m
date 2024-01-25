function image_RVB = ecriture_RVB(image_originale)

nb_lignes = size(image_originale,1)/2
nb_colonnes = size(image_originale,2)/2

image_RVB = zeros(nb_lignes,nb_colonnes,3)

image_RVB(:,:,1) = image_originale(1:2:end,2:2:end);
image_RVB(:,:,3) = image_originale(2:2:end,1:2:end);

V1 = zeros(nb_lignes,nb_colonnes)
V1 = image_originale(1:2:end,1:2:end)

V2 = zeros(nb_lignes,nb_colonnes)
V2 = image_originale(2:2:end,2:2:end)

image_RVB(:,:,2) = (V1+V2)/2



