library(hexSticker)

obj <- sticker(
       subplot = "man/figures/capt_1_screen.png",  # image de fond
       s_x = 1, s_y = 1.25, s_width = 0.55,    # position et taille de l'image

       package = "mix",          # nom du package affiché
       p_size = 11, p_y = 1.2, p_x=0.97,          # taille et position du texte
       p_color = "#FFFFFF",#,         # couleur du texte

       url = "WCE",               # Texte du bas
       u_size = 30,
       u_y = 0.47, u_x=0.58,
       u_color = "#222831",
       u_angle = 0,

       h_fill = "#F28C28",          # couleur de fond
       h_color = "#222831",         # couleur du contour

       filename = "man/figures/stickers_mixWCE.png"  # nom du fichier de sortie
       #,spotlight = TRUE
       #,l_x=	1,l_y =1,  l_width=5, l_alpha=0.6	#spotlight parameters
        )
obj

