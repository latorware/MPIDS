ITERACIONS <- read.delim("clipboard")

ITERACIONS_TRANSFORMADA <- ITERACIONS         #per fer la grafica
ITERACIONS_TRANSFORMADA$tempstransformat <- ITERACIONS_TRANSFORMADA$Temps  * 63



#COMPARACIO RESULTATS i TEMPS SEGONS ITERACIONS DE SIMULATED ANNEALING

ggplot(data=ITERACIONS_TRANSFORMADA, aes(x=Iteracions, y=Resultat))+
  geom_line(size=1.2, ) +
  
  theme_ft_rc() + 
  labs( x="Iteracions", y="Resultat (numero vertexs conjunt d.  de influencia positiva)", 
        title="Comparació dels resultats i temps segons el nombre d'iteracions de s annealing",
        subtitle="stiter = 50; k = 1; lambda = (1*10^(-9))",
        caption="Dades: results.csv. graf:ego-facebook; "
  ) +
  scale_color_manual(values = c("#579199"), ) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 4000), sec.axis = sec_axis(~ .*1/63, name = "Temps (en segons)")) + 
  geom_line(aes(y = tempstransformat), color="#ad7d8b",  size=1.2,) +
  theme(legend.title=element_blank(), axis.title.y = element_text(color = "#579199", size = 11), axis.title.y.right = element_text(color = "#ad7d8b",size = 11 ))


#COMPARACIO RESULTATS SEGONS LES ITERACIONS PER CANVI DE TEMPERATURA ANNEALING

STITER <- read.delim("clipboard")

ggplot(data=STITER, aes(y=Resultat, color=factor(Stiter), linetype=factor(Stiter)),)+
  geom_boxplot() +
  theme_ft_rc() + 
  labs(  y="Resultat (numero vertexs conjunt d.  de influencia positiva)", 
         title="Resultats segons nombre de iteracions per canvi de temperatura de s annealing",
         subtitle="iteracions = 5000; k = 1; lambda = (1*10^(-9))",
         caption="Dades: results.csv. graf:ego-facebook;"
  ) +
  scale_color_manual(values = c("#579199", "white", "#d6a69c"), ) +
  theme(legend.title=element_blank()) + theme(axis.title.x=element_blank(),
                                              axis.text.x=element_blank(),
                                              axis.ticks.x=element_blank())




