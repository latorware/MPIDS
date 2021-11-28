ITERACIONS <- read.delim("clipboard")

ITERACIONS_TRANSFORMADA <- ITERACIONS         #per fer la grafica
ITERACIONS_TRANSFORMADA$tempstransformat <- ITERACIONS_TRANSFORMADA$Temps  * 63



#COMPARACIO RESULTATS i TEMPS SEGONS ITERACIONS DE SIMULATED ANNEALING DE FACEBOOK

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
  )  +
  theme(legend.title=element_blank()) + theme(axis.title.x=element_blank(),
                                              axis.text.x=element_blank(),
                                              axis.ticks.x=element_blank())




#COMPARACIO RESULTATS i TEMPS SEGONS ITERACIONS DE SIMULATED ANNEALING DE CONDMAT

ITERACIONS2 <- read.delim("clipboard")

ITERACIONS2_TRANSFORMADA <- ITERACIONS2         #per fer la grafica
ITERACIONS2_TRANSFORMADA$tempstransformat <- ITERACIONS2_TRANSFORMADA$Temps  * 63


ggplot(data=ITERACIONS2_TRANSFORMADA, aes(x=Iteracions, y=Resultat))+
  geom_line(size=1.2, ) +
  
  theme_ft_rc() + 
  labs( x="Iteracions", y="Resultat (numero vertexs conjunt d.  de influencia positiva)", 
        title="Comparació dels resultats i temps segons el nombre d'iteracions de s annealing",
        subtitle="stiter = 100; k = 1; lambda = (1*10^(-9))",
        caption="Dades: results.csv. graf:graph_CA-CondMat.txt"
  ) +
  scale_color_manual(values = c("#579199"), ) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 25000), sec.axis = sec_axis(~ .*1/63, name = "Temps (en segons)")) + 
  geom_line(aes(y = tempstransformat), color="#ad7d8b",  size=1.2,) +
  theme(legend.title=element_blank(), axis.title.y = element_text(color = "#579199", size = 11), axis.title.y.right = element_text(color = "#ad7d8b",size = 11 ))




#COMPARACIO PAN OFICIAL, GREEDY, ANNEALING TRIAT, ANNEALING MAX


#COMPARACIO TEMPS TOTS

TOTCOMPARACIORESULTATS <- read.delim("clipboard")

ggplot(data=TOTCOMPARACIORESULTATS, aes( x = Graf,y = Temps, fill = Tipus),)+
  geom_bar(position="dodge", size=1, stat='identity') +
  
  theme_ft_rc() + 
  labs( x="Graf", y="Temps (en segons)", 
        title="Comparació del temps d'execució dels nostres alg. respecte el Pan's official",
        subtitle="Paràmetres comuns simulated: stiter = 100, k = 1; lambda = (1*10^(-9))",
        caption="Dades: results.csv."
  ) +
  theme(legend.title=element_blank()) + theme(axis.text.x = element_text(angle = 20, vjust = 1, hjust=1, size = 8))



#COMPARACIO TEMPS OFICIAL I GREEDY

COMPARACIORESULTATSGREDYOFICIAL <- read.delim("clipboard")

ggplot(data=COMPARACIORESULTATSGREDYOFICIAL, aes( x = Graf,y = Temps, fill = Tipus),)+
  geom_bar(position="dodge", size=1, stat='identity') +
  
  theme_ft_rc() + 
  labs( x="Graf", y="Temps (en segons)", 
        title="Comparació del temps d'execució del nostre greedy respecte el Pan's official",
        caption="Dades: results.csv."
  ) +
  theme(legend.title=element_blank()) + theme(axis.text.x = element_text(angle = 20, vjust = 1, hjust=1, size = 8))


#COMPARACIO RESULTATS TOTS

ggplot(data=TOTCOMPARACIORESULTATS, aes( x = Graf,y = Resultat, fill = Tipus),)+
  geom_bar(position="dodge", size=1, stat='identity') +
  
  theme_ft_rc() + 
  labs( x="Graf", y="Resultat (nombre de vertexs del conjunt d de influcencia pos)", 
        title="Comparació de els resultats dels nostres alg. respecte el Pan's official",
        subtitle="Paràmetres comuns simulated: stiter = 100, k = 1; lambda = (1*10^(-9))",
        caption="Dades: results.csv."
  ) +
  theme(legend.title=element_blank()) + theme(axis.text.x = element_text(angle = 20, vjust = 1, hjust=1, size = 8))
