# HPC-Mandelbrot

## Objectif:
- Optimiser la parallélisation OpenMPI
- Effectuer des essaies avec différentes tailles de cellule
- Comparer le temps d'exécution avec la version naïve

## Modifications apportés:
- Utilisation de la méthode master - worker et les 4 étapes :
- Envoie du travail à chaque tâche
- Réception des résultats et envoie d'un nouveau travail
- Réception des derniers résultats attendu
- Tuer les tâches
- La structure MandelbrotParams utilisé pour la génération de l'image utilise comme delta, la taille d'une cellule
- Envoie de la position de départ de la cellule à chaque tâche
- La fonction compute (ligne 108) utilise maintenant delta_i au lieu de NX pour la création du tableau data
- Possibilité d'envoyé en argument de l'application la taille de l'image et d'une cellule
- Gestions des erreurs suivantes:
- Si la taille d'une cellule est plus grande que la taille de l'image, on ajuste la taille de la cellule à la taille de l'image
- Si la division entière de la taille de l'image par la taille d'une cellule n'est pas 0, on retourne une erreur
- Si le nombre de processeurs alloué à l'application est supérieur au nombre de tâches réellement nécessaire, on utilise seulement la tâche nécessaire en calculant lors de l'initialisation le nombre de processeur réellement nécessaire (évite d'attendre en boucle une réponse d'une tâche ne travaillant pas)

## Conclusion:
- Une taille de cellule trop petite engendre plus d'envoi de message et de message en attente de réception.
  Il y a plus d'attente de la part des workers avant l'envoi d'un nouveau travail. Le master est surchargé de travail.
- Une taille trop grande engendre moins d'envoi de message et moins de message en attente de réception.
  Il y a moins d'attente de la part des workers et moins de travail pour le master mais les workers sont surchargé de travail en particulier. On rejoint le problème de l'algorithme naïf, des tâches beaucoup plus lourdes en calcul que d'autres.
- La taille de cellule idéal devrait permettre d'atteindre une moyenne de message en attente de réception d'au plus 1, afin qu'il y ait toujours une tâche en attente d'un travail pendant que les autres travaillent.
- Dans notre cas, l'équilibre n'est pas parfait mais est bien mieux partagé que dans la version naïve. Pour obtenir une version mieux équilibrée, il faudrait connaitre la charge de calcul probable ou exact de chaque travail donné et de créer une file de travaux à l'avance pour chaque tâche.