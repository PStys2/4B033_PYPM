print(x)  # affiche à l'écran une ou plusieurs variable
type(x)  # type de la variable

float()  # force à avoir un chiffre à virgule
int(x)  # force à avoir un chiffre entier
str(x)  # indique qu'il s'agit d'une liste de charactères

dir(" ")  #
help(" ".upper)

if __name__ == '__main__':  # si ce code est exécuté en tant que script principal (appelé directement avec Python et pas importé), alors exécuter cette fonction.
    test_2()

st = dataset.upper()  # passe toutes les données en majuscules
st.count('A')  #

fileObject.read()  # This method returns the bytes read in string
read().split('\n')  # enleve les retours à la ligne

f.readline().rstrip()

current_name = l[1:].strip()

" ".join(["%d" % x for x in counting_nucleotides(s)])

k, m, n = map(float, (k, m,
                      n))  # map() applique la fonction float (ou une autre) pour chaque valeur k, m et n  # c'est l'équivalent d'une boucle for   k, m, n = (float(i) for i in (k, m, n))
map(rna, s)
