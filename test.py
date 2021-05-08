def leer_grafo_archivo(nombreArchivo):
    archivo = open(nombreArchivo)
    grafo = ([],[])
    cantVert = int(archivo.readline().rstrip("\n"))
    lineasRestantes = archivo.readlines()
    # agrego los vertices a la tupla
    for i in range(cantVert):
        grafo[0].append(lineasRestantes[i].rstrip("\n"))
    # cantLineas es cantidad de vertices + cantidad de aristas
    cantLineas = len(lineasRestantes)
    for i in range(cantVert, cantLineas):
        grafo[1].append(tuple(lineasRestantes[i].rstrip("\n").split()))
    return grafo

def main():
    grafo = leer_grafo_archivo("grafotest")
    print(grafo)

if __name__ == '__main__':
    main()
