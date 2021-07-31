#! /usr/bin/python

# 6ta Practica Laboratorio 
# Complementos Matematicos I
# Ejemplo parseo argumentos

import argparse
import matplotlib.pyplot as plt
import numpy as np
from numpy.testing._private.utils import tempdir


def leer_grafo_archivo(nombreArchivo):
    archivo = open(nombreArchivo)
    grafo = ([],[])
    cantVert = int(archivo.readline().rstrip("\n"))
    lineasRestantes = archivo.readlines()
    # agrego los vertices a la primer lista de tupla
    for i in range(cantVert):
        grafo[0].append(lineasRestantes[i].rstrip("\n"))
    
    cantLineas = len(lineasRestantes)
    for i in range(cantVert, cantLineas):
        grafo[1].append(tuple(lineasRestantes[i].rstrip("\n").split()))
    return grafo

class LayoutGraph:

    def __init__(self, grafo, temp, iters, refresh, c1, c2, ctemp, ancho, grav, verbose=False):
        """
        Parámetros:
        grafo: grafo en formato lista
        iters: cantidad de iteraciones a realizar
        refresh: cada cuántas iteraciones graficar. Si su valor es cero, entonces debe graficarse solo al final.
        c1: constante de repulsión
        c2: constante de atracción
        verbose: si está encendido, activa los comentarios
        """

        # Guardo el grafo
        self.grafo = grafo

        # Inicializo estado
        self.posicion_X = {}
        self.posicion_Y = {}
        self.accum_X = {}
        self.accum_Y = {}
        self.ancho = ancho
        self.epsilon = 0.5
        self.ctemp = ctemp

        # Guardo opciones
        self.iters = iters
        self.verbose = verbose
        self.refresh = refresh
        self.temp = temp
        self.c1 = c1
        self.c2 = c2
        # las constantes k las guardamos para no calcularlas cada vez que se la necesite
        # consideramos el layout cuadrado (de forma que el area es ancho*ancho)
        self.k_atraccion = c2 * np.sqrt(self.ancho**2 / len(self.grafo[0]))
        self.k_repulsion = c1 * np.sqrt(self.ancho**2 / len(self.grafo[0]))

    #Si verbose está activado, se van imprimiendo los mensajes
    def mostrar_mensaje(self, mensaje):
    	if self.verbose:
    		print(mensaje)

    def randomize_positions(self):
        for vertice in self.grafo[0]:
            self.posicion_X[vertice] = np.random.uniform(0,self.ancho)
            self.posicion_Y[vertice] = np.random.uniform(0,self.ancho)

    def distancia(self, v0, v1):
        return np.sqrt((self.posicion_X[v0] - self.posicion_X[v1])**2 + (self.posicion_Y[v0] - self.posicion_Y[v1])**2)

    def initialize_accumulators(self):
        for vertice in self.grafo[0]:
            self.accum_X[vertice] = 0
            self.accum_Y[vertice] = 0

    def f_attraction(self, d):
        return d**2 / self.k_atraccion

    def f_repulsion(self, d):
        return self.k_atraccion**2 / d

    def division_por_cero(self, distancia, v0, v1):
        while (distancia < self.epsilon):
            fuerzaRand = np.random.uniform(0.0,1.0)
            self.posicion_X[v0] += fuerzaRand
            self.posicion_Y[v0] += fuerzaRand
            self.posicion_X[v1] -= fuerzaRand
            self.posicion_Y[v1] -= fuerzaRand
            distancia = self.distancia(v0, v1)
        return distancia

    def compute_attraction_forces(self):
        for v0,v1 in self.grafo[1]:
            distancia = self.distancia(v0,v1)

            #caso division por cero
            distancia = self.division_por_cero(distancia, v0, v1)
            
            mod_fa = self.f_attraction(distancia)
            fx = (mod_fa * (self.posicion_X[v1] - self.posicion_X[v0])) / distancia
            fy = (mod_fa * (self.posicion_Y[v1] - self.posicion_Y[v0])) / distancia
            self.accum_X[v0] += fx
            self.accum_X[v0] += fy
            self.accum_Y[v1] -= fx
            self.accum_Y[v1] -= fy

    def compute_repulsion_forces(self):
        for v0 in self.grafo[0]:
            for v1 in self.grafo[0]:
                if v0 != v1:
                    distancia = self.distancia(v0,v1)

                    #caso division por cero
                    distancia = self.division_por_cero(distancia, v0, v1)

                    mod_fr = self.f_repulsion(distancia)
                    fx = (mod_fr * (self.posicion_X[v1] - self.posicion_X[v0])) / distancia
                    fy = (mod_fr * (self.posicion_Y[v1] - self.posicion_Y[v0])) / distancia
                    self.accum_X[v0] -= fx
                    self.accum_Y[v0] -= fy
                    self.accum_X[v1] += fx
                    self.accum_Y[v1] += fy

    def update_positions(self):
        self.mostrar_mensaje("Inicio actualizacion posiciones")
        for v in self.grafo[0]:
            moduloFuerza = np.sqrt(self.accum_X[v]**2 + self.accum_Y[v]**2)
            if moduloFuerza > self.temp:
                self.accum_X[v] = (self.accum_X[v] / moduloFuerza) * self.temp
                self.accum_Y[v] = (self.accum_Y[v] / moduloFuerza) * self.temp
            
            self.posicion_X[v] += self.accum_X[v] 
            self.posicion_Y[v] += self.accum_Y[v]

            if self.posicion_X[v] < 0:
                self.mostrar_mensaje("Posición fuera del borde del layout. Cambiando posicion")
                self.posicion_X[v] = 0
            elif self.posicion_X[v] > self.ancho:
                self.mostrar_mensaje("Posición fuera del borde del layout. Cambiando posicion")
                self.posicion_X[v] = self.ancho
            
            if self.posicion_Y[v] < 0:
                self.mostrar_mensaje("Posición fuera del borde del layout. Cambiando posicion")
                self.posicion_Y[v] = 0
            elif self.posicion_Y[v] > self.ancho:
                self.mostrar_mensaje("Posición fuera del borde del layout. Cambiando posicion")
                self.posicion_Y[v] = self.ancho
        
        self.mostrar_mensaje("Fin actualizacion posiciones")

    def update_temperature(self):
        self.mostrar_mensaje("Inicio actualizacion temperatura")
        self.temp *= self.ctemp
        self.mostrar_mensaje("Fin actualizacion temperatura")

    def compute_gravity_forces(self):
        #self.mostrar_mensaje("Inicio calculo gravedad")
        #self.mostrar_mensaje("Fin calculo gravedad")
        pass

    def step(self):
        self.initialize_accumulators()
        self.compute_attraction_forces()
        self.compute_repulsion_forces()
        self.compute_gravity_forces()
        self.update_positions()
        self.update_temperature()

    def dibujar_grafo(self):
        plt.pause(0.1)
        x = [self.posicion_X[v] for v in self.grafo[0]]
        y = [self.posicion_Y[v] for v in self.grafo[0]]
        plt.clf()
        # ejes = plt.gca()
        # ejes.set_xlim([0, self.ancho])
        # ejes.set_ylim([0, self.ancho])
        plt.scatter(x, y)
        for v0, v1 in self.grafo[1]:
            plt.plot((self.posicion_X[v0], self.posicion_X[v1]), 
                     (self.posicion_Y[v0], self.posicion_Y[v1]))

    def layout(self):
        """
        Aplica el algoritmo de Fruchtermann-Reingold para obtener (y mostrar)
        un layout
        """
        self.randomize_positions()
        plt.ion()
        for i in range(1, self.iters):
            self.mostrar_mensaje("Iteracion nro " + str(i))
            self.step()
            self.dibujar_grafo()
        plt.ioff()
        plt.show()
    

def main():
    # Definimos los argumentos de linea de comando que aceptamos
    parser = argparse.ArgumentParser()

    # Verbosidad, opcional, False por defecto
    parser.add_argument(
        '-v', '--verbose',
        action='store_true',
        help='Muestra mas informacion al correr el programa'
    )
    # Cantidad de iteraciones, opcional, 50 por defecto
    parser.add_argument(
        '--iters',
        type=int,
        help='Cantidad de iteraciones a efectuar',
        default=50
    )
    # Temperatura inicial
    parser.add_argument(
        '--temp',
        type=float,
        help='Temperatura inicial',
        default=100.0
    )
    # Cantidad de refreshs
    parser.add_argument(
        '--refresh',
        type=int,
        help='Cada cuantas iteraciones graficar',
        default=0
    )
    # Archivo del cual leer el grafo
    parser.add_argument(
        'file_name',
        help='Archivo del cual leer el grafo a dibujar'
    )

    parser.add_argument(
        '--c1',
        type=float,
        help='Constante de repulsión',
        default=0.1
    )

    parser.add_argument(
        '--c2',
        type=float,
        help='Constante de atracción',
        default=5.0
    )

    # Constante para disminuir temp
    parser.add_argument(
        '--ctemp',
        type = float,
        help = 'Constante para multiplicar por la temperatura',
        default = 0.95
    )

    #Ancho layout
    parser.add_argument(
        '--ancho',
        type=int,
        help='Ancho del layout',
        default=1000
    )

    parser.add_argument(
        '--grav',
        type=float,
        help='Constante de gravedad',
        default=0.05
    )

    args = parser.parse_args()

    # Descomentar abajo para ver funcionamiento de argparse
    # print(args.verbose)
    # print(args.iters)    
    # print(args.file_name)
    # print(args.temp)
    # return

    # Creamos nuestro objeto LayoutGraph
    layout_gr = LayoutGraph(
        leer_grafo_archivo(args.file_name),
        temp =args.temp,
        iters=args.iters,
        refresh=args.refresh,
        c1=args.c1,
        c2=args.c2,
        ctemp=args.ctemp,
        ancho=args.ancho,
        grav=args.grav,
        verbose=args.verbose
    )

    # Ejecutamos el layout
    layout_gr.layout()
    return


if __name__ == '__main__':
    main()
