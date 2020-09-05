# Desenvolvido por: Luciano Soares <lpsoares@insper.edu.br>
# Disciplina: Computação Gráfica
# Data: 28 de Agosto de 2020

import argparse     # Para tratar os parâmetros da linha de comando
import x3d          # Faz a leitura do arquivo X3D, gera o grafo de cena e faz traversal
import interface    # Janela de visualização baseada no Matplotlib
import gpu          # Simula os recursos de uma GPU
import math


def polypoint2D(point, color):
    """ Função usada para renderizar Polypoint2D. """
    for i in range(0, len(point), 2):
        # transforma a posicao dos pontos em posicao dos pixeis
        x = math.floor(point[i])
        y = math.floor(point[i+1])
        # transforma as cores de x3d (0,1) para o framebuffer (0,255)
        r = int(255*color[0])
        g = int(255*color[1])
        b = int(255*color[2])
        gpu.GPU.set_pixel(x, y, r, g, b)  # altera um pixel da imagem
    # cuidado com as cores, o X3D especifica de (0,1) e o Framebuffer de (0,255)


def polyline2D(lineSegments, color):
    """ Função usada para renderizar Polyline2D. """
    # Pega as coordenadas do inicio e fim da linha
    u1 = lineSegments[0]
    v1 = lineSegments[1]
    u2 = lineSegments[2]
    v2 = lineSegments[3]
    # transforma as cores de x3d (0,1) para o framebuffer (0,255)
    r = int(255*color[0])
    g = int(255*color[1])
    b = int(255*color[2])
    # mult indica se vamos calcular na vertical ou na horizontal

    v = v1
    u = u1
    # calcula o coeficiente angular da reta
    if(u2-u1 != 0):
        s = (v2-v1)/(u2-u1)
    else:
        s = 10
    if(s <= 1 and s >= -1):
        if(u1 <= u2):
            while u <= u2 + 1:
                # altera um pixel da imagem
                gpu.GPU.set_pixel(int(u), int(v), r, g, b)
                v += s
                u += 1
        else:
            while u >= u2 - 1:
                # altera um pixel da imagem
                gpu.GPU.set_pixel(int(u), int(v), r, g, b)
                v -= s
                u -= 1
    else:
        # calcula o coeficiente angular da reta
        if(v2-v1 != 0):
            s = (u2-u1)/(v2-v1)
        else:
            return
        if(v1 <= v2):
            while v < v2 + 1:
                # altera um pixel da imagem
                gpu.GPU.set_pixel(int(u), int(v), r, g, b)
                u += s
                v += 1
        else:
            while v > v2 - 1:
                # altera um pixel da imagem
                gpu.GPU.set_pixel(int(u), int(v), r, g, b)
                u -= s
                v -= 1


def triangleSet2D(vertices, color):
    """ Função usada para renderizar TriangleSet2D. """
    print(vertices)
    # Pega as coordenadas do inicio e fim da linha
    u1 = [vertices[0], vertices[1]]
    u2 = [vertices[2], vertices[3]]
    u3 = [vertices[4], vertices[5]]
    # transforma as cores de x3d (0,1) para o framebuffer (0,255)
    r = int(255*color[0])
    g = int(255*color[1])
    b = int(255*color[2])
    # calcula o vetor transversal
    T1 = [u1[0] - u2[0], u1[1] - u2[1]]
    T2 = [u2[0] - u3[0], u2[1] - u3[1]]
    T3 = [u3[0] - u1[0], u3[1] - u1[1]]
    # calcula o vetor perpendicular (normal)
    N1 = [T1[1], -T1[0]]
    N2 = [T2[1], -T2[0]]
    N3 = [T3[1], -T3[0]]

    i = 0.25
    j = 0.25
    while i < 30:
        while j < 20:
            count = 0
            # calcula o vetor do ponto
            V0 = [i - u1[0], j - u1[1]]
            V1 = [i+0.25 - u1[0], j - u1[1]]
            V2 = [i - u1[0], j+0.25 - u1[1]]
            V3 = [i+0.25 - u1[0], j+0.25 - u1[1]]
            # calcula o produto escalar
            Pe0 = V0[0]*N1[0] + V0[1]*N1[1]
            Pe1 = V1[0]*N1[0] + V1[1]*N1[1]
            Pe2 = V2[0]*N1[0] + V2[1]*N1[1]
            Pe3 = V3[0]*N1[0] + V3[1]*N1[1]

            if Pe0 < 0:
                # calcula o vetor do ponto
                V0 = [i - u2[0], j - u2[1]]
                # calcula o produto escalar
                Pe0 = V0[0]*N2[0] + V0[1]*N2[1]

                if Pe0 < 0:
                    # calcula o vetor do ponto
                    V0 = [i - u3[0], j - u3[1]]
                    # calcula o produto escalar
                    Pe0 = V0[0]*N3[0] + V0[1]*N3[1]
                    if Pe0 < 0:
                        count += 0.25

            if Pe1 < 0:
                # calcula o vetor do ponto
                V1 = [i+0.25 - u2[0], j - u2[1]]
                # calcula o produto escalar
                Pe1 = V1[0]*N2[0] + V1[1]*N2[1]

                if Pe1 < 0:
                    # calcula o vetor do ponto
                    V1 = [i+0.25 - u3[0], j - u3[1]]
                    # calcula o produto escalar
                    Pe1 = V1[0]*N3[0] + V1[1]*N3[1]
                    if Pe1 < 0:
                        count += 0.25

            if Pe2 < 0:
                # calcula o vetor do ponto
                V2 = [i - u2[0], j+0.25 - u2[1]]
                # calcula o produto escalar
                Pe2 = V2[0]*N2[0] + V2[1]*N2[1]

                if Pe2 < 0:
                    # calcula o vetor do ponto
                    V2 = [i - u3[0], j+0.25 - u3[1]]
                    # calcula o produto escalar
                    Pe2 = V2[0]*N3[0] + V2[1]*N3[1]
                    if Pe2 < 0:
                        count += 0.25

            if Pe3 < 0:
                # calcula o vetor do ponto
                V3 = [i+0.25 - u2[0], j+0.25 - u2[1]]
                # calcula o produto escalar
                Pe3 = V3[0]*N2[0] + V3[1]*N2[1]

                if Pe3 < 0:
                    # calcula o vetor do ponto
                    V3 = [i+0.25 - u3[0], j+0.25 - u3[1]]
                    # calcula o produto escalar
                    Pe3 = V3[0]*N3[0] + V3[1]*N3[1]
                    if Pe3 < 0:
                        count += 0.25

            if count > 0:
                gpu.GPU.set_pixel(math.floor(i), math.floor(j), r*count, g*count,
                                  b*count)  # altera um pixel da imagem
            j += 1
        i += 1
        j = 0.5


def triangleSet(point, color):
    """ Função usada para renderizar TriangleSet. """
    print("TriangleSet : pontos = {0}".format(point)) # imprime no terminal pontos

def viewpoint(position, orientation, fieldOfView):
    """ Função usada para renderizar (na verdade coletar os dados) de Viewpoint. """
    print("Viewpoint : position = {0}, orientation = {0}, fieldOfView = {0}".format(position, orientation, fieldOfView)) # imprime no terminal

def transform(translation, scale, rotation):
    """ Função usada para renderizar (na verdade coletar os dados) de Transform. """
    print("Transform : ", end = '')
    if translation:
        print("translation = {0} ".format(translation), end = '') # imprime no terminal
    if scale:
        print("scale = {0} ".format(scale), end = '') # imprime no terminal
    if rotation:
        print("rotation = {0} ".format(rotation), end = '') # imprime no terminal
    print("")

def triangleStripSet(point, stripCount, color):
    """ Função usada para renderizar TriangleStripSet. """
    print("TriangleStripSet : pontos = {0} ".format(point), end = '') # imprime no terminal pontos
    for i, strip in enumerate(stripCount):
        print("strip[{0}] = {1} ".format(i, strip), end = '') # imprime no terminal
    print("")

def indexedTriangleStripSet(point, index, color):
    """ Função usada para renderizar IndexedTriangleStripSet. """
    print("IndexedTriangleStripSet : pontos = {0}, index = {1}".format(point, index)) # imprime no terminal pontos

def box(size, color):
    """ Função usada para renderizar Boxes. """
    print("Box : size = {0}".format(size)) # imprime no terminal pontos


LARGURA = 30
ALTURA = 20

if __name__ == '__main__':

    # Valores padrão da aplicação
    width = LARGURA
    height = ALTURA
    x3d_file = "exemplo4.x3d"
    image_file = "tela.png"

    # Tratando entrada de parâmetro
    # parser para linha de comando
    parser = argparse.ArgumentParser(add_help=False)
    parser.add_argument("-i", "--input", help="arquivo X3D de entrada")
    parser.add_argument("-o", "--output", help="arquivo 2D de saída (imagem)")
    parser.add_argument("-w", "--width", help="resolução horizonta", type=int)
    parser.add_argument("-h", "--height", help="resolução vertical", type=int)
    parser.add_argument(
        "-q", "--quiet", help="não exibe janela de visualização", action='store_true')
    args = parser.parse_args()  # parse the arguments
    if args.input:
        x3d_file = args.input
    if args.output:
        image_file = args.output
    if args.width:
        width = args.width
    if args.height:
        height = args.height

    # Iniciando simulação de GPU
    gpu.GPU(width, height, image_file)

    # Abre arquivo X3D
    scene = x3d.X3D(x3d_file)
    scene.set_resolution(width, height)

    # funções que irão fazer o rendering
    x3d.X3D.render["Polypoint2D"] = polypoint2D
    x3d.X3D.render["Polyline2D"] = polyline2D
    x3d.X3D.render["TriangleSet2D"] = triangleSet2D
    x3d.X3D.render["TriangleSet"] = triangleSet
    x3d.X3D.render["Viewpoint"] = viewpoint
    x3d.X3D.render["Transform"] = transform
    x3d.X3D.render["TriangleStripSet"] = triangleStripSet
    x3d.X3D.render["IndexedTriangleStripSet"] = indexedTriangleStripSet
    x3d.X3D.render["Box"] = box

    # Se no modo silencioso não configurar janela de visualização
    if not args.quiet:
        window = interface.Interface(width, height)
        scene.set_preview(window)

    scene.parse()  # faz o traversal no grafo de cena

    # Se no modo silencioso salvar imagem e não mostrar janela de visualização
    if args.quiet:
        gpu.GPU.save_image()  # Salva imagem em arquivo
    else:
        window.image_saver = gpu.GPU.save_image  # pasa a função para salvar imagens
        window.preview(gpu.GPU._frame_buffer)  # mostra janela de visualização
