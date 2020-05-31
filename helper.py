from matplotlib.patches import Rectangle, Polygon
import numpy as np
EPSILON = 0.0000000002



def show_balk(height, length, ax, fig):
    ax.add_patch(Rectangle((0, -height / 2), length, height))
    fig.show()


def create_mesh(H, L, H_step, L_step):
    def check_key_in_dict(k, d):
        for key in d.keys():
            if (np.array(np.array(k)-np.array(key))[0]**2+np.array(np.array(k)-np.array(key))[1]**2)**0.5<=EPSILON:
                return key
        return None
    H_steps = list(np.arange(0, H, H_step) - H / 2)
    L_steps = list(np.arange(0, L, L_step))
    elements = {}
    nodes = {}
    element_id = 0
    node_id = 0
    for y in H_steps:
        for x in L_steps:
            tmp_rectangle = {
                'ld': (x, y),
                'rd': (x + L_step, y),
                'lt': (x, y + H_step),
                'rt': (x + L_step, y + H_step),
            }
            rectangle = {}
            for key, i in tmp_rectangle.items():
                condition = check_key_in_dict(i, nodes)
                if condition is None:
                # if i not in nodes.keys():
                    nodes[i] = node_id
                    rectangle[key] = node_id
                    node_id += 1
                else:
                    rectangle[key] = nodes[condition]
            triangle1 = rectangle.copy()
            # triangle2 = rectangle.copy()
            triangle1.pop('rt')
            triangle2 = {
                'ld': rectangle['rt'],
                'rd': rectangle['lt'],
                'lt': rectangle['rd'],
            }
            # triangle2.pop('ld')
            elements[element_id] = triangle1
            element_id += 1
            elements[element_id] = triangle2
            element_id += 1

    nodes = {value: key for key, value in nodes.items()}
    return nodes, elements


def show_mesh(elements: dict, nodes: dict, ax):
    for i in elements.values():
        ax.add_patch(Polygon(
            [nodes[x] for x in i.values()],
            facecolor='yellow', edgecolor='violet'
        ))


def solve(nodes, elems, yung, puass, P, H):
    def make_zakrep(nodes, k_glob, f_glob):
        nodes_count = len(nodes)
        for key, val in nodes.items():
            if val[0]==0:
                k_glob[key*2:key*2+2, :] = np.zeros((2, nodes_count*2))
                k_glob[:, key*2:key*2+2] = np.zeros((nodes_count*2, 2))
                k_glob[key*2, key*2] = 1
                k_glob[key*2+1, key*2+1] = 1
                f_glob[key*2:key*2+2, :] = np.zeros((2,1))
                #f_glob[perenumber[key]*2] = 0
                # print(key, val)
        return k_glob, f_glob

    def triangle_square(el):
        def distance(p1, p2):
            return ((p1[0]-p2[0])**2+(p1[1]-p2[1])**2)**0.5
        a = distance(nodes[el['ld']], nodes[el['rd']])
        b = distance(nodes[el['ld']], nodes[el['lt']])
        c = distance(nodes[el['lt']], nodes[el['rd']])
        s = (a + b + c) / 2
        return (s*(s-a)*(s-b)*(s-c)) ** 0.5
    def get_f_loc(elem, P):
        f_loc = np.zeros((6,1))
        x1 = nodes[elem['ld']][0]
        x2 = nodes[elem['rd']][0]
        x3 = nodes[elem['lt']][0]

        y1 = nodes[elem['ld']][1]
        y2 = nodes[elem['rd']][1]
        y3 = nodes[elem['lt']][1]

        edge = []
        for i, point in enumerate([(x1,y1),(x2,y2),(x3,y3)]):
            if abs(point[1]-H/2)<EPSILON:
                # нормаль (nx,ny)=(0,1). Fj = (P*nx,P*ny, P*nx, P*ny, 0, 0)Т*Lj - если 1 и 2 узел треугольника под нагрузкой, а третий нет.
                f_loc[i*2+1] = P
                edge.append(point[0])
        if len(edge) < 2:
            return f_loc
        #print(edge)
        f_loc = f_loc * np.abs(edge[1] - edge[0])
        # print(np.abs(edge[1] - edge[0]))
        # print(f_loc)
        return f_loc
    k_glob = np.zeros((2*len(nodes), 2*len(nodes)))
    f_glob = np.zeros((2*len(nodes), 1))

    # i = 0
    for key, elem in elems.items():
        x1 = nodes[elem['ld']][0]
        x2 = nodes[elem['rd']][0]
        x3 = nodes[elem['lt']][0]

        y1 = nodes[elem['ld']][1]
        y2 = nodes[elem['rd']][1]
        y3 = nodes[elem['lt']][1]
        square = triangle_square(elem)
        B = np.array([
            [y2 - y3, 0, y3 - y1, 0, y1 - y2, 0],
            [0, x3 - x2, 0, x1 - x3, 0, x2 - x1],
            [x3 - x2, y2 - y3, x1 - x3, y3 - y1, x2 - x1, y1 - y2],
        ])/(2*square)
        D = yung*(1-puass)/((1+puass)*(1-2*puass))*np.array([
            [1, puass/(1-puass), 0],
            [puass/(1-puass), 1, 0],
            [0, 0, (1 - 2*puass)/(2*(1-puass))],
        ])
        k_mini = np.transpose(B).dot(D).dot(B) * square

        col1 = elem['ld']
        col2 = elem['rd']
        col3 = elem['lt']
        k_glob[col1*2:col1*2+2, col1*2:col1*2 + 2] += k_mini[0:2, 0:2]
        k_glob[col1*2:col1*2+2, col2*2:col2*2 + 2] += k_mini[0:2, 2:4]
        k_glob[col1*2:col1*2+2, col3*2:col3*2 + 2] += k_mini[0:2, 4:6]

        k_glob[col2*2:col2*2+2, col1*2:col1*2 + 2] += k_mini[2:4, 0:2]
        k_glob[col2*2:col2*2+2, col2*2:col2*2 + 2] += k_mini[2:4, 2:4]
        k_glob[col2*2:col2*2+2, col3*2:col3*2 + 2] += k_mini[2:4, 4:6]

        k_glob[col3*2:col3*2+2, col1*2:col1*2 + 2] += k_mini[4:6, 0:2]
        k_glob[col3*2:col3*2+2, col2*2:col2*2 + 2] += k_mini[4:6, 2:4]
        k_glob[col3*2:col3*2+2, col3*2:col3*2 + 2] += k_mini[4:6, 4:6]
        f_loc = get_f_loc(elem, P)
        f_glob[col1*2:col1*2+2] += f_loc[0:2]
        f_glob[col2*2:col2*2+2] += f_loc[2:4]
        f_glob[col3*2:col3*2+2] += f_loc[4:6]
    k_glob, f_glob = make_zakrep(nodes, k_glob, f_glob)
    # print(f_glob)

    return np.linalg.solve(k_glob, f_glob)


def get_napryazh_and_deform(nodes, elems, U, puass, yung):
    def triangle_square(el):
        def distance(p1, p2):
            return ((p1[0]-p2[0])**2+(p1[1]-p2[1])**2)**0.5
        a = distance(nodes[el['ld']], nodes[el['rd']])
        b = distance(nodes[el['ld']], nodes[el['lt']])
        c = distance(nodes[el['lt']], nodes[el['rd']])
        s = (a + b + c) / 2
        return (s*(s-a)*(s-b)*(s-c)) ** 0.5
    deform_vs_id = {}
    napryazh_vs_id = {}
    for key, elem in elems.items():

        x1 = nodes[elem['ld']][0]
        x2 = nodes[elem['rd']][0]
        x3 = nodes[elem['lt']][0]

        y1 = nodes[elem['ld']][1]
        y2 = nodes[elem['rd']][1]
        y3 = nodes[elem['lt']][1]
        square = triangle_square(elem)

        B = np.array([
            [y2 - y3, 0, y3 - y1, 0, y1 - y2, 0],
            [0, x3 - x2, 0, x1 - x3, 0, x2 - x1],
            [x3 - x2, y2 - y3, x1 - x3, y3 - y1, x2 - x1, y1 - y2],
        ])/(2*square)
        D = yung*(1-puass)/((1+puass)*(1-2*puass))*np.array([
        [1, puass/(1-puass), 0],
        [puass/(1-puass), 1, 0],
        [0, 0, (1 - 2*puass)/(2*(1-puass))],
    ])

        u_elem = np.zeros((6,1))


        u_elem[0] = U[elem['ld']*2][0]
        u_elem[2] = U[elem['rd']*2][0]
        u_elem[4] = U[elem['lt']*2][0]
        u_elem[1] = U[elem['ld']*2 + 1][0]
        u_elem[3] = U[elem['rd']*2 + 1][0]
        u_elem[5] = U[elem['lt']*2 + 1][0]

        deform = B.dot(u_elem)
        napryazh = D.dot(deform)

        napryazh_vs_id[key] = napryazh
        deform_vs_id[key] = deform

    return napryazh_vs_id, deform_vs_id


def show_color_mesh(elements: dict, nodes: dict, ax, color_elements):
    for key, i in elements.items():
        # print(color_elements[key])
        ax.add_patch(Polygon(
            [nodes[x] for x in i.values()],
            facecolor=color_elements[key]
        ))
