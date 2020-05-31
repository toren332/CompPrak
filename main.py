from matplotlib import pyplot as plt
from helper import *
from matplotlib import cm
from collections import OrderedDict


def heatMapColorforValue(value, maximum, minimum):
    h = (value - minimum) / (maximum - minimum)
    return cm.get_cmap('hot')(h)


def get_maximum_minimum(napryazh, kind):
    q = []
    for i in list(napryazh.values()):
        q.append(i[kind][0])
    return max(q), min(q)


def create_colored_elements(napr, kind):
    maximum, minimum = get_maximum_minimum(napr, kind)
    colored_elems = {

    }
    for key, value in napr.items():
        colored_elems[key] = heatMapColorforValue(value[kind][0], maximum, minimum)
    return colored_elems


def get_napr_from_elems(elems, napr, kind):
    new_napr = {}
    for key in elems:
        new_napr[key]=napr[key][kind][0]
    return new_napr


H = 2
L = 10
puass = 0.3
yung = 180

# all_graphs =
# fig, ax = plt.subplots(1)
# ax.axvline(x=0, c='black')
# ax.set_xlim([-L*0.05, L * 1.05])
# ax.set_ylim([-H, H])

# стриоим балку
# show_balk(H, L, ax, fig)

STEPS = np.arange(1, 20, 1)[:20:1]

plot_fig, axs = plt.subplots(3, figsize=(30, 25))

# axs[0].plot(x, y)
# axs[1].plot(x, -y)

for step_N in STEPS[:]:
    P = -0.5 / step_N
    H_step = H / step_N
    L_step = L / step_N

    # строим узлы и элементы
    nodes, elements = create_mesh(H, L, H_step, L_step)
    fig, ax = plt.subplots(1)
    ax.axvline(x=0, c='black')
    ax.set_xlim([-L * 0.05, L * 1.05])
    ax.set_ylim([-H, H])
    show_mesh(elements, nodes, ax)
    # fig.show()

    # Строим перемещения
    U = solve(nodes, elements, yung, puass, P, H)

    napryazh_vs_id, deform_vs_id = get_napryazh_and_deform(nodes, elements, U, puass, yung)
    fig, ax = plt.subplots(1)
    ax.axvline(x=0, c='black')
    ax.set_xlim([-L * 0.05, L * 1.05])
    ax.set_ylim([-H, H])
    new_nodes = {}

    for i in range(int(len(U) / 2)):
        new_nodes[i] = (nodes[i][0] + U[i * 2][0], nodes[i][1] + U[i * 2 + 1][0])
    color_elements = create_colored_elements(napryazh_vs_id, 0)
    new_color_elements = {}
    null_elements = set()
    for key in color_elements.keys():
        if step_N % 2 != 0:
            for j in range(2 * step_N):
                null_elements.add(step_N * step_N + j - step_N)
        else:
            for j in range(2 * step_N):
                null_elements.add(step_N * step_N + j)

    for key, value in color_elements.items():
        if key not in null_elements:
            new_color_elements[key] = (1.0, 0.0, 0.0, 1.0)
        else:
            new_color_elements[key] = (0.0, 0.0, 0.0, 1.0)


    for q_ in range(3):
    # q_ = 2
        new_napr = get_napr_from_elems(null_elements, napryazh_vs_id, q_)
        # print(new_napr.values())
        x = []
        # print(null_elements)
        for element_id in null_elements:
            x.append(nodes[elements[element_id]['ld']][0])
        # if step_N==11:
        #     print(x)
        z = dict(zip(x, list(new_napr.values())))
        z = OrderedDict(sorted(z.items()))
        axs[q_].plot(list(z.keys()), list(z.values()))

    # print(z)
    # show_color_mesh(elements, new_nodes, ax, color_elements)
    show_color_mesh(elements, new_nodes, ax, new_color_elements)
    fig.show()
    fig.savefig(f'steps/step_{step_N}.png')

plot_fig.show()
plot_fig.savefig(f'napr.png')

