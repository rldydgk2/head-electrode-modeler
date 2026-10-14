"""Experimental electrode editor adapted from the March 2021 research script.

Geometry routines retain the original dictionary-based implementation. This
module is import-safe; call run() once per process to initialize the editor.
See docs/algorithm.md for limitations and the original source revision.
"""

from pathlib import Path
import gzip


import numpy as np
from matplotlib import pyplot as plt
import threading
from scipy.linalg import norm
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
from tqdm import tqdm


def open_input(path):
    opener = gzip.open if Path(path).suffix == ".gz" else open
    return opener(path, "rt", encoding="utf-8-sig")


def head_surface():
    """Load marked head surfaces and isolate the skin boundary (1005)."""
    mni2 = open_input(INPUT_PATH)
    with mni2:
        lines = mni2.readlines()
    l = []
    for i in lines:
        l2 = i.split()
        l.append(l2)
    data = []
    node = {}
    facet = {}
    for i in range(len(l)):
        if not l[i][0] == "#":
            data.append(l[i])
    node_num = int(data[0][0])
    for i in range(len(data)):
        if 1 <= i <= node_num:
            node[i] = [
                int(data[i][4]),
                float(data[i][1]),
                float(data[i][2]),
                float(data[i][3]),
            ]
    face_num = int(data[node_num + 1][0])
    for i in range(face_num):
        facet[i + 1] = [
            int(data[node_num + 2 * i + 2][2]),
            int(data[node_num + 2 * i + 3][1]),
            int(data[node_num + 2 * i + 3][2]),
            int(data[node_num + 2 * i + 3][3]),
        ]
    region = {}
    skin = {}
    for i in range(len(node)):
        if node[i + 1][0] == 1005:
            skin[i + 1] = np.array(node[i + 1][1:])
    skin_facet = {}
    for i in facet:
        if facet[i][0] == 1005:
            skin_facet[i] = facet[i][1:]
    return (node, facet, region, skin, skin_facet)


def tumor_surface(radius, center_x, center_y, center_z):
    """Sample a spherical tumor layer for display and triangulation."""
    theta = np.linspace(0, 2 * np.pi, div_t)
    z = np.linspace(center_z - radius, center_z + radius, div_tz)
    radi = []
    z_array = []
    for i in range(len(z)):
        r = (radius**2 - abs(z[i] - center_z) ** 2) ** 0.5
        radi.append(r)
        z_array.append(np.array([z[i]]))
    radi = np.array(radi)
    z_array = np.array(z_array)
    x = np.outer(radi, np.cos(theta)) + center_x
    y = np.outer(radi, np.sin(theta)) + center_y
    ax.plot_surface(x, y, z_array, color="yellow", alpha=0.3)
    return (x, y, z_array)


def add_tumor(xi, yi, zi, boundary_num):
    """Append a tumor layer to the global node and facet tables."""
    many_node = max(node.keys())
    point = []
    for i in range(div_tz):
        if i == 0:
            point.append([2, xi[i][0], yi[i][0], zi[i][0]])
        elif 0 < i < div_tz - 1:
            for j in range(div_t - 1):
                point.append([2, xi[i][j], yi[i][j], zi[i][0]])
        elif i == div_tz - 1:
            point.append([2, xi[i][0], yi[i][0], zi[i][0]])
    for i in range(len(point)):
        node[many_node + i + 1] = point[i]
    facet_list = []
    for i in range(div_tz - 1):
        if i == 0:
            for j in range(div_t - 1):
                if j < div_t - 2:
                    facet_list.append(
                        [
                            boundary_num,
                            many_node + 1,
                            many_node + j + 2,
                            many_node + j + 3,
                        ]
                    )
                else:
                    facet_list.append(
                        [boundary_num, many_node + 1, many_node + 2, many_node + j + 2]
                    )
        elif 0 < i < div_tz - 2:
            for j in range(div_t - 1):
                if j < div_t - 2:
                    facet_list.append(
                        [
                            boundary_num,
                            many_node + (div_t - 1) * (i - 1) + j + 2,
                            many_node + (div_t - 1) * (i - 1) + j + 3,
                            many_node + (div_t - 1) * i + j + 3,
                        ]
                    )
                    facet_list.append(
                        [
                            boundary_num,
                            many_node + (div_t - 1) * i + j + 3,
                            many_node + (div_t - 1) * i + j + 2,
                            many_node + (div_t - 1) * (i - 1) + j + 2,
                        ]
                    )
                else:
                    facet_list.append(
                        [
                            boundary_num,
                            many_node + (div_t - 1) * (i - 1) + j + 2,
                            many_node + (div_t - 1) * (i - 1) + 0 + 2,
                            many_node + (div_t - 1) * i + 0 + 2,
                        ]
                    )
                    facet_list.append(
                        [
                            boundary_num,
                            many_node + (div_t - 1) * i + 0 + 2,
                            many_node + (div_t - 1) * i + j + 2,
                            many_node + (div_t - 1) * (i - 1) + j + 2,
                        ]
                    )
        elif i == div_tz - 2:
            for j in range(div_t - 1):
                if j < div_t - 2:
                    facet_list.append(
                        [
                            boundary_num,
                            many_node + (div_t - 1) * (i - 1) + j + 2,
                            many_node + (div_t - 1) * (i - 1) + j + 3,
                            many_node + (div_t - 1) * i + 0 + 2,
                        ]
                    )
                else:
                    facet_list.append(
                        [
                            boundary_num,
                            many_node + (div_t - 1) * (i - 1) + j + 2,
                            many_node + (div_t - 1) * (i - 1) + 0 + 2,
                            many_node + (div_t - 1) * i + 0 + 2,
                        ]
                    )
    many_facet = max(facet.keys())
    for i in range(len(facet_list)):
        facet[many_facet + i + 1] = facet_list[i]
    region[1000] = [(r_1 + r_2) / 2 + center_x, center_y, center_z, 4]
    region[1001] = [center_x, center_y, center_z, 5]


def facet_centroids():
    """Make mass point of facet."""
    mass_p_dic = {}
    for i in skin_facet:
        mass_p = (
            skin[skin_facet[i][0]] + skin[skin_facet[i][1]] + skin[skin_facet[i][2]]
        ) / 3
        mass_p_dic[i] = mass_p
    return mass_p_dic


def make_neighbor(interest_facet):
    """Find adjacent triangles sharing an edge within a local patch."""
    neighbor = {}
    for i in interest_facet:
        l = {}
        for j in interest_facet:
            if not i == j:
                k = set(interest_facet[i]) - set(interest_facet[j])
                if len(k) == 1:
                    if interest_facet[i][0] == list(k)[0]:
                        l[0] = j
                    if interest_facet[i][1] == list(k)[0]:
                        l[1] = j
                    if interest_facet[i][2] == list(k)[0]:
                        l[2] = j
        p = []
        if len(l) == 3:
            p.append(l[0])
            p.append(l[1])
            p.append(l[2])
            neighbor[i] = p
    return neighbor


def triangle(skin, skin_facet):
    """Triangle."""
    element = []
    for i in skin_facet:
        element.append(
            [skin[skin_facet[i][0]], skin[skin_facet[i][1]], skin[skin_facet[i][2]]]
        )
    ax.add_collection3d(
        Poly3DCollection(
            element, facecolors="cyan", linewidths=0.1, edgecolors="y", alpha=0.1
        )
    )
    plt.show()


def parse_pick_coordinates(p):
    """Parse Matplotlib coordinates and snap to the closest skin node."""
    if p.count("x") >= 1:
        s = p.split(",")
        b = []
        for i in range(len(s)):
            b.append(s[i].strip())
        c = []
        for i in range(len(b)):
            c.append(b[i][2:])
        print(c)
        if c[0].find("−") == 0:
            cx = float("-" + c[0][1:])
        else:
            cx = float(c[0])
        if c[1].find("−") == 0:
            cy = float("-" + c[1][1:])
        else:
            cy = float(c[1])
        if c[2].find("−") == 0:
            cz = float("-" + c[2][1:])
        else:
            cz = float(c[2])
        i = nearest_node(cx, cy, cz, skin)
        kx = skin[i][0]
        ky = skin[i][1]
        kz = skin[i][2]
        return (kx, ky, kz, i)
    else:
        s = p.split(",")
        b = []
        for i in range(len(s)):
            b.append(s[i].rstrip("ged "))
        d = []
        d.insert(0, b[0][8:])
        d.insert(1, b[1][11:])
        print(d)


def nearest_node(x, y, z, interesting_node):
    """Return the closest node ID in a selected set of vertices."""
    d = []
    indexx = []
    for i in interesting_node:
        c = (
            (x - interesting_node[i][0]) ** 2
            + (y - interesting_node[i][1]) ** 2
            + (z - interesting_node[i][2]) ** 2
        )
        d.append(c)
        indexx.append(i)
    e = d.index(min(d))
    return indexx[e]


def onclick(event):
    """Onclick."""
    time_interval = 0.25
    global click_timer
    if event.button == 3:
        p = ax.format_coord(event.xdata, event.ydata)
        kx, ky, kz, i = parse_pick_coordinates(p)
        print(p)
        if click_timer is None:
            click_timer = threading.Timer(
                time_interval, on_singleclick, [event, kx, ky, kz, d_1, d_2, a]
            )
            click_timer.start()
        if event.dblclick:
            click_timer.cancel()
            ax.scatter(kx, ky, kz, color="green")
            on_dblclick(event, switch, i, s, d_1, d_2, a)


def on_dblclick(event, switch, i, s, d_1, d_2, a):
    """On dblclick."""
    global click_timer
    print("You double-clicked", event.button, event.xdata, event.ydata)
    click_timer = None
    global click_num
    edge_xy, edge_z = array_spacing(click_num)
    if switch == 1:
        place_electrode(0, i, s, d_1, d_2, a)
    elif switch == "all":
        place_paired_arrays(i, edge_xy, edge_z)
    num = unused_node_ids(facet)
    export_poly(node, facet, num)
    export_contact_nodes(dic)
    click_num = click_num + 1
    print(f"Model: {len(node):,} nodes, {len(facet):,} facets")
    print("click_num=%d" % click_num)
    return click_num


def on_singleclick(event, x, y, z, d_1, d_2, a):
    """On singleclick."""
    global click_timer
    print("You single-clicked", event.button, event.xdata, event.ydata)
    click_timer = None
    pass


def surface_frame(i):
    """Estimate a surface normal and a tangent from incident skin facets."""
    for j in skin_facet:
        if i == skin_facet[j][0]:
            break
        if i == skin_facet[j][1]:
            break
        if i == skin_facet[j][2]:
            break
    facet_number = j
    p1 = skin[skin_facet[facet_number][0]]
    p2 = skin[skin_facet[facet_number][1]]
    p3 = skin[skin_facet[facet_number][2]]
    v1 = p2 - p1
    v2 = p3 - p1
    n_v = np.cross(v1, v2)
    mag = norm(n_v)
    n_v = n_v / mag
    mag2 = norm(v1)
    v1 = v1 / mag2
    return (n_v, v1)


def build_electrode_stack(i, s, d_1, d_2, r):
    """Sample a layered disk in the local surface tangent frame."""
    x = skin[i][0]
    y = skin[i][1]
    z = skin[i][2]
    t = np.array([0, d_1, d_1 + d_2])
    theta = np.linspace(0, 2 * np.pi, div_c)
    radi = np.linspace(0, r, 2)
    v, n1 = surface_frame(i)
    n2 = np.cross(v, n1)
    mag2 = norm(n2)
    n2 = n2 / mag2
    p0 = np.array([x, y, z]) + s * v
    radi, theta1 = np.meshgrid(radi, theta)
    t, theta2 = np.meshgrid(t, theta)
    X1, Y1, Z1 = [
        p0[i] + v[i] * t + r * np.sin(theta2) * n1[i] + r * np.cos(theta2) * n2[i]
        for i in [0, 1, 2]
    ]
    X2, Y2, Z2 = [
        p0[i] + radi[i] * np.sin(theta1) * n1[i] + radi[i] * np.cos(theta1) * n2[i]
        for i in [0, 1, 2]
    ]
    X3, Y3, Z3 = [
        p0[i]
        + v[i] * d_1
        + radi[i] * np.sin(theta1) * n1[i]
        + radi[i] * np.cos(theta1) * n2[i]
        for i in [0, 1, 2]
    ]
    X4, Y4, Z4 = [
        p0[i]
        + v[i] * (d_1 + d_2)
        + radi[i] * np.sin(theta1) * n1[i]
        + radi[i] * np.cos(theta1) * n2[i]
        for i in [0, 1, 2]
    ]
    ax.plot_surface(X1, Y1, Z1, color="red", alpha=0.5)
    ax.plot_surface(X2, Y2, Z2, color="red", alpha=0.5)
    ax.plot_surface(X3, Y3, Z3, color="red", alpha=0.5)
    ax.plot_surface(X4, Y4, Z4, color="red", alpha=0.5)
    return (X1, Y1, Z1, X2, Y2, Z2, X3, Y3, Z3, X4, Y4, Z4)


def place_electrode(j, i, s, d_1, d_2, r):
    """Build an electrode stack and stitch its hydrogel to the skin."""
    interest_facet = nearby_facets(i, detec_l=15)
    x = skin[i][0]
    y = skin[i][1]
    z = skin[i][2]
    x1, y1, z1, x2, y2, z2, x3, y3, z3, x4, y4, z4 = build_electrode_stack(
        i, s, d_1, d_2, r
    )
    x_list = ring_coordinates(x2)
    y_list = ring_coordinates(y2)
    z_list = ring_coordinates(z2)
    v, n1 = surface_frame(i)
    p = []
    for i in range(len(x_list)):
        p.append(np.array([x_list[i], y_list[i], z_list[i]]))
    p = np.array(p)
    new_node_in_circle = append_electrode_geometry(
        x1, y1, z1, x2, y2, z2, x3, y3, z3, x4, y4, z4
    )
    append_electrode_regions(j, x, y, z, v)
    neighbor = make_neighbor(interest_facet)
    many_dic = len(dic)
    del_facet = make_hydrogel(interest_facet, neighbor, new_node_in_circle, v, many_dic)
    delete_facet(del_facet)


def ring_coordinates(a):
    """Pick p."""
    b = []
    for i in range(len(a)):
        b.append(a[i][1])
    b = np.array(b)
    return b


def append_electrode_geometry(x1, y1, z1, x2, y2, z2, x3, y3, z3, x4, y4, z4):
    """Append the sampled electrode stack to the global mesh tables."""
    point = []
    for i in range(div_c - 1):
        point.append([x1[i][0], y1[i][0], z1[i][0]])
    for i in range(div_c - 1):
        point.append([x1[i][1], y1[i][1], z1[i][1]])
    for i in range(div_c - 1):
        point.append([x1[i][2], y1[i][2], z1[i][2]])
    point.append([x2[0][0], y2[0][0], z2[0][0]])
    point.append([x3[0][0], y3[0][0], z3[0][0]])
    point.append([x4[0][0], y4[0][0], z4[0][0]])
    many = max(node.keys())
    for i in range(len(point)):
        node[many + i + 1] = [1] + point[i]
    many2 = max(facet.keys())
    facet_list = []
    for i in range(2):
        for j in range(div_c - 1):
            if j < div_c - 2:
                facet_list.append(
                    [
                        -1,
                        many + i * (div_c - 1) + j + 1,
                        many + i * (div_c - 1) + j + 2,
                        many + 3 * (div_c - 1) + i + 1,
                    ]
                )
            else:
                facet_list.append(
                    [
                        -1,
                        many + i * (div_c - 1) + j + 1,
                        many + i * (div_c - 1) + 0 + 1,
                        many + 3 * (div_c - 1) + i + 1,
                    ]
                )
    for j in range(div_c - 1):
        if j < div_c - 2:
            facet_list.append(
                [
                    click_num + 1,
                    many + 2 * (div_c - 1) + j + 1,
                    many + 2 * (div_c - 1) + j + 2,
                    many + 3 * (div_c - 1) + 2 + 1,
                ]
            )
        else:
            facet_list.append(
                [
                    click_num + 1,
                    many + 2 * (div_c - 1) + j + 1,
                    many + 2 * (div_c - 1) + 0 + 1,
                    many + 3 * (div_c - 1) + 2 + 1,
                ]
            )
    for i in range(2):
        for j in range(div_c - 1):
            if j < div_c - 2:
                facet_list.append(
                    [
                        -1,
                        many + i * (div_c - 1) + j + 1,
                        many + i * (div_c - 1) + j + 2,
                        many + (i + 1) * (div_c - 1) + j + 2,
                    ]
                )
                facet_list.append(
                    [
                        -1,
                        many + (i + 1) * (div_c - 1) + j + 2,
                        many + (i + 1) * (div_c - 1) + j + 1,
                        many + i * (div_c - 1) + j + 1,
                    ]
                )
            else:
                facet_list.append(
                    [
                        -1,
                        many + i * (div_c - 1) + j + 1,
                        many + i * (div_c - 1) + 0 + 1,
                        many + (i + 1) * (div_c - 1) + 0 + 1,
                    ]
                )
                facet_list.append(
                    [
                        -1,
                        many + (i + 1) * (div_c - 1) + 0 + 1,
                        many + (i + 1) * (div_c - 1) + j + 1,
                        many + i * (div_c - 1) + j + 1,
                    ]
                )
    for i in range(len(facet_list)):
        facet[many2 + i + 1] = facet_list[i]
    new_node_in_circle = range(many + 1, many + div_c)
    return new_node_in_circle


def append_electrode_regions(j, x, y, z, v):
    """Add hydrogel, ceramic, and electrode region seed points."""
    many_region = len(region)
    point = np.array([x, y, z])
    point1 = point + s / 2 * v
    point2 = point + (s + d_1 / 2) * v
    point3 = point + (s + d_1 + d_2 / 2) * v
    region[many_region + 1] = [point1[0], point1[1], point1[2], 1]
    region[many_region + 2] = [point2[0], point2[1], point2[2], 2]
    region[many_region + 3] = [point3[0], point3[1], point3[2], 3]


def near_point(p, interest_skin):
    """Near point."""
    many = max(node.keys())
    many2 = max(facet.keys())
    point_num = []
    for i in range(len(p) - 1):
        e = nearest_node(p[i][0], p[i][1], p[i][2], interest_skin)
        point_num.append(e)
    for i in range(len(point_num)):
        if i < len(point_num) - 1:
            facet[many2 + 2 * i + 1] = [
                -1,
                many - 3 * div_c + i + 1,
                point_num[i],
                many - 3 * div_c + i + 2,
            ]
            facet[many2 + 2 * i + 2] = [
                -1,
                many - 3 * div_c + i + 2,
                point_num[i],
                point_num[i + 1],
            ]
        else:
            facet[many2 + 2 * i + 1] = [
                -1,
                many - 3 * div_c + i + 1,
                point_num[i],
                many - 3 * div_c + 0 + 1,
            ]
            facet[many2 + 2 * i + 2] = [
                -1,
                many - 3 * div_c + 0 + 1,
                point_num[i],
                point_num[0],
            ]


def unused_node_ids(facet):
    """Identify node IDs absent from the current facet table."""
    many_list = []
    for i in facet:
        many_list.append(i)
    prb = []
    for i in facet:
        info = []
        info.append(facet[i][1])
        info.append(facet[i][2])
        info.append(facet[i][3])
        info2 = set(info)
        info2 = list(info2)
        if not len(info) == len(info2):
            prb.append(i)
    a_sub_b = [x for x in many_list if x not in prb]
    return a_sub_b


def export_poly(node, face, num):
    """Write the current piecewise linear complex for external TetGen meshing."""
    f = open(OUTPUT_DIR / "head-electrodes.poly", "w")
    f.write("# Part 1 - node list\n")
    f.write("# node count, 3 dim, no attribute, no boundary marker\n")
    f.write("{} 3 0 0\n".format(len(node)))
    f.write("# Node index, node coordinates\n")
    for i in range(len(node)):
        f.write(
            "{} {} {} {}\n".format(
                i + 1, node[i + 1][1], node[i + 1][2], node[i + 1][3]
            )
        )
    f.write("# Part 2 - facet list\n")
    f.write("# facet count, have boundary marker\n")
    f.write("{} 1\n".format(len(facet)))
    f.write("# facets\n")
    for i in facet:
        if len(facet[i]) == 4:
            f.write(
                "1 0 {}\n3 {} {} {}\n".format(
                    facet[i][0], facet[i][1], facet[i][2], facet[i][3]
                )
            )
        elif len(facet[i]) == 5:
            f.write(
                "1 0 {}\n4 {} {} {} {}\n".format(
                    facet[i][0], facet[i][1], facet[i][2], facet[i][3], facet[i][4]
                )
            )
    f.write("# Part 3 - hole list\n")
    f.write("0 # no hole\n")
    f.write("# Part 4 - region list\n")
    f.write("{} # number of region\n".format(len(region)))
    for i in region:
        f.write(
            "{} {} {} {} {}\n".format(
                i, region[i][0], region[i][1], region[i][2], region[i][3]
            )
        )
    f.close()


def array_node_ids(i, edge_xy, edge_z):
    """Find nine skin nodes for the historical electrode layout."""
    n_v, v1 = surface_frame(i)
    v_1 = np.cross(np.array([0, 0, 1]), n_v)
    if norm(v_1) == 0:
        v_1 = np.cross(np.array([1, 0, 0]), n_v)
    v_1 = v_1 / norm(v_1)
    v_2 = np.cross(n_v, v_1)
    v_2 = v_2 / norm(v_2)
    p = []
    p.append(skin[i] + edge_xy * v_1)
    p.append(skin[i] + edge_xy * v_1 + edge_z * v_2)
    p.append(skin[i] + edge_z * v_2)
    p.append(skin[i] + edge_xy * -v_1 + edge_z * v_2)
    p.append(skin[i] + edge_xy * -v_1)
    p.append(skin[i] + edge_xy * -v_1 - edge_z * v_2)
    p.append(skin[i] + edge_z * -v_2)
    p.append(skin[i] + edge_xy * v_1 - edge_z * v_2)
    result = []
    result.append(i)
    interest_skin2 = nearby_skin_nodes(i, 100)
    for i in p:
        result.append(nearest_node(i[0], i[1], i[2], interest_skin2))
    return result


def opposite_surface_node(i, n_v):
    """Locate a candidate node on the opposite side of the head."""
    point_vec = {}
    for j in skin:
        v = skin[j] - skin[i]
        if not norm(v) == 0:
            v = v / norm(v)
            point_vec[j] = v
        else:
            point_vec[j] = v
    dot = {}
    for j in point_vec:
        if not norm(point_vec[j]) == 0:
            dot[j] = abs(np.dot(n_v, point_vec[j]))
    k = list(dot.values()).index(max(dot.values()))
    i2 = list(dot.keys())[k]
    i3 = nearest_node(skin[i2][0], skin[i2][1], skin[i][2], skin)
    return i3


def place_paired_arrays(i, edge_xy, edge_z):
    """Place the original paired electrode arrays around a selected node."""
    result1 = array_node_ids(i, edge_xy, edge_z)
    n_v, v1 = surface_frame(i)
    i3 = opposite_surface_node(i, n_v)
    result2 = array_node_ids(i3, edge_xy, edge_z)
    result = result1 + result2
    global click_num
    for j in tqdm(range(len(result)), desc="plate", mininterval=1):
        if j <= 8:
            place_electrode(j + 9 * click_num, result[j], s, d_1, d_2, a)
        elif j == 9:
            click_num = click_num + 1
            place_electrode(j + 9 * (click_num - 1), result[j], s, d_1, d_2, a)
        else:
            place_electrode(j + 9 * (click_num - 1), result[j], s, d_1, d_2, a)
    return click_num


def find_facet(interest_facet, i, n_v):
    """Find facet."""
    for j in interest_facet:
        v1 = skin[interest_facet[j][0]]
        v2 = skin[interest_facet[j][1]]
        v3 = skin[interest_facet[j][2]]
        x1 = v2 - v1
        x2 = v3 - v2
        x3 = v1 - v3
        v = np.array(node[i][1:])
        k1 = np.dot(n_v, np.cross(x1, v - v1))
        k2 = np.dot(n_v, np.cross(x2, v - v2))
        k3 = np.dot(n_v, np.cross(x3, v - v3))
        if k1 > 0:
            k1 = 1
        else:
            k1 = 0
        if k2 > 0:
            k2 = 1
        else:
            k2 = 0
        if k3 > 0:
            k3 = 1
        else:
            k3 = 0
        if k1 == k2 == k3 == 1:
            k = j
            break
    return k


def cross_point(i, j, n_v):
    """Cross point."""
    p0 = np.array(node[i][1:])
    p1 = skin[skin_facet[j][0]]
    p2 = skin[skin_facet[j][1]]
    p3 = skin[skin_facet[j][2]]
    v1 = p2 - p1
    v2 = p3 - p1
    v = np.cross(v1, v2)
    l = np.dot(v, p0 - p1) / np.dot(v, n_v)
    p = p0 - l * n_v
    ax.scatter(p[0], p[1], p[2], color="blue")
    return p


def export_contact_nodes(dic):
    """Write intersection points in TetGen node format."""
    f = open(OUTPUT_DIR / "contact-points.node", "w")
    f.write("{}  3  0  0\n".format(len(dic)))
    for i in range(len(dic)):
        f.write(
            "   {}    {}  {}  {}\n".format(
                i + 1, dic[i + 1][0], dic[i + 1][1], dic[i + 1][2]
            )
        )
    f.close()


def determine_direction(neighbor, i, p, p0, n_v):
    """Determine direction."""
    p1 = skin[skin_facet[i][0]]
    p2 = skin[skin_facet[i][1]]
    p3 = skin[skin_facet[i][2]]
    v1 = p2 - p1
    v2 = p3 - p1
    x1 = p2 - p1
    x2 = p3 - p2
    x3 = p1 - p3
    v = np.cross(v1, v2)
    k1 = np.dot(n_v, np.cross(x1, p - p1))
    k2 = np.dot(n_v, np.cross(x2, p - p2))
    k3 = np.dot(n_v, np.cross(x3, p - p3))
    if k1 >= 0:
        k1 = 1
    else:
        k1 = -1
    if k2 >= 0:
        k2 = 1
    else:
        k2 = -1
    if k3 >= 0:
        k3 = 1
    else:
        k3 = -1
    k = [k1, k2, k3]
    if k == [1, -1, 1]:
        out = neighbor[i][0]
        case = 1
    elif k == [1, 1, -1]:
        out = neighbor[i][1]
        case = 2
    elif k == [-1, 1, 1]:
        out = neighbor[i][2]
        case = 3
    elif k == [1, -1, -1]:
        line = np.dot(n_v, np.cross(p - p0, p3 - p0))
        if line >= 0:
            out = neighbor[i][0]
            case = 1
        else:
            out = neighbor[i][1]
            case = 2
    elif k == [-1, 1, -1]:
        line = np.dot(n_v, np.cross(p - p0, p1 - p0))
        if line >= 0:
            out = neighbor[i][1]
            case = 2
        else:
            out = neighbor[i][2]
            case = 3
    elif k == [-1, -1, 1]:
        line = np.dot(n_v, np.cross(p - p0, p2 - p0))
        if line >= 0:
            out = neighbor[i][2]
            case = 3
        else:
            out = neighbor[i][0]
            case = 1
    elif k == [1, 1, 1] or k == [-1, -1, -1]:
        out = i
        case = 0
    x = cross_point2(i, out, p, p0, v, case, n_v)
    return (out, case, x)


def cross_point2(i, out, p, p0, v, case, n_v):
    """Cross point2."""
    nomal_v = np.cross(p - p0, n_v)
    if case == 1:
        t1 = skin[skin_facet[i][1]]
        t2 = skin[skin_facet[i][2]]
    elif case == 2:
        t1 = skin[skin_facet[i][2]]
        t2 = skin[skin_facet[i][0]]
    elif case == 3:
        t1 = skin[skin_facet[i][0]]
        t2 = skin[skin_facet[i][1]]
    elif case == 0:
        t1 = np.array([0, 0, 1])
        t2 = np.array([0, 0, 0])
    t_v = t1 - t2
    if case == 1 or case == 2 or case == 3:
        l = np.dot(nomal_v, t1 - p0) / np.dot(nomal_v, t_v)
        x = t1 - l * t_v
        if not 0 <= l <= 1:
            print("error")
    elif case == 0:
        l = 0
        x = p
    return x


def make_hydrogel(interest_facet, neighbor, new_node_in_circle, v, many_dic):
    """Trace the electrode footprint through neighboring skin triangles."""
    face = []
    point = []
    for i in range(len(new_node_in_circle)):
        number_facet = find_facet(interest_facet, new_node_in_circle[i], v)
        face.append(number_facet)
        p = cross_point(new_node_in_circle[i], number_facet, v)
        dic[many_dic + 1 + i] = p
        point.append(p)
    change = []
    for u in range(len(face)):
        if u < len(face) - 1:
            if not face[u] == face[u + 1]:
                change.append(u)
        elif not face[u] == face[0]:
            change.append(u)
    face_point = []
    for u in change:
        if u == change[0]:
            p0 = point[u]
            p = point[u + 1]
            out = face[u]
            l = [out, [], point[: u + 1]]
            out1, case, p0 = determine_direction(neighbor, out, p, p0, v)
            l[1].append(case)
            l[2].append(p0)
            face_point.append(l)
            while not out1 == face[u + 1]:
                case = find_neighbor(neighbor, out1, out)
                l2 = [out1, [case], [p0]]
                out2, case, p0 = determine_direction(neighbor, out1, p, p0, v)
                l2[1].append(case)
                l2[2].append(p0)
                face_point.append(l2)
                out = out1
                out1 = out2
                ax.scatter(p0[0], p0[1], p0[2], color="black")
        elif change[0] < u < change[len(change) - 1]:
            case = find_neighbor(neighbor, out1, out)
            l2 = [out1, [case], [p0]]
            p0 = point[u]
            p = point[u + 1]
            out = face[u]
            index = change.index(u)
            l2[2] = l2[2] + point[change[index - 1] + 1 : u + 1]
            out1, case, p0 = determine_direction(neighbor, out, p, p0, v)
            l2[1].append(case)
            l2[2].append(p0)
            face_point.append(l2)
            while not out1 == face[u + 1]:
                case = find_neighbor(neighbor, out1, out)
                l2 = [out1, [case], [p0]]
                out2, case, p0 = determine_direction(neighbor, out1, p, p0, v)
                l2[1].append(case)
                l2[2].append(p0)
                face_point.append(l2)
                out = out1
                out1 = out2
                ax.scatter(p0[0], p0[1], p0[2], color="black")
        elif u == change[len(change) - 1]:
            case = find_neighbor(neighbor, out1, out)
            l2 = [out1, [case], [p0]]
            p0 = point[u]
            p = point[0]
            out = face[u]
            index = change.index(u)
            l2[2] = l2[2] + point[change[index - 1] + 1 : u + 1]
            out1, case, p0 = determine_direction(neighbor, out, p, p0, v)
            l2[1].append(case)
            l2[2].append(p0)
            face_point.append(l2)
            while not out1 == face[0]:
                case = find_neighbor(neighbor, out1, out)
                l2 = [out1, [case], [p0]]
                out2, case, p0 = determine_direction(neighbor, out1, p, p0, v)
                l2[1].append(case)
                l2[2].append(p0)
                face_point.append(l2)
                out = out1
                out1 = out2
                ax.scatter(p0[0], p0[1], p0[2], color="black")
            case = find_neighbor(neighbor, out1, out)
            face_point[0][1] = [case] + face_point[0][1]
            if not u == len(point) - 1:
                face_point[0][2] = [p0] + point[u + 1 :] + face_point[0][2]
                revice_k = len(point[u + 1 :])
            else:
                face_point[0][2] = [p0] + face_point[0][2]
                revice_k = 0
    new_node, del_facet = split_facet(face_point)
    side_face(new_node_in_circle, new_node, revice_k)
    return del_facet


def find_neighbor(neighbor, i, j):
    """Find neighbor."""
    if neighbor[i][0] == j:
        case = 1
    if neighbor[i][1] == j:
        case = 2
    if neighbor[i][2] == j:
        case = 3
    return case


def split_facet(face_point):
    """Split intersected skin facets along the electrode contact boundary."""
    many_node = max(node.keys())
    many_facet = max(facet.keys())
    new_node = []
    i_list = []
    face_point2 = []
    k = 1
    for i in face_point:
        i_list.append(i[0])
        if i == face_point[0]:
            for j in range(len(i[2]) - 1):
                node[many_node + k + j] = [1005] + list(i[2][j + 1])
                if not j == len(i[2]) - 2:
                    new_node.append(many_node + k + j)
            k += len(i[2]) - 1
            face_point2.append([i[0], i[1], list(range(many_node + 1, many_node + k))])
        else:
            for j in range(len(i[2]) - 1):
                node[many_node + k + j] = [1005] + list(i[2][j + 1])
                if not j == len(i[2]) - 2:
                    new_node.append(many_node + k + j)
            face_point2.append(
                [
                    i[0],
                    i[1],
                    list(range(many_node + k - 1, many_node + k + len(i[2]) - 1)),
                ]
            )
            k += len(i[2]) - 1
    last_k = k - 1
    face_point2[0][2] = [many_node + last_k] + face_point2[0][2]
    for j in range(last_k):
        skin[many_node + 1 + j] = np.array(node[many_node + 1 + j][1:])
    prb_i = []
    count = {}
    for i in i_list:
        count[i] = 0
    for i in i_list:
        count[i] += 1
    for i in count:
        if not count[i] == 1:
            prb_i.append(i)
    face_point2_overlab = {}
    for i in prb_i:
        face_point2_overlab[i] = []
    g = 1
    for i in range(len(i_list)):
        if i_list[i] not in prb_i:
            if not face_point2[i][1][0] == face_point2[i][1][1]:
                point = [x for x in [1, 2, 3] if x not in face_point2[i][1]]
                point2 = face_point2[i][1][1]
                for j in range(len(face_point2[i][2]) - 1):
                    if j == 0:
                        facet[many_facet + g] = [
                            1005,
                            face_point2[i][2][j],
                            face_point2[i][2][j + 1],
                            skin_facet[face_point2[i][0]][point[0] - 1],
                        ]
                        facet[many_facet + g + 1] = [
                            1005,
                            face_point2[i][2][j],
                            face_point2[i][2][j + 1],
                            skin_facet[face_point2[i][0]][point2 - 1],
                        ]
                        g += 2
                        v1 = skin[face_point2[i][2][j + 1]] - skin[face_point2[i][2][j]]
                        v1 = v1 / norm(v1)
                        v2 = (
                            skin[skin_facet[face_point2[i][0]][point2 - 1]]
                            - skin[face_point2[i][2][j]]
                        )
                        v2 = v2 / norm(v2)
                        v0 = np.cross(v1, v2)
                    else:
                        facet[many_facet + g] = [
                            1005,
                            face_point2[i][2][j],
                            face_point2[i][2][j + 1],
                            skin_facet[face_point2[i][0]][point[0] - 1],
                        ]
                        v1 = skin[face_point2[i][2][j + 1]] - skin[face_point2[i][2][j]]
                        v1 = v1 / norm(v1)
                        v2 = (
                            skin[skin_facet[face_point2[i][0]][point2 - 1]]
                            - skin[face_point2[i][2][j]]
                        )
                        v2 = v2 / norm(v2)
                        v3 = np.cross(v1, v2)
                        if np.dot(v0, v3) > 0:
                            facet[many_facet + g + 1] = [
                                1005,
                                face_point2[i][2][j],
                                face_point2[i][2][j + 1],
                                skin_facet[face_point2[i][0]][point2 - 1],
                            ]
                            g += 2
                        else:
                            facet[many_facet + g + 1] = [
                                1005,
                                face_point2[i][2][j],
                                skin_facet[face_point2[i][0]][point2 - 1],
                                skin_facet[face_point2[i][0]][face_point2[i][1][0] - 1],
                            ]
                            point2 = face_point2[i][1][0]
                            facet[many_facet + g + 2] = [
                                1005,
                                face_point2[i][2][j],
                                face_point2[i][2][j + 1],
                                skin_facet[face_point2[i][0]][point2 - 1],
                            ]
                            g += 3
                if point2 == face_point2[i][1][1]:
                    facet[many_facet + g] = [
                        1005,
                        face_point2[i][2][-1],
                        skin_facet[face_point2[i][0]][point2 - 1],
                        skin_facet[face_point2[i][0]][face_point2[i][1][0] - 1],
                    ]
                    g += 1
            else:
                point = [x for x in [1, 2, 3] if x not in face_point2[i][1]]
                length = []
                l = [0, -1]
                for j in l:
                    length.append(
                        norm(
                            skin[skin_facet[face_point2[i][0]][point[0] - 1]]
                            - skin[face_point2[i][2][j]]
                        )
                    )
                small = length.index(min(length))
                if small == 0:
                    point2_list = [point[0], face_point2[i][1][0], point[1]]
                else:
                    point2_list = [point[1], face_point2[i][1][0], point[0]]
                count_p = 0
                point2 = point2_list[count_p]
                for j in range(len(face_point2[i][2]) - 1):
                    if j == 0:
                        facet[many_facet + g] = [
                            1005,
                            face_point2[i][2][j],
                            face_point2[i][2][j + 1],
                            face_point2[i][2][-1],
                        ]
                        facet[many_facet + g + 1] = [
                            1005,
                            face_point2[i][2][j],
                            face_point2[i][2][j + 1],
                            skin_facet[face_point2[i][0]][point2 - 1],
                        ]
                        g += 2
                        v1 = skin[face_point2[i][2][j + 1]] - skin[face_point2[i][2][j]]
                        v1 = v1 / norm(v1)
                        v2 = (
                            skin[skin_facet[face_point2[i][0]][point2 - 1]]
                            - skin[face_point2[i][2][j]]
                        )
                        v2 = v2 / norm(v2)
                        v0 = np.cross(v1, v2)
                    elif j == len(face_point2[i][2]) - 2:
                        v1 = skin[face_point2[i][2][j + 1]] - skin[face_point2[i][2][j]]
                        v1 = v1 / norm(v1)
                        v2 = (
                            skin[skin_facet[face_point2[i][0]][point2 - 1]]
                            - skin[face_point2[i][2][j]]
                        )
                        v2 = v2 / norm(v2)
                        v3 = np.cross(v1, v2)
                        if np.dot(v0, v3) > 0:
                            facet[many_facet + g] = [
                                1005,
                                face_point2[i][2][j],
                                face_point2[i][2][j + 1],
                                skin_facet[face_point2[i][0]][point2 - 1],
                            ]
                            g += 1
                        elif point2 == point2_list[1]:
                            facet[many_facet + g] = [
                                1005,
                                face_point2[i][2][j],
                                skin_facet[face_point2[i][0]][point2 - 1],
                                skin_facet[face_point2[i][0]][
                                    point2_list[count_p + 1] - 1
                                ],
                            ]
                            count_p += 1
                            point2 = point2_list[count_p]
                            facet[many_facet + g + 1] = [
                                1005,
                                face_point2[i][2][j],
                                face_point2[i][2][j + 1],
                                skin_facet[face_point2[i][0]][point2 - 1],
                            ]
                            g += 2
                        elif point2 == point2_list[0]:
                            facet[many_facet + g] = [
                                1005,
                                face_point2[i][2][j],
                                skin_facet[face_point2[i][0]][point2 - 1],
                                skin_facet[face_point2[i][0]][
                                    point2_list[count_p + 1] - 1
                                ],
                            ]
                            count_p += 1
                            point2 = point2_list[count_p]
                            facet[many_facet + g + 1] = [
                                1005,
                                face_point2[i][2][j],
                                skin_facet[face_point2[i][0]][point2 - 1],
                                skin_facet[face_point2[i][0]][
                                    point2_list[count_p + 1] - 1
                                ],
                            ]
                            count_p += 1
                            point2 = point2_list[count_p]
                            facet[many_facet + g + 2] = [
                                1005,
                                face_point2[i][2][j],
                                face_point2[i][2][j + 1],
                                skin_facet[face_point2[i][0]][point2 - 1],
                            ]
                            g += 3
                    else:
                        facet[many_facet + g] = [
                            6,
                            face_point2[i][2][j],
                            face_point2[i][2][j + 1],
                            face_point2[i][2][-1],
                        ]
                        v1 = skin[face_point2[i][2][j + 1]] - skin[face_point2[i][2][j]]
                        v1 = v1 / norm(v1)
                        v2 = (
                            skin[skin_facet[face_point2[i][0]][point2 - 1]]
                            - skin[face_point2[i][2][j]]
                        )
                        v2 = v2 / norm(v2)
                        v3 = np.cross(v1, v2)
                        if np.dot(v0, v3) > 0:
                            facet[many_facet + g + 1] = [
                                1005,
                                face_point2[i][2][j],
                                face_point2[i][2][j + 1],
                                skin_facet[face_point2[i][0]][point2 - 1],
                            ]
                            g += 2
                        else:
                            facet[many_facet + g + 1] = [
                                1005,
                                face_point2[i][2][j],
                                skin_facet[face_point2[i][0]][point2 - 1],
                                skin_facet[face_point2[i][0]][
                                    point2_list[count_p + 1] - 1
                                ],
                            ]
                            count_p += 1
                            point2 = point2_list[count_p]
                            facet[many_facet + g + 2] = [
                                1005,
                                face_point2[i][2][j],
                                face_point2[i][2][j + 1],
                                skin_facet[face_point2[i][0]][point2 - 1],
                            ]
                            g += 3
                if point2 == point2_list[1]:
                    facet[many_facet + g] = [
                        1005,
                        face_point2[i][2][-1],
                        skin_facet[face_point2[i][0]][point2 - 1],
                        skin_facet[face_point2[i][0]][point2_list[2] - 1],
                    ]
                    g += 1
        elif len(face_point2_overlab[i_list[i]]) == 0:
            face_point2_overlab[i_list[i]] = face_point2[i][1:]
        else:
            face_point2_overlab[i_list[i]] = (
                face_point2_overlab[i_list[i]] + face_point2[i][1:]
            )
    for i in face_point2_overlab:
        if face_point2_overlab[i][0][1] == face_point2_overlab[i][2][0]:
            if not face_point2_overlab[i][0][0] == face_point2_overlab[i][2][1]:
                point = [x for x in [1, 2, 3] if x not in face_point2_overlab[i][0]]
                for j in range(len(face_point2_overlab[i][1]) - 1):
                    if j < len(face_point2_overlab[i][1]) - 2:
                        facet[many_facet + g] = [
                            1005,
                            face_point2_overlab[i][1][j],
                            face_point2_overlab[i][1][j + 1],
                            skin_facet[i][point[0] - 1],
                        ]
                        facet[many_facet + g + 1] = [
                            1005,
                            face_point2_overlab[i][1][j],
                            face_point2_overlab[i][1][j + 1],
                            face_point2_overlab[i][1][-1],
                        ]
                        g += 2
                    else:
                        facet[many_facet + g] = [
                            1005,
                            face_point2_overlab[i][1][j],
                            face_point2_overlab[i][1][j + 1],
                            skin_facet[i][point[0] - 1],
                        ]
                        g += 1
                point = [x for x in [1, 2, 3] if x not in face_point2_overlab[i][2]]
                for j in range(len(face_point2_overlab[i][3]) - 1):
                    facet[many_facet + g] = [
                        1005,
                        face_point2_overlab[i][3][j],
                        face_point2_overlab[i][3][j + 1],
                        skin_facet[i][point[0] - 1],
                    ]
                    facet[many_facet + g + 1] = [
                        1005,
                        face_point2_overlab[i][3][j],
                        face_point2_overlab[i][3][j + 1],
                        face_point2_overlab[i][1][-1],
                    ]
                    g += 2
                facet[many_facet + g] = [
                    1005,
                    face_point2_overlab[i][1][-1],
                    face_point2_overlab[i][3][-1],
                    skin_facet[i][face_point2_overlab[i][0][1] - 1],
                ]
                facet[many_facet + g + 1] = [
                    1005,
                    face_point2_overlab[i][1][-1],
                    face_point2_overlab[i][1][0],
                    skin_facet[i][face_point2_overlab[i][0][1] - 1],
                ]
                g += 2
            else:
                point = [x for x in [1, 2, 3] if x not in face_point2_overlab[i][0]]
                length = []
                length.append(
                    norm(
                        skin[skin_facet[i][point[0] - 1]]
                        - skin[face_point2_overlab[i][1][-1]]
                    )
                )
                length.append(
                    norm(
                        skin[skin_facet[i][point[0] - 1]]
                        - skin[face_point2_overlab[i][3][0]]
                    )
                )
                small = length.index(min(length))
                if small == 0:
                    pass
                else:
                    face_point2_overlab[i] = [
                        face_point2_overlab[i][2],
                        face_point2_overlab[i][3],
                        face_point2_overlab[i][0],
                        face_point2_overlab[i][1],
                    ]
                for j in range(len(face_point2_overlab[i][1]) - 1):
                    if j < len(face_point2_overlab[i][1]) - 2:
                        facet[many_facet + g] = [
                            1005,
                            face_point2_overlab[i][1][j],
                            face_point2_overlab[i][1][j + 1],
                            skin_facet[i][point[0] - 1],
                        ]
                        facet[many_facet + g + 1] = [
                            1005,
                            face_point2_overlab[i][1][j],
                            face_point2_overlab[i][1][j + 1],
                            face_point2_overlab[i][1][-1],
                        ]
                        g += 2
                    else:
                        facet[many_facet + g] = [
                            1005,
                            face_point2_overlab[i][1][j],
                            face_point2_overlab[i][1][j + 1],
                            skin_facet[i][point[0] - 1],
                        ]
                        g += 1
                for j in range(len(face_point2_overlab[i][3]) - 1):
                    facet[many_facet + g] = [
                        1005,
                        face_point2_overlab[i][3][j],
                        face_point2_overlab[i][3][j + 1],
                        face_point2_overlab[i][1][-1],
                    ]
                    g += 1
                facet[many_facet + g] = [
                    1005,
                    face_point2_overlab[i][1][-1],
                    face_point2_overlab[i][3][-1],
                    face_point2_overlab[i][1][0],
                ]
                g += 1
                point2 = face_point2_overlab[i][2][1]
                for j in range(len(face_point2_overlab[i][3]) - 1):
                    if j == 0:
                        facet[many_facet + g] = [
                            1005,
                            face_point2_overlab[i][3][j],
                            face_point2_overlab[i][3][j + 1],
                            skin_facet[i][point2 - 1],
                        ]
                        g += 1
                        v1 = (
                            skin[face_point2_overlab[i][3][j + 1]]
                            - skin[face_point2_overlab[i][3][j]]
                        )
                        v1 = v1 / norm(v1)
                        v2 = (
                            skin[skin_facet[i][point2 - 1]]
                            - skin[face_point2_overlab[i][3][j]]
                        )
                        v2 = v2 / norm(v2)
                        v0 = np.cross(v1, v2)
                    else:
                        v1 = (
                            skin[face_point2_overlab[i][3][j + 1]]
                            - skin[face_point2_overlab[i][3][j]]
                        )
                        v1 = v1 / norm(v1)
                        v2 = (
                            skin[skin_facet[i][point2 - 1]]
                            - skin[face_point2_overlab[i][3][j]]
                        )
                        v2 = v2 / norm(v2)
                        v3 = np.cross(v1, v2)
                        if np.dot(v0, v3) > 0:
                            facet[many_facet + g] = [
                                1005,
                                face_point2_overlab[i][3][j],
                                face_point2_overlab[i][3][j + 1],
                                skin_facet[i][point2 - 1],
                            ]
                            g += 1
                        else:
                            facet[many_facet + g] = [
                                1005,
                                face_point2_overlab[i][3][j],
                                skin_facet[i][point2 - 1],
                                skin_facet[i][face_point2_overlab[i][2][0] - 1],
                            ]
                            point2 = face_point2_overlab[i][2][0]
                            facet[many_facet + g + 1] = [
                                1005,
                                face_point2_overlab[i][3][j],
                                face_point2_overlab[i][3][j + 1],
                                skin_facet[i][point2 - 1],
                            ]
                            g += 2
                if point2 == face_point2_overlab[i][2][1]:
                    facet[many_facet + g] = [
                        1005,
                        face_point2_overlab[i][3][-1],
                        skin_facet[i][face_point2_overlab[i][2][0] - 1],
                        skin_facet[i][point2 - 1],
                    ]
                    g += 1
        else:
            point = [x for x in [1, 2, 3] if x not in face_point2_overlab[i][0]]
            for j in range(len(face_point2_overlab[i][1]) - 1):
                if j > 0:
                    facet[many_facet + g] = [
                        1005,
                        face_point2_overlab[i][1][j],
                        face_point2_overlab[i][1][j + 1],
                        skin_facet[i][point[0] - 1],
                    ]
                    facet[many_facet + g + 1] = [
                        1005,
                        face_point2_overlab[i][1][j],
                        face_point2_overlab[i][1][j + 1],
                        face_point2_overlab[i][1][0],
                    ]
                    g += 2
                else:
                    facet[many_facet + g] = [
                        1005,
                        face_point2_overlab[i][1][j],
                        face_point2_overlab[i][1][j + 1],
                        skin_facet[i][point[0] - 1],
                    ]
                    g += 1
            point = [x for x in [1, 2, 3] if x not in face_point2_overlab[i][2]]
            for j in range(len(face_point2_overlab[i][3]) - 1):
                facet[many_facet + g] = [
                    1005,
                    face_point2_overlab[i][3][j],
                    face_point2_overlab[i][3][j + 1],
                    skin_facet[i][point[0] - 1],
                ]
                facet[many_facet + g + 1] = [
                    1005,
                    face_point2_overlab[i][3][j],
                    face_point2_overlab[i][3][j + 1],
                    face_point2_overlab[i][1][0],
                ]
                g += 2
            facet[many_facet + g] = [
                1005,
                face_point2_overlab[i][1][0],
                face_point2_overlab[i][3][0],
                skin_facet[i][face_point2_overlab[i][0][0] - 1],
            ]
            facet[many_facet + g + 1] = [
                1005,
                face_point2_overlab[i][1][-1],
                face_point2_overlab[i][1][0],
                skin_facet[i][face_point2_overlab[i][0][0] - 1],
            ]
            g += 2
    for j in range(g - 1):
        skin_facet[many_facet + 1 + j] = facet[many_facet + 1 + j][1:]
        mass_p = (
            skin[skin_facet[many_facet + 1 + j][0]]
            + skin[skin_facet[many_facet + 1 + j][1]]
            + skin[skin_facet[many_facet + 1 + j][2]]
        ) / 3
        mass_p_dic[many_facet + 1 + j] = mass_p
    for j in range(many_facet + 1, many_facet + g):
        orient_facet(j)
    return (new_node, i_list)


def orient_facet(i):
    """Orient a skin facet using its centroid and outward normal."""
    v1 = skin[skin_facet[i][1]] - skin[skin_facet[i][0]]
    v2 = skin[skin_facet[i][2]] - skin[skin_facet[i][0]]
    n_v = np.cross(v1, v2)
    if not np.dot(n_v, mass_p_dic[i]) > 0:
        if len(facet[i]) == 4:
            facet[i] = [facet[i][0], facet[i][3], facet[i][2], facet[i][1]]
            skin_facet[i] = [skin_facet[i][2], skin_facet[i][1], skin_facet[i][0]]
        elif len(facet[i]) == 5:
            facet[i] = [facet[i][0], facet[i][4], facet[i][3], facet[i][2], facet[i][1]]
            skin_facet[i] = [
                skin_facet[i][3],
                skin_facet[i][2],
                skin_facet[i][1],
                skin_facet[i][0],
            ]


def side_face(new_node_in_circle, new_node, revice_k):
    """Side face."""
    many_node = max(node.keys())
    many_facet = max(facet.keys())
    g = 0
    revice_new_node_in_circle = list(new_node_in_circle[-revice_k:]) + list(
        new_node_in_circle[:-revice_k]
    )
    for i in range(len(revice_new_node_in_circle)):
        if i < len(revice_new_node_in_circle) - 1:
            k = 0
            while not new_node[i] + k == new_node[i + 1]:
                facet[many_facet + g + 1] = [
                    -1,
                    revice_new_node_in_circle[i],
                    new_node[i] + k,
                    new_node[i] + k + 1,
                ]
                k += 1
                g += 1
            facet[many_facet + g + 1] = [
                -1,
                revice_new_node_in_circle[i],
                new_node[i + 1],
                revice_new_node_in_circle[i + 1],
            ]
            g += 1
        else:
            k = 0
            while not new_node[i] + k == many_node:
                facet[many_facet + g + 1] = [
                    -1,
                    revice_new_node_in_circle[i],
                    new_node[i] + k,
                    new_node[i] + k + 1,
                ]
                k += 1
                g += 1
            facet[many_facet + g + 1] = [
                -1,
                revice_new_node_in_circle[i],
                new_node[i] + k,
                new_node[0],
            ]
            facet[many_facet + g + 2] = [
                -1,
                revice_new_node_in_circle[i],
                new_node[0],
                revice_new_node_in_circle[0],
            ]


def delete_facet(del_facet_list):
    """Delete facet."""
    del_facet2 = set(del_facet_list)
    for i in list(del_facet2):
        del facet[i]
        del skin_facet[i]
        del mass_p_dic[i]


def nearby_facets(i, detec_l):
    """Select facets with centroids inside a local search radius."""
    point = skin[i]
    interest_facet = {}
    for j in mass_p_dic:
        if norm(point - mass_p_dic[j]) ** 2 <= detec_l**2:
            interest_facet[j] = skin_facet[j]
    return interest_facet


def nearby_skin_nodes(i, detec_l):
    """Select skin vertices inside a local search radius."""
    point = skin[i]
    interest_skin = {}
    for j in skin:
        if norm(point - skin[j]) ** 2 <= detec_l**2:
            interest_skin[j] = skin[j]
    return interest_skin


def array_spacing(click_num):
    """Alternate array spacing for successive pairs of placements."""
    plat_type = 1 + click_num // 2 % 2
    if plat_type == 2:
        edge_xy = 25
        edge_z = 36
    elif plat_type == 1:
        edge_xy = 36
        edge_z = 25
    return (edge_xy, edge_z)


def _initialize():
    global a, ax, center_x, center_y, center_z, cid, click_num, d_1, d_2, dic, div_c, div_t, div_tz, facet, fig, mass_p_dic, new_facet, node, r_1, r_2, region, s, skin, skin_facet, switch, click_timer, xi, xi2, yi, yi2, zi, zi2
    switch = "all"
    # Layer thicknesses and electrode radius, in millimeters.
    s = 3.0
    d_1 = 0.2
    d_2 = 1.0
    a = 9.0
    div_c = 10
    # Sampling and concentric tumor geometry.
    div_t = 30
    div_tz = 30
    r_1 = 5.0
    r_2 = 0.8 * r_1
    center_x, center_y, center_z = (33, -45, 23)
    click_num = 0
    fig = plt.figure()
    ax = fig.add_subplot(111, projection="3d")
    ax.set_xlabel("X")
    ax.set_ylabel("Y")
    ax.set_zlabel("Z")
    plt.xlim(-100, 100)
    plt.ylim(-100, 100)
    ax.set_zlim(-100, 100)
    node, facet, region, skin, skin_facet = head_surface()
    dic = {}
    new_facet = {}
    xi, yi, zi = tumor_surface(r_1, center_x, center_y, center_z)
    xi2, yi2, zi2 = tumor_surface(r_2, center_x, center_y, center_z)
    add_tumor(xi, yi, zi, 1007)
    add_tumor(xi2, yi2, zi2, 1008)
    mass_p_dic = facet_centroids()
    click_timer = None
    cid = fig.canvas.mpl_connect("button_press_event", onclick)
    print(f"Model: {len(node):,} nodes, {len(facet):,} facets")


def run(input_path, output_dir, *, single=False, show=True):
    """Initialize the research editor and optionally enter the GUI event loop."""
    global INPUT_PATH, OUTPUT_DIR, switch
    INPUT_PATH = Path(input_path)
    OUTPUT_DIR = Path(output_dir)
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    _initialize()
    switch = 1 if single else "all"
    faces = list(skin_facet.values())
    stride = max(1, len(faces) // 40000)
    polygons = [[skin[index] for index in face] for face in faces[::stride]]
    ax.add_collection3d(
        Poly3DCollection(
            polygons, facecolors="lightsteelblue", edgecolors="none", alpha=0.4
        )
    )
    ax.set_box_aspect((1, 1, 1))
    ax.set_title("Research editor: right-double-click to place electrodes")
    if show:
        plt.show()
